/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Include the corresponding header */
#include "logger_index.h"

/* Include the standard headers */
#include <errno.h>
#include <fcntl.h>

/* Include local headers */
#include "logger_loader_io.h"

/**
 * @brief Read the index file header.
 *
 * @param index The #logger_index.
 * @param filename The filename of the index file.
 */
void logger_index_read_header(struct logger_index *index, const char *filename) {

  /* Open the file */
  FILE *fd = fopen(filename, O_RDONLY);

  if (fd == NULL)
    error("Unable to open file %s (%s).", filename, strerror(errno));

  /* Read times */
  fread(&index->time, sizeof(double), 1, fd);
  fread(&index->integer_time, sizeof(long long), 1, fd);

  /* Read the number of particles */
  fread(index->nparts, sizeof(long long), swift_type_count, fd);

  /* Read if the file is sorted */
  fread(&index->is_sorted, sizeof(char), 1, fd);

  /* Close the file */
  fclose(fd);

  /* Set the mapped file to NULL */
  index->log.map = NULL;
}

/**
 * @brief Map the file and if required sort it.
 *
 * @param index The #logger_index.
 * @param filename The index filename.
 * @param sorted Does the file needs to be sorted?
 */
void logger_index_map_file(struct logger_index *index, const char *filename, int sorted) {

  /* Map the index file */
  index->log.map = logger_loader_io_mmap_file(filename, &index->log.file_size,
                                            /* read_only */ 1);

  /* Sort the file if required */
  if (sorted && !index->is_sorted) {
    logger_index_sort_file(index);
  }

}
