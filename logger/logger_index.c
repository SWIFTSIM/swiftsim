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
#include "radix_sort.h"

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
 * @brief Get the #index_data of a particle type.
 *
 * @param index The #logger_index.
 * @param type The particle type.
 */
struct index_data *logger_index_get_data(struct logger_index *index, int type) {
  /* Compute the header size */
  const size_t header = sizeof(double) + 2 * sizeof(long long)
    + sizeof(char);

  /* Count the offset due to the previous types */
  size_t count = 0;
  for(int i = 0; i < type; i++) {
    count += index->nparts[i];
  }
  count *= sizeof(struct index_data);

  return index->log.map + count + header;
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
    for(int i = 0; i < swift_type_count; i++) {
      struct index_data *data = logger_index_get_data(index, i);
      radix_sort(data, index->nparts[i]);
    }
  }

}

/**
 * @brief Cleanup the memory of a logger_index
 *
 * @param index The #logger_index.
 */
void logger_index_free(struct logger_index *index) {
  logger_loader_io_munmap_file(index->log.map, index->log.file_size);

  index->log.map = NULL;
}

/**
 * @brief Get the offset of a given particle
 *
 * @param index The #logger_index.
 * @param id The ID of the particle.
 * @param type The type of the particle.
 *
 * @return The offset of the particle or 0 if not found.
 */
size_t logger_index_get_particle(struct logger_index *index, long long id, int type) {
  /* Define a few variables */
  struct index_data *data = logger_index_get_data(index, type);
  size_t left = 0;
  size_t right = index->nparts[type] - 1;

  /* Search for the value (binary search) */
  while (left <= right) {
    size_t m = (left + right) / 2;
    if (data[m].id < id) {
      left = m + 1;
    }
    else if (data[m].id > id) {
      right = m - 1;
    }
    else {
      return data[m].offset;
    }

  }

  return 0;
}
