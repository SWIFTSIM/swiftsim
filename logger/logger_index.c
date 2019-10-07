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
#include "logger_reader.h"
#include "radix_sort.h"

/**
 * @brief Read the index file header.
 *
 * @param index The #logger_index.
 * @param filename The filename of the index file.
 */
void logger_index_read_header(struct logger_index *index,
                              const char *filename) {

  /* Open the file */
  message("Reading %s", filename);
  FILE *fd = fopen(filename, "r");

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
}

/**
 * @brief Write that the file is sorted.
 *
 * WARNING The file must be mapped.
 *
 * @param index The #logger_index.
 */
void logger_index_write_sorted(struct logger_index *index) {
  /* Get the offset */
  size_t offset = sizeof(double) + (1 + swift_type_count) * sizeof(long long);

  /* Set the value */
  char is_sorted = 1;

  /* Write the value */
  memcpy(index->index.map + offset, &is_sorted, sizeof(char));
}
/**
 * @brief Get the #index_data of a particle type.
 *
 * @param index The #logger_index.
 * @param type The particle type.
 */
struct index_data *logger_index_get_data(struct logger_index *index, int type) {
  /* Compute the header size */
  const size_t header = sizeof(double) +
                        (1 + swift_type_count) * sizeof(long long) +
                        sizeof(char);

  /* Count the offset due to the previous types */
  size_t count = 0;
  for (int i = 0; i < type; i++) {
    count += index->nparts[i];
  }
  count *= sizeof(struct index_data);

  return index->index.map + count + header;
}

/**
 * @brief Map the file and if required sort it.
 *
 * @param index The #logger_index.
 * @param filename The index filename.
 * @param sorted Does the file needs to be sorted?
 */
void logger_index_map_file(struct logger_index *index, const char *filename,
                           int sorted) {
  /* Un-map previous file */
  if (index->index.map != NULL) {
    logger_index_free(index);
  }

  /* Read header */
  logger_index_read_header(index, filename);

  /* Check if need to sort the file */
  if (sorted && !index->is_sorted) {
    if (index->reader->verbose > 0) {
      message("Sorting the index file.");
    }
    /* Map the index file */
    logger_loader_io_mmap_file(&index->index, filename,
                               /* read_only */ 0);
    /* Sort the file */
    for (int i = 0; i < swift_type_count; i++) {
      struct index_data *data = logger_index_get_data(index, i);
      radix_sort(data, index->nparts[i]);
    }

    /* Write that the file is sorted */
    logger_index_write_sorted(index);

    /* Free the index file before opening it again in read only */
    logger_index_free(index);

    if (index->reader->verbose > 0) {
      message("Sorting done.");
    }
  }

  /* Map the index file */
  logger_loader_io_mmap_file(&index->index, filename,
                             /* read_only */ 1);
}

/**
 * @brief Cleanup the memory of a logger_index
 *
 * @param index The #logger_index.
 */
void logger_index_free(struct logger_index *index) {
  logger_loader_io_munmap_file(&index->index);
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
size_t logger_index_get_particle(struct logger_index *index, long long id,
                                 int type) {
  /* Define a few variables */
  struct index_data *data = logger_index_get_data(index, type);
  size_t left = 0;
  size_t right = index->nparts[type] - 1;

  /* Search for the value (binary search) */
  while (left <= right) {
    size_t m = (left + right) / 2;
    if (data[m].id < id) {
      left = m + 1;
    } else if (data[m].id > id) {
      right = m - 1;
    } else {
      return data[m].offset;
    }
  }

  return 0;
}

/**
 * @brief Initialize the #logger_index.
 *
 * @param index The #logger_index.
 * @param reader The #logger_reader.
 */
void logger_index_init(struct logger_index *index,
                       struct logger_reader *reader) {
  /* Set the mapped file to NULL */
  index->index.map = NULL;

  /* Set the pointer to the reader */
  index->reader = reader;
}
