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

#ifndef LOGGER_LOGGER_INDEX_H
#define LOGGER_LOGGER_INDEX_H

#include "logger_tools.h"

/* predefine the structure */
struct logger_reader;

/**
 * @brief Data structure contained in the logger files.
 */
struct index_data {
  /* Id of the particle. */
  long long id;

  /* Offset of the particle in the file. */
  size_t offset;
};

/**
 * @brief Structure dealing with the index files.
 *
 * The structure is initialized with #logger_index_init and
 * then a file can be read with #logger_index_read_header and
 * #logger_index_map_file.
 *
 * The functions #logger_index_get_particle_offset and
 * #logger_index_get_data should be used to access the element
 * stored inside the index file.
 * The first one access a particle through its ids and the second one
 * gives a pointer to the first element that can be looped through.
 */
struct logger_index {
  /* Pointer to the reader */
  struct logger_reader *reader;

  /* Time of the index file */
  double time;

  /* Integer time of the index file */
  long long integer_time;

  /* Number of particles in the file */
  long long nparts[swift_type_count];

  /* Is the file sorted ? */
  char is_sorted;

  /* The index file's variables. */
  struct {
    /* Mapped data. */
    void *map;

    /* File size. */
    size_t file_size;

  } log;
};

void logger_index_write_sorted(struct logger_index *index);
void logger_index_init(struct logger_index *index,
                       struct logger_reader *reader);
void logger_index_read_header(struct logger_index *index, const char *filename);
void logger_index_map_file(struct logger_index *index, const char *filename,
                           int sorted);
size_t logger_index_get_particle_offset(struct logger_index *index,
                                        long long id, int type);
void logger_index_free(struct logger_index *index);
void logger_index_sort_file(struct logger_index *index);
struct index_data *logger_index_get_data(struct logger_index *index, int type);

#endif  // LOGGER_LOGGER_INDEX_H
