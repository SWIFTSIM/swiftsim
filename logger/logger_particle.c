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
#include "logger_particle.h"
#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_reader.h"
#include "logger_time.h"
#include "logger_tools.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Read a particle (of any type) record in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param output Array of buffer where the data are written (size given by
 * local_count).
 * @param local_to_global The list converting local to global id for the fields
 * (e.g. gravity_logger_local_to_global)
 * @param local_count The number of element in logger_mask_id (e.g.
 * gravity_logger_field_count)
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 *
 * @return position after the record.
 */
__attribute__((always_inline)) INLINE size_t
logger_particle_read(const struct logger_reader *reader, size_t offset,
                     void **output, const int *local_to_global,
                     const int local_count, size_t *mask, size_t *h_offset) {

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  *mask = 0;
  *h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, (char *)map + offset, mask, h_offset);

  /* Check that the mask is meaningful */
  if (*mask > (unsigned int)(1 << h->masks_count)) {
    error("Found an unexpected mask %zi", *mask);
  }

  /* Check if it is not a time record. */
  if (*mask == h->timestamp_mask) {
    error("Unexpected timestamp while reading a particle: %lu.", *mask);
  }

  /* Read the record and copy it to a particle. */
  for (int local = 0; local < local_count; local++) {
    const int global = local_to_global[local];

    /* Is the mask present? */
    if (!(*mask & h->masks[global].mask)) {
      continue;
    }

    /* Read the data if needed and update the buffer position? */
    map = logger_loader_io_read_data(map, h->masks[global].size, output[local]);
  }
  return map - reader->log.log.map;
}
