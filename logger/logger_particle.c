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

int gravity_logger_mask_id[gravity_logger_field_count];
int stars_logger_mask_id[stars_logger_field_count];
int hydro_logger_mask_id[hydro_logger_field_count];

/**
 * @brief Read an hydro-part record in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param fields_wanted Indices (in #hydro_logger_fields) of the fields to
 * return. They are assumed to be sorted according to #hydro_logger_fields.
 * @param n_fields_wanted Number of fields requested.
 * @param output Array of buffer where the data are written (size given by
 * n_fields).
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 *
 * @return position after the record.
 */
size_t logger_particle_read(const struct logger_reader *reader, size_t offset,
                            const int *fields_wanted, const int n_fields_wanted,
                            void **output, size_t *mask, size_t *h_offset) {
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
  int current_field = 0;
  for (int i = 0; i < hydro_logger_field_count; i++) {
    const int field_index = hydro_logger_mask_id[i];

    /* Is the mask present? */
    if (!(*mask & h->masks[field_index].mask)) {
      /* Requested mask not present, go to the next. */
      if (current_field < n_fields_wanted &&
          fields_wanted[current_field] == i) {
        current_field += 1;
      }
      continue;
    }

    /* Read the data if needed and update the buffer position? */
    if (current_field < n_fields_wanted && fields_wanted[current_field] == i) {
      map = logger_loader_io_read_data(map, h->masks[field_index].size,
                                       output[current_field]);
      current_field += 1;
    } else {
      map += h->masks[field_index].size;
    }
  }
  return map - reader->log.log.map;
}

/**
 * @brief Read a gpart record in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param output Array of buffer where the data are written (size given by
 * gravity_logger_field_count).
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 *
 * @return position after the record.
 */
size_t logger_gparticle_read(const struct logger_reader *reader, size_t offset,
                             void **output, size_t *mask, size_t *h_offset) {
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
  for (int i = 0; i < gravity_logger_field_count; i++) {
    const int field_index = gravity_logger_mask_id[i];

    /* Is the mask present? */
    if (!(*mask & h->masks[field_index].mask)) {
      continue;
    }

    /* Read the data if needed and update the buffer position? */
    map =
        logger_loader_io_read_data(map, h->masks[field_index].size, output[i]);
  }
  return map - reader->log.log.map;
}

/**
 * @brief Read a spart record in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param fields_wanted Indices (in #stars_logger_fields) of the fields to
 * return. They are assumed to be sorted according to #stars_logger_fields.
 * @param n_fields_wanted Number of fields requested.
 * @param output Array of buffer where the data are written (size given by
 * n_fields).
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 *
 * @return position after the record.
 */
size_t logger_sparticle_read(const struct logger_reader *reader, size_t offset,
                             const int *fields_wanted,
                             const int n_fields_wanted, void **output,
                             size_t *mask, size_t *h_offset) {
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
  int current_field = 0;
  for (int i = 0; i < stars_logger_field_count; i++) {
    const int field_index = stars_logger_mask_id[i];

    /* Is the mask present? */
    if (!(*mask & h->masks[field_index].mask)) {
      /* Requested mask not present, go to the next. */
      if (current_field < n_fields_wanted &&
          fields_wanted[current_field] == i) {
        current_field += 1;
      }
      continue;
    }

    /* Read the data if needed and update the buffer position? */
    if (current_field < n_fields_wanted && fields_wanted[current_field] == i) {
      map = logger_loader_io_read_data(map, h->masks[field_index].size,
                                       output[current_field]);
      current_field += 1;
    } else {
      map += h->masks[field_index].size;
    }
  }
  return map - reader->log.log.map;
}
