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

/**
 * @brief Read a particle entry in the log file.
 *
 * @param reader The #logger_reader.
 * @param part The #logger_particle to update.
 * @param offset offset of the record to read.
 * @param time time to interpolate (not used if constant interpolation).
 * @param reader_type #logger_reader_type.
 *
 * @return position after the record.
 */
size_t logger_particle_read(struct logger_particle *part,
                            const struct logger_reader *reader, size_t offset,
                            const double time,
                            const enum logger_reader_type reader_type) {

  /* /\* Save the offset *\/ */
  /* part->offset = offset; */

  /* /\* Get a few pointers. *\/ */
  /* const struct header *h = &reader->log.header; */
  /* void *map = reader->log.log.map; */

  /* const struct time_array *times = &reader->log.times; */

  /* size_t mask = 0; */
  /* size_t h_offset = 0; */

  /* logger_particle_init(part); */

  /* /\* Read the record's mask. *\/ */
  /* map = logger_loader_io_read_mask(h, (char *)map + offset, &mask, &h_offset); */

  /* /\* Check that the mask is meaningful *\/ */
  /* if (mask > (unsigned int)(1 << h->masks_count)) { */
  /*   error("Found an unexpected mask %zi", mask); */
  /* } */

  /* /\* Check if it is not a time record. *\/ */
  /* if (mask == h->timestamp_mask) { */
  /*   error("Unexpected timestamp while reading a particle: %lu.", mask); */
  /* } */

  /* /\* Compute the size of the record. *\/ */
  /* size_t record_size = 0; */
  /* for (int i = 0; i < h->masks_count; i++) { */
  /*   if (mask & h->masks[i].mask) { */
  /*     record_size += h->masks[i].size; */
  /*   } */
  /* } */
  /* /\* Read the record and copy it to a particle *\/ */
  /* char *buff = (char *)malloc(record_size); */
  /* map = logger_loader_io_read_data(map, record_size, buff); */
  /* for (int i = 0; i < h->masks_count; i++) { */
  /*   if (mask & h->masks[i].mask) { */
  /*     logger_particle_read_single_field(part, buff, h->masks[i].name, */
  /*                                       h->masks[i].size); */
  /*     buff += h->masks[i].size; */
  /*   } */
  /* } */
  /* /\* Free the memory *\/ */
  /* buff -= record_size; */
  /* free(buff); */

  /* /\* Get the time of current record. */
  /*    This check is required for the manipulation of the file before */
  /*    the initialization of the time_array. *\/ */
  /* if (times->size != 0) { */
  /*   part->time = time_array_get_time(times, offset); */
  /* } else */
  /*   part->time = -1; */

  /* /\* update the offset. *\/ */
  /* offset = (size_t)((char *)map - (char *)reader->log.log.map); */

  /* /\* Check if an interpolation is required. *\/ */
  /* if (reader_type == logger_reader_const) return offset; */

  /* /\* Start reading next record. *\/ */
  /* struct logger_particle part_next; */

  /* /\* Check that the offset are in the correct direction. *\/ */
  /* if (!header_is_forward(h)) { */
  /*   error("Cannot read a particle with non forward offsets."); */
  /* } */

  /* /\* No next particle. *\/ */
  /* if (h_offset == 0) return (size_t)((char *)map - (char *)reader->log.log.map); */

  /* /\* get absolute offset of next particle. *\/ */
  /* h_offset += offset - header_get_record_size_from_mask(h, mask) - */
  /*             LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE; */

  /* /\* Get time of next record. *\/ */
  /* part_next.time = time_array_get_time(times, h_offset); */

  /* /\* Read next record. *\/ */
  /* h_offset = logger_particle_read(&part_next, reader, h_offset, part_next.time, */
  /*                                 logger_reader_const); */

  /* /\* Interpolate the two particles. *\/ */
  /* *part = logger_particle_interpolate(part, &part_next, time); */

  /* return offset; */
}

/**
 * @brief Read a particle entry in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param id_masks_wanted Indices (in #gravity_logger_fields) of the fields to return.
 * They are assumed to be sorted.
 * @param n_field Number of fields requested.
 * @param output Buffer where the data are written.
 *
 * @return position after the record.
 */
size_t logger_gparticle_read(const struct logger_reader *reader, size_t offset,
                             const int *id_masks_wanted,
                             const int n_mask_wanted,
                             void **output, size_t *mask,
                             size_t *h_offset) {
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

  /* Read the record and copy it to a particle */
  int current_mask = 0;
  for (int i = 0; i < gravity_logger_field_count; i++) {
    const int mask_index = gravity_logger_mask_id[i];

    /* Is the mask present? */
    if (!(*mask & h->masks[mask_index].mask)) {
      /* Requested mask not present, go to the next. */
      if (current_mask < n_mask_wanted && id_masks_wanted[current_mask] == i) {
        current_mask += 1;
      }
      continue;
    }

    /* Read the data if needed and update the buffer position? */
    if (current_mask < n_mask_wanted && id_masks_wanted[current_mask] == i) {
      map = logger_loader_io_read_data(map, h->masks[mask_index].size, output[current_mask]);
      current_mask += 1;
    }
    else {
      map += h->masks[mask_index].size;
    }
  }
  return map - reader->log.log.map;
}

/**
 * @brief Read a particle entry in the log file.
 *
 * @param reader The #logger_reader.
 * @param part The #logger_sparticle to update.
 * @param offset offset of the record to read.
 * @param time time to interpolate (not used if constant interpolation).
 * @param reader_type #logger_reader_type.
 *
 * @return position after the record.
 */
size_t logger_sparticle_read(struct logger_sparticle *part,
                             const struct logger_reader *reader, size_t offset,
                             const double time,
                             const enum logger_reader_type reader_type) {

  /* /\* Save the offset *\/ */
  /* part->offset = offset; */

  /* /\* Get a few pointers. *\/ */
  /* const struct header *h = &reader->log.header; */
  /* void *map = reader->log.log.map; */

  /* const struct time_array *times = &reader->log.times; */

  /* size_t mask = 0; */
  /* size_t h_offset = 0; */

  /* logger_sparticle_init(part); */

  /* /\* Read the record's mask. *\/ */
  /* map = logger_loader_io_read_mask(h, (char *)map + offset, &mask, &h_offset); */

  /* /\* Check that the mask is meaningful *\/ */
  /* if (mask > (unsigned int)(1 << h->masks_count)) { */
  /*   error("Found an unexpected mask %zi", mask); */
  /* } */

  /* /\* Check if it is not a time record. *\/ */
  /* if (mask == h->timestamp_mask) { */
  /*   error("Unexpected timestamp while reading a particle: %lu.", mask); */
  /* } */

  /* /\* Compute the size of the record. *\/ */
  /* size_t record_size = 0; */
  /* for (int i = 0; i < h->masks_count; i++) { */
  /*   if (mask & h->masks[i].mask) { */
  /*     record_size += h->masks[i].size; */
  /*   } */
  /* } */
  /* /\* Read the record and copy it to a particle *\/ */
  /* char *buff = (char *)malloc(record_size); */
  /* map = logger_loader_io_read_data(map, record_size, buff); */
  /* for (int i = 0; i < h->masks_count; i++) { */
  /*   if (mask & h->masks[i].mask) { */
  /*     logger_sparticle_read_single_field(part, buff, h->masks[i].name, */
  /*                                        h->masks[i].size); */
  /*     buff += h->masks[i].size; */
  /*   } */
  /* } */
  /* /\* Free the memory *\/ */
  /* buff -= record_size; */
  /* free(buff); */

  /* /\* Get the time of current record. */
  /*    This check is required for the manipulation of the file before */
  /*    the initialization of the time_array. *\/ */
  /* if (times->size != 0) { */
  /*   part->time = time_array_get_time(times, offset); */
  /* } else */
  /*   part->time = -1; */

  /* /\* update the offset. *\/ */
  /* offset = (size_t)((char *)map - (char *)reader->log.log.map); */

  /* /\* Check if an interpolation is required. *\/ */
  /* if (reader_type == logger_reader_const) return offset; */

  /* /\* Start reading next record. *\/ */
  /* struct logger_sparticle part_next; */

  /* /\* Check that the offset are in the correct direction. *\/ */
  /* if (!header_is_forward(h)) { */
  /*   error("Cannot read a particle with non forward offsets."); */
  /* } */

  /* /\* No next particle. *\/ */
  /* if (h_offset == 0) return (size_t)((char *)map - (char *)reader->log.log.map); */

  /* /\* get absolute offset of next particle. *\/ */
  /* h_offset += offset - header_get_record_size_from_mask(h, mask) - */
  /*             LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE; */

  /* /\* Get time of next record. *\/ */
  /* part_next.time = time_array_get_time(times, h_offset); */

  /* /\* Read next record. *\/ */
  /* h_offset = logger_sparticle_read(&part_next, reader, h_offset, part_next.time, */
  /*                                  logger_reader_const); */

  /* /\* Interpolate the two particles. *\/ */
  /* *part = logger_sparticle_interpolate(part, &part_next, time); */

  /* return offset; */
}
