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
#include "logger_io.h"
#include "logger_time.h"
#include "logger_tools.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Print the properties of a logger_particle.
 *
 * @param p The #logger_particle to print
 */
void logger_particle_print(const struct logger_particle *p) {
  message("ID:            %lu\n", p->id);
  message("Mass:          %g\n", p->mass);
  message("Time:          %g\n", p->time);
  message("Cutoff Radius: %g\n", p->h);
  message("Positions:     (%g, %g, %g)\n", p->pos[0], p->pos[1], p->pos[2]);
  message("Velocities:    (%g, %g, %g)\n", p->vel[0], p->vel[1], p->vel[2]);
  message("Accelerations: (%g, %g, %g)\n", p->acc[0], p->acc[1], p->acc[2]);
  message("Entropy:       %g\n", p->entropy);
  message("Density:       %g\n", p->density);
}

/**
 * @brief Initialize a logger_particle.
 *
 * @param part The #logger_particle to initialize.
 */
void logger_particle_init(struct logger_particle *part) {
  for (size_t k = 0; k < DIM; k++) {
    part->pos[k] = 0;
    part->vel[k] = 0;
    part->acc[k] = 0;
  }

  part->entropy = -1;
  part->density = -1;
  part->h = -1;
  part->mass = -1;
  part->id = SIZE_MAX;
}

/**
 * @brief Read a single field for a particle
 *
 * @param part The #logger_particle to update
 * @param data Pointer to the data to read.
 * @param offset In: read position, Out: input shifted by the required amount of
 * data
 * @param field field to read
 * @param size number of bits to read
 *
 */
void logger_particle_read_field(struct logger_particle *part, void *map,
                                size_t *offset, const char *field,
                                const size_t size) {
  void *p = NULL;

  if (strcmp("positions", field) == 0) {
    p = &part->pos;
  } else if (strcmp("velocities", field) == 0) {
    p = &part->vel;
  } else if (strcmp("accelerations", field) == 0) {
    p = &part->acc;
  } else if (strcmp("entropy", field) == 0) {
    p = &part->entropy;
  } else if (strcmp("smoothing length", field) == 0) {
    p = &part->h;
  } else if (strcmp("density", field) == 0) {
    p = &part->density;
  } else if (strcmp("consts", field) == 0) {
    p = malloc(size);
  } else {
    error("Type %s not defined", field);
  }

  io_read_data(map, size, p, offset);

  if (strcmp("consts", field) == 0) {
    part->mass = 0;
    part->id = 0;
    memcpy(&part->mass, p, sizeof(float));
    p += sizeof(float);
    memcpy(&part->id, p, sizeof(size_t));
    p -= sizeof(float);
    free(p);
  }
}

/**
 * @brief Read a particle in the dump file.
 *
 * @param part The #logger_particle to update.
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset offset of the chunk to read (update it to the end of the chunk)
 * @param time time to interpolate (if #logger_reader_type is an interpolating
 * one)
 * @param reader #logger_reader_type
 * @param times #time_array times in the dump
 *
 */
void logger_particle_read(struct logger_particle *part, const struct header *h,
                          void *map, size_t *offset, const double time,
                          const int reader, struct time_array *times) {
  size_t mask = 0;
  size_t h_offset = 0;

  logger_particle_init(part);

  io_read_mask(h, map, offset, &mask, &h_offset);

  if (mask != 127) error("Unexpected mask: %lu", mask);

  for (size_t i = 0; i < h->number_mask; i++) {
    if (mask & h->masks[i].mask) {
      logger_particle_read_field(part, map, offset, h->masks[i].name,
                                 h->masks[i].size);
    }
  }

  if (times) /* move offset by 1 in order to be in the required chunk */
    part->time = time_array_get_time(times, *offset - 1);
  else
    part->time = -1;

  /* end of const case */
  if (reader == logger_reader_const) return;

  /* read next particle */
  struct logger_particle part_next;

  if (!h->forward_offset) error("TODO");

  if (h_offset == 0) return;
  /* get absolute offset of next particle */
  h_offset += *offset - header_get_mask_size(h, mask) - LOGGER_MASK_SIZE -
              LOGGER_OFFSET_SIZE;

  part_next.time = time_array_get_time(times, h_offset);

  /* previous part exists */
  logger_particle_read(&part_next, h, map, &h_offset, part_next.time,
                       logger_reader_const, times);

  logger_particle_interpolate(part, &part_next, time);
}

/**
 * @brief interpolate two particles at a given time
 *
 * @param part_curr #logger_particle In: current particle (before time), Out:
 * interpolated particle
 * @param part_next #logger_particle next particle (after time)
 * @param time interpolation time
 *
 */
void logger_particle_interpolate(struct logger_particle *part_curr,
                                 const struct logger_particle *part_next,
                                 const double time) {

  if (!part_curr) error("part_curr is NULL");
  if (!part_next) error("part_next is NULL");

#ifdef SWIFT_DEBUG_CHECKS
  if (part_next->time <= part_curr->time)
    error("Wrong particle order (next before current)");
  if ((time < part_curr->time) || (part_next->time < time))
    error(
        "Interpolating, not extrapolating (particle time: %f, "
        "interpolating time: %f, next particle time: %f)",
        part_curr->time, time, part_next->time);
#endif

  double scaling = part_next->time - part_curr->time;

  scaling = (time - part_curr->time) / scaling;

  double tmp;
  float ftmp;

  /* interpolate vectors */
  for (size_t i = 0; i < DIM; i++) {
    tmp = (part_next->pos[i] - part_curr->pos[i]);
    part_curr->pos[i] += tmp * scaling;

    ftmp = (part_next->vel[i] - part_curr->vel[i]);
    part_curr->vel[i] += ftmp * scaling;

    ftmp = (part_next->acc[i] - part_curr->acc[i]);
    part_curr->acc[i] += ftmp * scaling;
  }

  /* interpolate scalars */
  ftmp = (part_next->entropy - part_curr->entropy);
  part_curr->entropy += ftmp * scaling;

  /* set time */
  part_curr->time = time;
}
