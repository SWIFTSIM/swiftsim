/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_LOGGER_STARS_H
#define SWIFT_DEFAULT_LOGGER_STARS_H

#include "../config.h"

/* local includes */
#include "logger_loader_io.h"
#include "logger_python_tools.h"

/**
 * @brief Store the data for a record (stars particle).
 *
 * This structure contains all the required fields
 * present in a file.
 *
 * The particle is initialized with #logger_sparticle_init
 * and can be updated with a record through #logger_sparticle_read.
 *
 * In #logger_sparticle_read, we use #logger_sparticle_read_field on
 * each field and #logger_sparticle_interpolate if an interpolation is required.
 */
struct logger_sparticle {

  /* Particle ID. */
  long long id;

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v[3];

  /* Particle acceleration. */
  float a[3];

  /* Particle mass. */
  float mass;

  /* Particle cutoff radius. */
  float h;

  /* time of the record. */
  double time;

  /* offset of the particle */
  size_t offset;

  /* Special flag */
  uint32_t flag;
};

/**
 * @brief Print the properties of a #logger_sparticle.
 *
 * @param p The #logger_sparticle to print
 */
__attribute__((always_inline)) INLINE static void logger_sparticle_print(
    const struct logger_sparticle *p) {
  message("ID:              %lli.", p->id);
  message("Mass:            %g", p->mass);
  message("SmoothingLength: %g", p->h);
  message("Time:            %g.", p->time);
  message("Positions:       (%g, %g, %g).", p->x[0], p->x[1], p->x[2]);
  message("Velocities:      (%g, %g, %g).", p->v[0], p->v[1], p->v[2]);
  message("Accelerations:   (%g, %g, %g).", p->a[0], p->a[1], p->a[2]);
}

/**
 * @brief Initialize a #logger_sparticle.
 *
 * @param part The #logger_sparticle to initialize.
 */
__attribute__((always_inline)) INLINE static void logger_sparticle_init(
    struct logger_sparticle *part) {
  for (size_t k = 0; k < 3; k++) {
    part->x[k] = 0;
    part->v[k] = 0;
    part->a[k] = 0;
  }

  part->mass = -1;
  part->h = -1;
  part->id = SIZE_MAX;

  part->flag = 0;
}

/**
 * @brief Read a single named entry for a particle.
 *
 * @param part The #logger_sparticle to update.
 * @param buff The buffer containing the particle.
 * @param field field to read.
 * @param size number of bits to read.
 *
 * @return mapped data after the block read.
 */
__attribute__((always_inline)) INLINE static void logger_sparticle_read_field(
    struct logger_sparticle *part, char *buff, const char *field,
    const size_t size) {

  /* Copy the buffer to the particle. */
  if (strcmp("Coordinates", field) == 0) {
    memcpy(&part->x, buff, size);
  } else if (strcmp("Velocities", field) == 0) {
    memcpy(&part->v, buff, size);
  } else if (strcmp("Accelerations", field) == 0) {
    memcpy(&part->a, buff, size);
    // TODO link mass and ids together?
  } else if (strcmp("Masses", field) == 0) {
    memcpy(&part->mass, buff, size);
  } else if (strcmp("SmoothingLengths", field) == 0) {
    memcpy(&part->h, buff, size);
  } else if (strcmp("ParticleIDs", field) == 0) {
    memcpy(&part->id, buff, size);
  } else if (strcmp("SpecialFlags", field) == 0) {
    memcpy(&part->flag, buff, size);
  } else {
    error("Type %s not defined.", field);
  }
}

/**
 * @brief Check if the required fields are correct.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void logger_sparticle_check_fields(
    struct header *head) {

  struct logger_sparticle part;
  for(int i = 0; i < head->masks_count; i++) {
    int size = -1;
    if (strcmp(head->masks[i].name, "Positions") == 0) {
      size = sizeof(part.x);
    }
    else if (strcmp(head->masks[i].name, "Velocities") == 0) {
      size = sizeof(part.v);
    }
    else if (strcmp(head->masks[i].name, "Accelerations") == 0) {
      size = sizeof(part.a);
    }
    else if (strcmp(head->masks[i].name, "Masses") == 0) {
      size = sizeof(part.mass);
    }
    else if (strcmp(head->masks[i].name, "SmoothingLengths") == 0) {
      size = sizeof(part.h);
    }
    else if (strcmp(head->masks[i].name, "ParticleIDs") == 0) {
      size = sizeof(part.id);
    }

    if (size != -1 && size != head->masks[i].size) {
      error("The field size is not compatible for %s (%i %i)",
            head->masks[i].name, size, head->masks[i].size);
    }

  }
}

/**
 * @brief interpolate two particles at a given time
 *
 * @param part_curr #logger_sparticle In: current particle (before time), Out:
 * interpolated particle
 * @param part_next #logger_sparticle next particle (after time)
 * @param time interpolation time
 *
 */
__attribute__((always_inline)) INLINE static void logger_sparticle_interpolate(
    struct logger_sparticle *part_curr,
    const struct logger_sparticle *part_next, const double time) {

  /* Check that a particle is provided. */
  if (!part_curr) error("part_curr is NULL.");
  if (!part_next) error("part_next is NULL.");

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the particle order. */
  if (part_next->time < part_curr->time)
    error("Wrong particle order (next before current): %g, %g", part_next->time,
          part_curr->time);
  if ((time < part_curr->time) || (part_next->time < time))
    error(
        "Cannot extrapolate (particle time: %f, "
        "interpolating time: %f, next particle time: %f).",
        part_curr->time, time, part_next->time);
#endif

  /* Compute the interpolation scaling. */
  double scaling = part_next->time - part_curr->time;

  scaling = (time - part_curr->time) / scaling;

  /* interpolate vectors. */
  for (size_t i = 0; i < 3; i++) {
    part_curr->x[i] = logger_tools_quintic_hermite_spline(
      part_curr->time, part_curr->x[i], part_curr->v[i], part_curr->a[i],
      part_next->time, part_next->x[i], part_next->v[i], part_next->a[i], time);

    part_curr->v[i] = logger_tools_cubic_hermite_spline(
      part_curr->time, part_curr->v[i], part_curr->a[i],
      part_next->time, part_next->v[i], part_next->a[i], time);

    const float tmp = (part_next->a[i] - part_curr->a[i]);
    part_curr->a[i] += tmp * scaling;
  }

  float tmp = (part_next->mass - part_curr->mass);
  part_curr->mass += tmp * scaling;

  tmp = (part_next->h - part_curr->h);
  part_curr->h += tmp * scaling;

  /* set time. */
  part_curr->time = time;
}

#ifdef HAVE_PYTHON
/**
 * @brief Specifies which particle fields to give access to in python
 *
 * @param list The list of structure used to store the information.
 *
 * @param num_fields The number of i/o fields to give access to.
 */
INLINE static int logger_sparticles_generate_python(
    struct logger_python_field *list) {
  struct logger_sparticle *part;

  /* List what we want to use in python */
  list[0] = logger_loader_python_field("Coordinates", part, x, 3, NPY_DOUBLE);

  list[1] = logger_loader_python_field("Velocities", part, v, 3, NPY_FLOAT32);

  list[2] = logger_loader_python_field("Accelerations", part, a, 3, NPY_FLOAT32);

  list[3] = logger_loader_python_field("Masses", part, mass, 1, NPY_FLOAT32);

  list[4] =
      logger_loader_python_field("SmoothingLengths", part, h, 1, NPY_FLOAT32);

  list[5] =
      logger_loader_python_field("ParticleIDs", part, id, 1, NPY_LONGLONG);

  list[6] = logger_loader_python_field("Times", part, time, 1, NPY_DOUBLE);

  return 7;
}
#endif /* HAVE_PYTHON */

#endif  // SWIFT_DEFAULT_LOGGER_GRAVITY_H
