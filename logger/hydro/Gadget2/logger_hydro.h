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
#ifndef SWIFT_GADGET2_LOGGER_HYDRO_H
#define SWIFT_GADGET2_LOGGER_HYDRO_H

#include "../config.h"

/* local includes */
#include "logger_loader_io.h"
#include "logger_python_tools.h"

#define logger_part_field_coordinates "Coordinates"
#define logger_part_field_velocities "Velocities"
#define logger_part_field_accelerations "Accelerations"
#define logger_part_field_masses "Masses"
#define logger_part_field_ids "ParticleIDs"
#define logger_part_field_flags "SpecialFlags"
#define logger_part_field_h "SmoothingLengths"
#define logger_part_field_entropies "Entropies"
#define logger_part_field_densities "Densities"

/**
 * @brief Store the data from a record (hydro particle).
 *
 * This structure contains all the required fields
 * present in a file.
 *
 * The particle is initialized with #logger_particle_init
 * and can be updated with a record through #logger_particle_read.
 *
 * In #logger_particle_read, we use #logger_particle_read_single_field on
 * each field and #logger_particle_interpolate if an interpolation is required.
 */
struct logger_particle {

  /* Particle ID. */
  uint64_t id;

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v[3];

  /* Particle acceleration. */
  float a[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle mass. */
  float mass;

  /* Particle density. */
  float rho;

  /* Particle entropy. */
  float entropy;

  /* time of the record. */
  double time;

  /* offset of the particle in the logfile. */
  size_t offset;

  /* Special flag */
  uint32_t flag;
};

/**
 * @brief Print the properties of a #logger_particle.
 *
 * @param p The #logger_particle to print
 */
__attribute__((always_inline)) INLINE static void logger_particle_print(
    const struct logger_particle *p) {
  message("ID:            %lu.", p->id);
  message("Mass:          %g", p->mass);
  message("Time:          %g.", p->time);
  message("Cutoff Radius: %g.", p->h);
  message("Position:     (%g, %g, %g).", p->x[0], p->x[1], p->x[2]);
  message("Velocity:    (%g, %g, %g).", p->v[0], p->v[1], p->v[2]);
  message("Acceleration: (%g, %g, %g).", p->a[0], p->a[1], p->a[2]);
  message("Entropy:       %g.", p->entropy);
  message("Density:       %g.", p->rho);
}

/**
 * @brief Initialize a logger_particle.
 *
 * @param part The #logger_particle to initialize.
 */
__attribute__((always_inline)) INLINE static void logger_particle_init(
    struct logger_particle *part) {
  for (size_t k = 0; k < 3; k++) {
    part->x[k] = 0;
    part->v[k] = 0;
    part->a[k] = 0;
  }

  part->entropy = -1;
  part->rho = -1;
  part->h = -1;
  part->mass = -1;
  part->id = UINT64_MAX;

  part->flag = 0;
}

/**
 * @brief Read a single field for a particle.
 * This function is called multiple times for each record and the buffer is
 * updated each time in order to start at the position of the current field (no
 * need to shift it).
 *
 * @param part The #logger_particle where to write the field.
 * @param buff The buffer to read (starting directly with the current field).
 * @param field Name of the field to read.
 * @param size Number of bits to read.
 */
__attribute__((always_inline)) INLINE static void
logger_particle_read_single_field(struct logger_particle *part, char *buff,
                                  const char *field, const size_t size) {

  /* Copy the buffer to the particle. */
  if (strcmp(logger_part_field_coordinates, field) == 0) {
    memcpy(&part->x, buff, size);
  } else if (strcmp(logger_part_field_velocities, field) == 0) {
    memcpy(&part->v, buff, size);
  } else if (strcmp(logger_part_field_accelerations, field) == 0) {
    memcpy(&part->a, buff, size);
  } else if (strcmp(logger_part_field_entropies, field) == 0) {
    memcpy(&part->entropy, buff, size);
  } else if (strcmp(logger_part_field_h, field) == 0) {
    memcpy(&part->h, buff, size);
  } else if (strcmp(logger_part_field_densities, field) == 0) {
    memcpy(&part->rho, buff, size);
    // TODO link mass and ids together?
  } else if (strcmp(logger_part_field_masses, field) == 0) {
    memcpy(&part->mass, buff, size);
  } else if (strcmp(logger_part_field_ids, field) == 0) {
    memcpy(&part->id, buff, size);
  } else if (strcmp(logger_part_field_flags, field) == 0) {
    memcpy(&part->flag, buff, size);
  } else {
    error("Field %s not defined.", field);
  }
}

/**
 * @brief When starting to read a logfile, check the required fields in the
 * logfile's header.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void logger_particle_check_fields(
    struct header *head) {

  struct logger_particle part;
  for (int i = 0; i < head->masks_count; i++) {
    int size = -1;
    if (strcmp(head->masks[i].name, logger_part_field_coordinates) == 0) {
      size = sizeof(part.x);
    } else if (strcmp(head->masks[i].name, logger_part_field_velocities) == 0) {
      size = sizeof(part.v);
    } else if (strcmp(head->masks[i].name, logger_part_field_accelerations) ==
               0) {
      size = sizeof(part.a);
    } else if (strcmp(head->masks[i].name, logger_part_field_masses) == 0) {
      size = sizeof(part.mass);
    } else if (strcmp(head->masks[i].name, logger_part_field_h) == 0) {
      size = sizeof(part.h);
    } else if (strcmp(head->masks[i].name, logger_part_field_entropies) == 0) {
      size = sizeof(part.entropy);
    } else if (strcmp(head->masks[i].name, logger_part_field_ids) == 0) {
      size = sizeof(part.id);
    } else if (strcmp(head->masks[i].name, logger_part_field_densities) == 0) {
      size = sizeof(part.rho);
    }

    if (size != -1 && size != head->masks[i].size) {
      error("The field size is not compatible for %s (%i %i)",
            head->masks[i].name, size, head->masks[i].size);
    }
  }
}

/**
 * @brief Interpolate the location of the particle at the given time.
 * Here we use a linear interpolation for most of the fields.
 * For the position (velocity), we use a quintic (cubic) hermite interpolation
 * based on the positions, velocities and accelerations at the time of the two
 * particles
 *
 * @param part_bef #logger_particle current particle (before time)
 * @param part_next #logger_particle next particle (after time)
 * @param time interpolation time
 *
 * @return The interpolated particle.
 *
 */
__attribute__((always_inline)) INLINE static struct logger_particle
logger_particle_interpolate(const struct logger_particle *part_bef,
                            const struct logger_particle *part_next,
                            const double time) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the particle order. */
  if (part_next->time < part_bef->time)
    error("Wrong particle order (next before current): %g, %g", part_next->time,
          part_bef->time);
  if ((time < part_bef->time) || (part_next->time < time))
    error(
        "Cannot extrapolate (particle time: %f, "
        "interpolating time: %f, next particle time: %f).",
        part_bef->time, time, part_next->time);
#endif

  /* Copy everything into the return particle */
  struct logger_particle ret = *part_bef;

  /* Compute the interpolation scaling. */
  const double wa =
      (time - part_bef->time) / (part_next->time - part_bef->time);
  const double wb =
      (part_next->time - time) / (part_next->time - part_bef->time);

  /* interpolate vectors. */
  for (size_t i = 0; i < 3; i++) {
    /* position */
    ret.x[i] = logger_tools_quintic_hermite_spline(
        part_bef->time, part_bef->x[i], part_bef->v[i], part_bef->a[i],
        part_next->time, part_next->x[i], part_next->v[i], part_next->a[i],
        time);

    /* velocity */
    ret.v[i] = logger_tools_cubic_hermite_spline(
        part_bef->time, part_bef->v[i], part_bef->a[i], part_next->time,
        part_next->v[i], part_next->a[i], time);

    /* acceleration */
    ret.a[i] = wa * part_bef->a[i] + wb * part_next->a[i];
  }

  /* mass */
  ret.mass = wa * part_bef->mass + wb * part_next->mass;

  /* density */
  ret.rho = wa * part_bef->rho + wb * part_next->rho;

  /* entropy */
  ret.entropy = wa * part_bef->entropy + wb * part_next->entropy;

  /* smoothing length */
  ret.h = wa * part_bef->h + wb * part_next->h;

  /* set time. */
  ret.time = time;

  return ret;
}

#ifdef HAVE_PYTHON
/**
 * @brief Specifies which particle fields to give access to in python
 *
 * @param list The list of structure used to store the information.
 *
 * @return The number of i/o fields to give access to.
 */
INLINE static int logger_particles_generate_python(
    struct logger_python_field *list) {

  struct logger_particle *part;

  /* List what we want to use in python */
  list[0] =
      logger_loader_python_field(logger_part_field_coordinates, part, x, "3f8");

  list[1] =
      logger_loader_python_field(logger_part_field_velocities, part, v, "3f4");

  // TODO sum the grav + hydro accelerations
  list[2] = logger_loader_python_field(logger_part_field_accelerations, part, a,
                                       "3f4");

  list[3] =
      logger_loader_python_field(logger_part_field_masses, part, mass, "f4");

  list[4] = logger_loader_python_field(logger_part_field_h, part, h, "f4");

  list[5] = logger_loader_python_field(logger_part_field_entropies, part,
                                       entropy, "f4");

  list[6] = logger_loader_python_field(logger_part_field_ids, part, id, "i8");

  list[7] =
      logger_loader_python_field(logger_part_field_densities, part, rho, "f4");

  list[8] = logger_loader_python_field("Times", part, time, "f8");

  return 9;
}
#endif /* HAVE_PYTHON */

#endif /* SWIFT_GADGET2_LOGGER_HYDRO_H */
