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
#ifndef SWIFT_MULTISOFTENING_LOGGER_GRAVITY_H
#define SWIFT_MULTISOFTENING_LOGGER_GRAVITY_H

#include "../config.h"
/**
   TODO

   Merge everything into writer
   Use a single enum for both reader / writer
   Use particle from writer
   Check_fields -> populate
   read_gparticle on a single field at a time
   Use NaN to flag available fields
   Remove available_fields
 */

/* local includes */
#include "logger_loader_io.h"
#include "logger_python_tools.h"

#define logger_gpart_field_coordinates "Coordinates"
#define logger_gpart_field_velocities "Velocities"
#define logger_gpart_field_accelerations "Accelerations"
#define logger_gpart_field_masses "Masses"
#define logger_gpart_field_ids "ParticleIDs"
#define logger_gpart_field_flags "SpecialFlags"

#define gravity_mask_count 5

/* Index of the mask in the header mask array */
static int gravity_mask_id[gravity_mask_count];
/* Size of each field */
static int gravity_mask_size[gravity_mask_count];

enum gravity_mask_index {
  gravity_coordinate = 0,
  gravity_velocity = 1,
  gravity_acceleration = 2,
  gravity_mass = 3,
  gravity_id = 4,
};

struct logger_gravity_particle {
  /* Particle ID. */
  uint64_t id;

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v[3];

  /* Particle acceleration. */
  float a[3];

  /* Particle mass. */
  float mass;
};

/**
 * @brief When starting to read a logfile, check the required fields in the
 * logfile's header.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void logger_gparticle_check_fields(
    struct header *head) {

  /* This loop will be moved outside of the module */
  for(int i = 0; i < gravity_mask_count; i++) {
    gravity_mask_id[i] = -1;
  }

  for (int i = 0; i < head->masks_count; i++) {
    if (strcmp(head->masks[i].name, logger_gpart_field_coordinates) == 0) {
      gravity_mask_size[gravity_coordinate] = 3 * sizeof(double);
      gravity_mask_id[gravity_coordinate] = i;

    } else if (strcmp(head->masks[i].name, logger_gpart_field_velocities) ==
               0) {
      gravity_mask_size[gravity_velocity] = 3 * sizeof(float);
      gravity_mask_id[gravity_velocity] = i;

    } else if (strcmp(head->masks[i].name, logger_gpart_field_accelerations) ==
               0) {
      gravity_mask_size[gravity_acceleration] = 3 * sizeof(float);
      gravity_mask_id[gravity_acceleration] = i;

    } else if (strcmp(head->masks[i].name, logger_gpart_field_masses) == 0) {
      gravity_mask_size[gravity_mass] = sizeof(float);
      gravity_mask_id[gravity_mass] = i;

    } else if (strcmp(head->masks[i].name, logger_gpart_field_ids) == 0) {
      gravity_mask_size[gravity_id] = sizeof(uint64_t);
      gravity_mask_id[gravity_id] = i;
    }

  }

  /* This loop will be moved outside of the module */
  for(int i = 0; i < gravity_mask_count; i++) {
    int id = gravity_mask_id[i];
    if (head->masks[id].size != gravity_mask_size[i]) {
      error("Size are not compatible for the field %s.",
            head->masks[i].name);
    }
  }
}


/**
 * @param mask_id Mask to read (index in header)
 * @param n_mask number of mask to read
 * @param buffer The file
 * @param fields Output: will see format later
 * @param mask Mask of the current record
 */
void read_gparticle(const int *mask_id, const int n_mask, char *buffer, void *fields,
                    unsigned int mask_record) {
  for(int i = 0; i < head->masks_count; i++) {
    unsigned int mask = 1 << i;

    /* Can we skip this mask? */
    if (!(mask_record & mask)) {
      continue;
    }

    /* Check if we want to read it */
    int j;
    for(j = 0; j < n_mask; j++) {
      if (mask_id[j] == i) {
        break;
      }
    }

    /* Copy the field if we need it. */
    if (j != n_mask) {
      memcpy(fields, buffer, head->masks[i].size);
    }

    /* Move the buffer */
    buffer += head->masks.size[i]
  }
}

/**
 * @brief Interpolate the location of the particle at the given time.
 * Here we use a linear interpolation for most of the fields.
 * For the position (velocity), we use a quintic (cubic) hermite interpolation
 * based on the positions, velocities and accelerations at the time of the two
 * particles
 *
 * @param part_bef #logger_gparticle current particle (before time)
 * @param part_next #logger_gparticle next particle (after time)
 * @param time interpolation time
 *
 * @return The interpolated particle.
 *
 */
__attribute__((always_inline)) INLINE static void
logger_gparticle_interpolate_field(
    const logger_gravity_particle *before, const logger_gravity_particle *after,
    void *output, const double t_before, const double t_after,
    const double t, const int field, const int *available_fields) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the times */
#endif

  /* Compute the interpolation scaling. */
  const double wa =
      (t - t_before) / (t_after - t_before);
  const double wb =
      (t_after - t) / (t_after - t_before);

  switch(field) {
    case gravity_mask_id[gravity_coordinate]:
      if (available_fields[gravity_velocity] &&
          available_fields[gravity_acceleration]) {
        /* interpolate vectors. */
        for (int i = 0; i < 3; i++) {
          /* position */
          ((double *)output)[i] = logger_tools_quintic_hermite_spline(
            t_before, before->x[i], before->v[i], before->a[i],
            t_after->time, after->x[i], after->v[i], after->a[i],
            t);
        }
      }
      else if (available_fields[gravity_velocity]) {
        // TODO
      };
      break;
  }
}

#endif  // SWIFT_MULTISOFTENING_LOGGER_GRAVITY_H
