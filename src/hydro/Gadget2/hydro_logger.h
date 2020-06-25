/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GADGET2_HYDRO_LOGGER_H
#define SWIFT_GADGET2_HYDRO_LOGGER_H

#include "logger_io.h"

/**
 * @brief Initialize the logger.
 *
 * WARNING: this should be done in the same order than
 * #hydro_logger_write_particle.
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int hydro_logger_init(struct mask_data *mask_data) {
  mask_data[0] = logger_create_mask_entry("Coordinates", 3 * sizeof(double));
  mask_data[1] = logger_create_mask_entry("Velocities", 3 * sizeof(float));
  mask_data[2] = logger_create_mask_entry("Accelerations", 3 * sizeof(float));
  mask_data[3] = logger_create_mask_entry("Masses", sizeof(float));
  mask_data[4] = logger_create_mask_entry("SmoothingLengths", sizeof(float));
  mask_data[5] = logger_create_mask_entry("Entropies", sizeof(float));
  mask_data[6] = logger_create_mask_entry("ParticleIDs", sizeof(long long));
  mask_data[7] = logger_create_mask_entry("Densities", sizeof(float));

  return 8;
}

/**
 * @brief Generates the mask and compute the size of the record.
 *
 * @param masks The list of masks (same order than in #hydro_logger_init).
 * @param part The #part that will be written.
 * @param xpart The #xpart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void hydro_logger_prepare_to_write_particle(
    const struct mask_data *masks, const struct part *part,
    const struct xpart *xpart, const int write_all, size_t *buffer_size,
    unsigned int *mask) {

  /* Here you can decide your own writing logic */

  /* Add the coordinates. */
  *mask |= logger_add_field_to_mask(masks[0], "Coordinates", buffer_size);

  /* Add the velocities. */
  *mask |= logger_add_field_to_mask(masks[1], "Velocities", buffer_size);

  /* Add the accelerations. */
  *mask |= logger_add_field_to_mask(masks[2], "Accelerations", buffer_size);

  /* Add the masses. */
  *mask |= logger_add_field_to_mask(masks[3], "Masses", buffer_size);

  /* Add the smoothing lengths. */
  *mask |= logger_add_field_to_mask(masks[4], "SmoothingLengths", buffer_size);

  /* Add the entropies. */
  *mask |= logger_add_field_to_mask(masks[5], "Entropies", buffer_size);

  /* Add the ID. */
  *mask |= logger_add_field_to_mask(masks[6], "ParticleIDs", buffer_size);

  /* Add the density. */
  *mask |= logger_add_field_to_mask(masks[7], "Densities", buffer_size);
}

/**
 * @brief Write a particle to the logger.
 *
 * @param masks The list of masks (same order than in #hydro_logger_init).
 * @param p The #part to write.
 * @param xp The #xpart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *hydro_logger_write_particle(
    const struct mask_data *mask_data, const struct part *p,
    const struct xpart *xp, unsigned int *mask, char *buff) {
#ifdef WITH_LOGGER

  /* Write the coordinate. */
  if (logger_should_write_field(mask_data[0], mask, "Coordinates")) {
    memcpy(buff, p->x, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Write the velocity. */
  if (logger_should_write_field(mask_data[1], mask, "Velocities")) {
    memcpy(buff, p->v, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the acceleration. */
  if (logger_should_write_field(mask_data[2], mask, "Accelerations")) {
    float acc[3] = {
        p->a_hydro[0] + xp->a_grav[0],
        p->a_hydro[1] + xp->a_grav[1],
        p->a_hydro[2] + xp->a_grav[2],
    };
    memcpy(buff, acc, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the mass. */
  if (logger_should_write_field(mask_data[3], mask, "Masses")) {
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the smoothing length. */
  if (logger_should_write_field(mask_data[4], mask, "SmoothingLengths")) {
    memcpy(buff, &p->h, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the entropy. */
  if (logger_should_write_field(mask_data[5], mask, "Entropies")) {
    memcpy(buff, &p->entropy, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the Id. */
  if (logger_should_write_field(mask_data[6], mask, "ParticleIDs")) {
    memcpy(buff, &p->id, sizeof(long long));
    buff += sizeof(long long);
  }

  /* Write the density. */
  if (logger_should_write_field(mask_data[7], mask, "Densities")) {
    memcpy(buff, &p->rho, sizeof(float));
    buff += sizeof(float);
  }

  return buff;

#else
  error("Should not be called without the logger.");
#endif /* WITH_LOGGER */
}

#endif  // SWIFT_GADGET2_HYDRO_LOGGER_H
