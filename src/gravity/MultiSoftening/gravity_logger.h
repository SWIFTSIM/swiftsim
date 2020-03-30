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
#ifndef SWIFT_MULTISOFTENING_GRAVITY_LOGGER_H
#define SWIFT_MULTISOFTENING_GRAVITY_LOGGER_H


/**
 * @brief Initialize the logger.
 *
 * WARNING: this should be done in the same order than #gravity_logger_write_particle.
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int gravity_logger_init(struct mask_data *mask_data) {
  mask_data[0] = logger_add_field_to_logger("Coordinates", 3 * sizeof(double));
  mask_data[1] = logger_add_field_to_logger("Velocities", 3 * sizeof(float));
  mask_data[2] = logger_add_field_to_logger("Accelerations", 3 * sizeof(float));
  mask_data[3] = logger_add_field_to_logger("Masses", sizeof(float));
  mask_data[4] = logger_add_field_to_logger("ParticleIDs", sizeof(long long));

  return 5;
}

/**
 * @brief Generates the mask and compute the size of the record.
 *
 * @param masks The list of masks (same order than in #gravity_logger_init).
 * @param gpart The #gpart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 *
 * @return The maximal number of masks that can be written.
 */
INLINE static void gravity_logger_prepare_to_write_particle(
    const struct mask_data *masks, const struct gpart *gpart,
    const int write_all,
    size_t *buffer_size, unsigned int *mask) {

  /* Here you can decide you own writing logic */

  /* Add the coordinates. */
  *mask |= logger_add_field_to_mask(masks[0], "Coordinates", buffer_size);

  /* Add the velocities. */
  *mask |= logger_add_field_to_mask(masks[1], "Velocities", buffer_size);

  /* Add the accelerations. */
  *mask |= logger_add_field_to_mask(masks[2], "Accelerations", buffer_size);

  /* Add the masses. */
  *mask |= logger_add_field_to_mask(masks[3], "Masses", buffer_size);

  /* Add the ID. */
  *mask |= logger_add_field_to_mask(masks[4], "ParticleIDs", buffer_size);
}

/**
 * @brief Write a particle to the logger.
 *
 * @param masks The list of masks (same order than in #gravity_logger_init).
 * @param gp The #gpart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *gravity_logger_write_particle(
    const struct mask_data *mask_data, const struct gpart* gp,
    unsigned int *mask, char *buff) {
#ifdef WITH_LOGGER

  /* Write the coordinate. */
  if (logger_should_write_field(mask_data[0], mask, "Coordinates")) {
    memcpy(buff, gp->x, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Write the velocity. */
  if (logger_should_write_field(mask_data[1], mask, "Velocities")) {
    memcpy(buff, gp->v_full, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the acceleration. */
  if (logger_should_write_field(mask_data[2], mask, "Accelerations")) {
    memcpy(buff, gp->a_grav, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the mass. */
  if (logger_should_write_field(mask_data[3], mask, "Masses")) {
    memcpy(buff, &gp->mass, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the Id. */
  if (logger_should_write_field(mask_data[4], mask, "ParticleIDs")) {
    memcpy(buff, &gp->id_or_neg_offset, sizeof(long long));
    buff += sizeof(long long);
  }

  return buff;

#else
  error("Should not be called without the logger.");
#endif /* WITH_LOGGER */
}

#endif // SWIFT_MULTISOFTENING_GRAVITY_LOGGER_H
