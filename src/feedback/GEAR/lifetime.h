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
#ifndef SWIFT_LIFETIME_GEAR_H
#define SWIFT_LIFETIME_GEAR_H

#include "hdf5_functions.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Compute the lifetime of a star.
 *
 * @param life The #lifetime model.
 * @param log_mass The star's mass (in log10).
 * @param metallicity The star's metallicity.
 *
 * @return The star's lifetime (in log10).
 */
__attribute__((always_inline)) INLINE static float lifetime_get_log_lifetime_from_mass(
    const struct lifetime *life, float log_mass, float metallicity) {

  /* Compute quadratic term */
  const float quadratic = (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity + life->quadratic[2];
  /* Compute linear term */
  const float linear = (life->linear[0] * metallicity + life->linear[1]) * metallicity + life->linear[2];
  /* Compute constant term */
  const float constant = (life->constant[0] * metallicity + life->constant[1]) * metallicity + life->constant[2];

  /* Apply unit change */
  const float constant_internal_units = constant + (quadratic * life->log_unit_mass - linear) * life->log_unit_mass;
  const float linear_internal_units = linear - 2. * quadratic * life->log_unit_mass;

  /* Compute lifetime */
  return (quadratic * log_mass + linear_internal_units) * log_mass + constant_internal_units;
}

/**
 * @brief Compute the mass of a star with a given lifetime
 * Returns -1 if out of range.
 *
 * @param life The #lifetime model.
 * @param log_time The star's lifetime (in log10).
 * @param metallicity The star's metallicity.
 *
 * @return The star's mass (in log10) or -1.
 */
__attribute__((always_inline)) INLINE static float lifetime_get_log_mass_from_lifetime(
    const struct lifetime *life, float log_time, float metallicity) {

  /* Compute quadratic term */
  const float quadratic = (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity + life->quadratic[2];
  /* Compute linear term */
  const float linear = (life->linear[0] * metallicity + life->linear[1]) * metallicity + life->linear[2];
  /* Compute constant term */
  const float constant = (life->constant[0] * metallicity + life->constant[1]) * metallicity + life->constant[2];

  /* Apply unit change */
  const float constant_internal_units = constant + (quadratic * life->log_unit_mass - linear) * life->log_unit_mass;
  const float linear_internal_units = linear - 2. * quadratic * life->log_unit_mass;

  /* Compute the "c" with the time */
  const float c_t = constant_internal_units - log_time;

  /* Use the quadratic formula to find the mass */
  if (quadratic != 0) {
    const float delta = linear_internal_units * linear_internal_units - 4 * quadratic * c_t;

    /* Avoid complex number should not happen in real simulation */
    if (delta < 0) {
      return - linear_internal_units / (2. * quadratic);
    }
    else {
      return (-linear_internal_units - sqrt(delta)) / (2. * quadratic);
    }
  }
  else {
    return - c_t / linear_internal_units;
  }
}


/**
 * @brief Read lifetime parameters from tables.
 *
 * @param lt The #lifetime.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void lifetime_read_from_tables(
    struct lifetime* lt, struct swift_params* params) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(params, "Data/LiveTimes", &file_id, &group_id);

  /* Allocate the temporary array */
  float *tmp;
  if ((tmp = (float *)malloc(sizeof(float) * 9)) == NULL)
    error("Failed to allocate the temporary array.");

  /* Read the coefficients */
  io_read_array_dataset(group_id, "coeff_z", FLOAT, tmp, 9);
  
  /* Copy the coefficents */
  const int dim = 3;
  for(int i = 0; i < dim; i++) {
    lt->quadratic[i] = tmp[i];
    lt->linear[i] = tmp[i+dim];
    lt->constant[i] = tmp[i+2*dim];
  }

  /* Cleanup everything */
  free(tmp);
  h5_close_group(file_id, group_id);

}

/**
 * @brief Read lifetime parameters from params.
 *
 * @param lt The #lifetime.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void lifetime_read_from_params(
    struct lifetime* lt, struct swift_params* params) {

  /* Read quadratic terms */
  parser_get_opt_param_float_array(params, "GEARLifetime:quadratic", 3, lt->quadratic);

  /* Read linear terms */
  parser_get_opt_param_float_array(params, "GEARLifetime:linear", 3, lt->linear);

  /* Read constant terms */
  parser_get_opt_param_float_array(params, "GEARLifetime:constant", 3, lt->constant);

}

/**
 * @brief Inititialize the Lifetime.
 *
 * @param lt The #lifetime.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void lifetime_init(
    struct lifetime* lt, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

  /* Read params from yields table */
  lifetime_read_from_tables(lt, params);

  /* overwrite the parameters if found in the params */
  lifetime_read_from_params(lt, params);

  /* Change the time unit (mass cannot be done here) */
  lt->constant[2] += log10(phys_const->const_year);

  /* Compute the variable for the change of mass unit */
  lt->log_unit_mass = log10(phys_const->const_solar_mass);
  
}

#endif // SWIFT_LIFETIME_GEAR_H
