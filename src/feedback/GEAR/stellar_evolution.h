/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STELLAR_EVOLUTION_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_GEAR_H

#include "stellar_evolution_struct.h"

#include <math.h>
#include <stddef.h>


/**
 * @brief Compute the lifetime of a star.
 *
 * @param life The #lifetime model.
 * @param log_mass The star's mass (in log10).
 * @param metallicity The star's metallicity.
 *
 * @return The star's lifetime (in log10).
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_log_lifetime(
    const struct lifetime *life, float log_mass, float metallicity) {
  // TODO units

  /* Compute quadratic term */
  const float a_m2 = (life->a_m2[0] * metallicity + life->a_m2[1]) * metallicity + life->a_m2[2];
  /* Compute linear term */
  const float b_m = (life->b_m[0] * metallicity + life->b_m[1]) * metallicity + life->b_m[2];
  /* Compute constant term */
  const float c = (life->c[0] * metallicity + life->c[1]) * metallicity + life->c[2];

  /* Compute lifetime */
  return (a_m2 * log_mass + b_m) * log_mass + c;
}

/**
 * @brief Compute the mass of a star with a given lifetime
 *
 * @param life The #lifetime model.
 * @param log_time The star's lifetime (in log10).
 * @param metallicity The star's metallicity.
 *
 * @return The star's mass (in log10).
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_log_mass_from_lifetime(
    const struct lifetime *life, float log_time, float metallicity) {

  // TODO units

  /* Compute quadratic term */
  const float a_m2 = (life->a_m2[0] * metallicity + life->a_m2[1]) * metallicity + life->a_m2[2];
  /* Compute linear term */
  const float b_m = (life->b_m[0] * metallicity + life->b_m[1]) * metallicity + life->b_m[2];
  /* Compute constant term */
  const float c = (life->c[0] * metallicity + life->c[1]) * metallicity + life->c[2];

  /* Compute the "c" with the time */
  const float c_t = c - log_time;

  /* Use the quadratic formula to find the mass */
  if (a_m2 != 0) {
    return (-b_m - sqrt(b_m * b_m - 4 * a_m2 * c_t)) / (2. * a_m2);
  }
  else {
    return - c_t / b_m;
  }
}

__attribute__((always_inline)) INLINE static float stellar_evolution_get_imf(void) {
  return 0.;
};

__attribute__((always_inline)) INLINE static float stellar_evolution_get_number_integrated_imf(void) {
  return 0.;
};

__attribute__((always_inline)) INLINE static float stellar_evolution_get_mass_integrated_imf(void) {
  return 0.;
};

__attribute__((always_inline)) INLINE static float stellar_evolution_get_supernovae_ia_rate(void) {
  return 0.;
};

__attribute__((always_inline)) INLINE static float stellar_evolution_get_supernovae_ii_rate(void) {
  return 0.;
};

__attribute__((always_inline)) INLINE static float *stellar_evolution_get_supernovae_ia_yields(void) {
  return (float*)NULL;
};

__attribute__((always_inline)) INLINE static float *stellar_evolution_get_supernovae_ii_yields(void) {
  return (float*)NULL;
};

#endif // SWIFT_STELLAR_EVOLUTION_GEAR_H
