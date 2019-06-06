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

/**
 * @brief Compute the lifetime of a star.
 *
 * @param life The #lifetime model.
 * @param log_mass The star's mass (in log10).
 * @param metallicity The star's metallicity.
 *
 * @return The star's lifetime (in log10).
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_log_lifetime_from_mass(
    const struct lifetime *life, float log_mass, float metallicity) {
  // TODO units

  /* Compute quadratic term */
  const float quadratic = (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity + life->quadratic[2];
  /* Compute linear term */
  const float linear = (life->linear[0] * metallicity + life->linear[1]) * metallicity + life->linear[2];
  /* Compute constant term */
  const float constant = (life->constant[0] * metallicity + life->constant[1]) * metallicity + life->constant[2];

  /* Compute lifetime */
  return (quadratic * log_mass + linear) * log_mass + constant;
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
__attribute__((always_inline)) INLINE static float stellar_evolution_get_log_mass_from_lifetime(
    const struct lifetime *life, float log_time, float metallicity) {

  // TODO units

  /* Compute quadratic term */
  const float quadratic = (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity + life->quadratic[2];
  /* Compute linear term */
  const float linear = (life->linear[0] * metallicity + life->linear[1]) * metallicity + life->linear[2];
  /* Compute constant term */
  const float constant = (life->constant[0] * metallicity + life->constant[1]) * metallicity + life->constant[2];

  /* Compute the "c" with the time */
  const float c_t = constant - log_time;

  /* Use the quadratic formula to find the mass */
  if (quadratic != 0) {
    const float delta = linear * linear - 4 * quadratic * c_t;

    /* Avoid complex number should not happen in real simulation */
    if (delta < 0) {
      return - linear / (2. * quadratic);
    }
    else {
      return (-linear - sqrt(delta)) / (2. * quadratic);
    }
  }
  else {
    return - c_t / linear;
  }
}

/**
 * @brief Inititialize the Lifetime.
 *
 * @param lt The #lifetime.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_init_lifetime(
    struct lifetime* lt, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

  /* Read quadratic terms */
  parser_get_param_float_array(params, "GEARLifetime:quadratic", 3, lt->quadratic);

  /* Read linear terms */
  parser_get_param_float_array(params, "GEARLifetime:linear", 3, lt->linear);
  // b -> b - 2 a log(Msun)

  /* Read constant terms */
  parser_get_param_float_array(params, "GEARLifetime:constant", 3, lt->constant);
  // c -> c + a log(Msun)^2 - b log(Msun) + log(Tyr)
  
}
#endif // SWIFT_LIFETIME_GEAR_H
