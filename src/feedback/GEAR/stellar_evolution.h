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
  const float quadratic = (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity + life->quadratic[2];
  /* Compute linear term */
  const float linear = (life->linear[0] * metallicity + life->linear[1]) * metallicity + life->linear[2];
  /* Compute constant term */
  const float constant = (life->constant[0] * metallicity + life->constant[1]) * metallicity + life->constant[2];

  /* Compute the "c" with the time */
  const float c_t = constant - log_time;

  /* Use the quadratic formula to find the mass */
  if (quadratic != 0) {
    return (-linear - sqrt(linear * linear - 4 * quadratic * c_t)) / (2. * quadratic);
  }
  else {
    return - c_t / linear;
  }
}

/**
 * @brief Compute the mass fraction of the initial mass function.
 *
 * @param imf The #initial_mass_function.
 * @param m The mass to evaluate.
 *
 * @return The mass fraction.
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_imf(
    const struct initial_mass_function *imf, float m) {

#ifdef SWIFT_DEBUG_CHECKS
  if (m > imf->mass_max || m < imf->mass_min)
    error("Mass below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m, imf->mass_max);
#endif

  for(int i = 0; i < imf->n_parts; i++) {
    if (m <= imf->mass_limits[i+1]) {
      return imf->coef[i] * pow(m, imf->exp[i]);
    }
  }

  error("Failed to find correct function part: %g larger than mass max %g.",
	m, imf->mass_max);
};

/**
 * @brief Compute the fraction of stars (in number) in a given range of mass.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lowest mass.
 * @param m2 The largest mass.
 *
 * @return The fraction stars in the interval (in number).
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_imf_number(
    const struct initial_mass_function *imf, float m1, float m2) {
  error("This has not been tested");
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > imf->mass_max || m1 < imf->mass_min)
    error("Mass 1 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m1, imf->mass_max);

  if (m2 > imf->mass_max || m2 < imf->mass_min)
    error("Mass 2 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m2, imf->mass_max);
  
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  float n = 0.;

  for(int i = 0; i < imf->n_parts; i++) {
    /* Are we above the lowest mass? */
    if (m1 > imf->mass_limits[i+1])
      continue;

    /* Are we after the above the largest mass? */
    if (m2 < imf->mass_limits[i+1])
      break;

    /* Get integral limits. */
    const float mass_min = (m1 > imf->mass_limits[i]) ? m1: imf->mass_limits[i];
    const float mass_max = (m2 < imf->mass_limits[i+1]) ? m2: imf->mass_limits[i+1];
    const float exp = imf->exp[i];

    /* Compute the contribution of the current part. */
    n += imf->coef[i] * (pow(mass_max, exp) - pow(mass_min, exp)) / exp;
  }

  return n;
};

/**
 * @brief Compute the fraction of stars (in mass) in a given range of mass.
 *
 * @param imf The #initial_mass_function.
 * @param m1 The lowest mass.
 * @param m2 The largest mass.
 *
 * @return The fraction stars in the interval (in mass).
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_imf_mass(
    const struct initial_mass_function *imf, float m1, float m2) {
  error("This has not been tested");
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > imf->mass_max || m1 < imf->mass_min)
    error("Mass 1 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m1, imf->mass_max);

  if (m2 > imf->mass_max || m2 < imf->mass_min)
    error("Mass 2 below or above limits expecting %g < %g < %g.",
	  imf->mass_min, m2, imf->mass_max);
  
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  float mass = 0.;

  for(int i = 0; i < imf->n_parts; i++) {
    /* Are we above the lowest mass? */
    if (m1 > imf->mass_limits[i+1])
      continue;

    /* Are we after the above the largest mass? */
    if (m2 < imf->mass_limits[i+1])
      break;

    /* Get integral limits. */
    const float mass_min = (m1 > imf->mass_limits[i]) ? m1: imf->mass_limits[i];
    const float mass_max = (m2 < imf->mass_limits[i+1]) ? m2: imf->mass_limits[i+1];
    const float exp = imf->exp[i] + 1.;

    /* Compute the contribution of the current part. */
    mass += imf->coef[i] * (pow(mass_max, exp) - pow(mass_min, exp)) / exp;
  }

  return mass;
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
