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
 * @brief Get the IMF exponent in between mass_min and mass_max.
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_imf_exponent(
    const struct initial_mass_function* imf, float mass_min, float mass_max) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mass_max > imf->mass_max)
    error("Cannot have mass larger than the largest one in the IMF");
  if (mass_min < imf->mass_min)
    error("Cannot have mass smaller than the smallest one in the IMF");
  if (mass_max < mass_min)
    error("Cannot have mass_min larger than mass_max");
#endif
  
  for(int i = 0; i < imf->n_parts; i++) {

    /* Check if in the correct part of the IMF */
    if (mass_min < imf->mass_limits[i+1]) {

      /* Check if in only one segment */
      if (mass_max > imf->mass_limits[i+1]) {
	error("The code is not implemented to deal with two different IMF part in the supernovae IMF");
      }

      return imf->exp[i];
    }
  }

  error("Masses outside IMF ranges");

  return -1;
}

/**
 * @brief Get the IMF coefficient in between mass_min and mass_max.
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_imf_coefficient(
    const struct initial_mass_function* imf, float mass_min, float mass_max) {

  for(int i = 0; i < imf->n_parts; i++) {

    /* Check if in the correct part of the IMF */
    if (mass_min < imf->mass_limits[i+1]) {

      /* Check if in only one segment */
      if (mass_max > imf->mass_limits[i+1]) {
	error("The code is not implemented to deal with two different exponent for the IMF in supernovae");
      }

      return imf->coef[i];
    }
  }
  
  error("Masses outside IMF ranges");

  return -1;
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

/**
 * @brief Compute the companion integral (second integral in equation 3.46 in Poirier 2004)
 *
 * @param snia The #supernovae_ia model.
 * @param m1 The lower mass limit.
 * @param m2 The upper mass limit.
 * @param companion_type The type of companion (e.g. index of snia->companion).
 *
 * @return The fraction of companion.
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_companion_fraction(
    struct supernovae_ia *snia, float m1, float m2, int companion_type) {
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  const float tmp = pow(m2, snia->companion_exponent) - pow(m1, snia->companion_exponent);
  return snia->companion[companion_type].coef * tmp / snia->companion_exponent;
  
}

/**
 * @brief Compute the number of supernovae Ia per unit of mass (equation 3.46 in Poirier 2004).
 *
 * @param snia The #supernovae_ia model.
 * @param m1 The lower mass limit.
 * @param m2 The upper mass limit.
 *
 * @return The number of supernovae Ia per unit of mass.
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_number_supernovae_ia(
    struct supernovae_ia *snia, float m1, float m2) {

#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  /* Do we have white dwarf? */
  if (m1 > snia->mass_max_progenitor) {
    return 0.;
  }

  float number_companion = 0.;
  for(int i = 0; i < NUMBER_TYPE_OF_COMPANION; i++) {
    /* Check if we are in the possible interval */
    if (m1 > snia->companion[i].mass_max ||
	m2 < snia->companion[i].mass_min)
      continue;

    /* Get mass limits */
    const float mass_min = max(m1, snia->companion[i].mass_min);
    const float mass_max = min(m2, snia->companion[i].mass_max);

    /* Compute number of companions */
    number_companion += stellar_evolution_get_companion_fraction(snia, mass_min, mass_max, i);
  }

  /* Use only the white dwarf already created */
  const float mass_min = max(m1, snia->mass_min_progenitor);

  /* Compute number of white dwarf */
  float number_white_dwarf = pow(snia->mass_max_progenitor, snia->progenitor_exponent);
  number_white_dwarf -= pow(mass_min, snia->progenitor_exponent);
  number_white_dwarf *= snia->progenitor_coef_exp;
  
  return number_companion * number_white_dwarf;
};

/**
 * @brief Compute the number of supernovae II per unit of mass (equation 3.47 in Poirier 2004).
 *
 * @param snii The #supernovae_ii model.
 * @param m1 The lower mass limit.
 * @param m2 The upper mass limit.
 *
 * @return The number of supernovae II per unit of mass.
 */
__attribute__((always_inline)) INLINE static float stellar_evolution_get_number_supernovae_ii(
    struct supernovae_ii *snii, float m1, float m2) {
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > m2)
    error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  /* Can we explode SNII? */
  if (m1 > snii->mass_max || m2 < snii->mass_min) {
    return 0.;
  }

  const float mass_min = max(m1, snii->mass_min);
  const float mass_max = min(m2, snii->mass_max);

  const float pow_mass = pow(mass_max, snii->exponent) - pow(mass_min, snii->exponent);

  return snii->coef_exp * pow_mass;
  
};

__attribute__((always_inline)) INLINE static float *stellar_evolution_get_supernovae_ia_yields(void) {
  return (float*)NULL;
};

__attribute__((always_inline)) INLINE static float *stellar_evolution_get_supernovae_ii_yields(void) {
  return (float*)NULL;
};

#endif // SWIFT_STELLAR_EVOLUTION_GEAR_H
