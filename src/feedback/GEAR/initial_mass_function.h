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
#ifndef SWIFT_INITIAL_MASS_FUNCTION_GEAR_H
#define SWIFT_INITIAL_MASS_FUNCTION_GEAR_H

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
  error("This has not been tested. Need to check the units");
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
  error("This has not been tested. Need to check the units");
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
 * @brief Compute the coefficients of the initial mass function.
 *
 * @param imf The #initial_mass_function.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_compute_initial_mass_function_coefficients(
    struct initial_mass_function *imf) {

  /* Allocate memory */
  if ((imf->coef = (float *) malloc(sizeof(float) * imf->n_parts)) == NULL)
    error("Failed to allocate the IMF coefficients.");

  /* Suppose that the first coefficients is 1 (will be corrected later) */
  imf->coef[0] = 1.;

  /* Use the criterion of continuity for the IMF */
  for(int i = 1; i < imf->n_parts; i++) {
    float exp = imf->exp[i-1] - imf->exp[i];
    imf->coef[i] = imf->coef[i-1] * pow(imf->mass_limits[i], exp);
  }

  /* Use the criterion on the integral = 1 */
  float integral = 0;
  for(int i = 0; i < imf->n_parts; i++) {
    const float exp = imf->exp[i] + 1.;
    const float m_i = pow(imf->mass_limits[i], exp);
    const float m_i1 = pow(imf->mass_limits[i+1], exp);
    integral += imf->coef[i] * (m_i1 - m_i)/ exp;
  }

  /* Normalize the coefficients (fix initial supposition) */
  for(int i = 0; i < imf->n_parts; i++) {
    imf->coef[i] /= integral;
  }

}

/**
 * @brief Initialize the initial mass function.
 *
 * @param imf The #initial_mass_function.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The #swift_params.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_init_initial_mass_function(
    struct initial_mass_function* imf, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

  /* Read the number of elements */
  imf->n_parts = parser_get_param_int(params, "GEARInitialMassFunction:number_function_part");

  /* Read the exponents */
  if ((imf->exp = (float *)malloc(sizeof(float) * imf->n_parts)) ==
      NULL)
    error("Failed to allocate the IMF exponents.");

  parser_get_param_float_array(params, "GEARInitialMassFunction:exponents", imf->n_parts, imf->exp);

  /* Read the mass limits */
  if ((imf->mass_limits = (float *)malloc(sizeof(float) * (imf->n_parts + 1))) ==
      NULL)
    error("Failed to allocate the IMF temporary masses.");

  parser_get_param_float_array(params, "GEARInitialMassFunction:mass_limits_msun",
			       imf->n_parts + 1, imf->mass_limits);

  /* change the mass limits to internal units */
  for(int i = 0; i < imf->n_parts + 1; i++) {
    imf->mass_limits[i] *= phys_const->const_solar_mass;
  }

  /* Write the masses in the correct attributes */
  imf->mass_min = imf->mass_limits[0];
  imf->mass_max = imf->mass_limits[imf->n_parts];

  /* Compute the coefficients */
  stellar_evolution_compute_initial_mass_function_coefficients(imf);

}

#endif // SWIFT_INITIAL_MASS_FUNCTION_GEAR_H
