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
#ifndef SWIFT_SUPERNOVAE_IA_GEAR_H
#define SWIFT_SUPERNOVAE_IA_GEAR_H


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

__attribute__((always_inline)) INLINE static float *stellar_evolution_get_supernovae_ia_yields(void) {
  return (float*)NULL;
};

/**
 * @brief Initialize the companion structure in the #supernovae_ia.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_init_supernovae_ia_companion(
    struct supernovae_ia *snia) {

  for(int i = 0; i < NUMBER_TYPE_OF_COMPANION; i++) {
    /* Compute the integral */
    float integral = stellar_evolution_get_companion_fraction(
        snia, snia->companion[i].mass_min, snia->companion[i].mass_max, i);

    /* Update the coefficient for a normalization to 1 of the IMF */
    snia->companion[i].coef *= snia->companion[i].coef / integral;
  }
  
}

/**
 * @brief Initialize the #supernovae_ia structure.
 *
 * @param snia The #supernovae_ia model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The simulation parameters.
 * @param imf The #initial_mass_function model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_init_supernovae_ia(
    struct supernovae_ia *snia, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct initial_mass_function *imf) {

  /* Read the exponent of the IMF for companion */
  snia->companion_exponent = parser_get_param_float(params, "GEARSupernovaeIa:exponent");
  
  /* Read the minimal mass for a white dwarf */
  snia->mass_min_progenitor = parser_get_param_float(params, "GEARSupernovaeIa:min_mass_white_dwarf_progenitor");

  /* Read the maximal mass for a white dwarf */
  snia->mass_max_progenitor = parser_get_param_float(params, "GEARSupernovaeIa:max_mass_white_dwarf_progenitor");

  /* Get the IMF parameters */
  snia->progenitor_exponent = stellar_evolution_get_imf_exponent(imf,
      snia->mass_min_progenitor, snia->mass_max_progenitor);
  snia->progenitor_coef_exp = stellar_evolution_get_imf_coefficient(imf,
      snia->mass_min_progenitor, snia->mass_max_progenitor);
  snia->progenitor_coef_exp /= snia->progenitor_exponent;

  /* Read the maximal mass of a red giant companion */
  snia->companion[0].mass_max = parser_get_param_float(params, "GEARSupernovaeIa:max_mass_red_giant");

  /* Read the minimal mass of a red giant companion */
  snia->companion[0].mass_min = parser_get_param_float(params, "GEARSupernovaeIa:min_mass_red_giant");

  /* Read the coefficient of the main sequence companion */
  snia->companion[0].coef = parser_get_param_float(params, "GEARSupernovaeIa:coef_red_giant");

  /* Read the maximal mass of a main sequence companion */
  snia->companion[1].mass_max = parser_get_param_float(params, "GEARSupernovaeIa:max_mass_main_sequence");

  /* Read the minimal mass of a main sequence companion */
  snia->companion[1].mass_min = parser_get_param_float(params, "GEARSupernovaeIa:min_mass_main_sequence");

  /* Read the coefficient of the main sequence companion */
  snia->companion[1].coef = parser_get_param_float(params, "GEARSupernovaeIa:coef_main_sequence");

  /* Compute the normalization coefficients of the companion IMF */
  stellar_evolution_init_supernovae_ia_companion(snia);
  
}

#endif // SWIFT_SUPERNOVAE_IA_GEAR_H
