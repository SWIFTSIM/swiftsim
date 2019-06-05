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
#ifndef SWIFT_STELLAR_EVOLUTION_IO_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_IO_GEAR_H

#include "stellar_evolution_struct.h"
#include "stellar_evolution.h"


/**
 * @brief Compute the coefficients of the initial mass function.
 *
 * @param imf The #initial_mass_function.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_compute_initial_mass_function_coefficients(
    struct initial_mass_function *imf) {
  /* Allocate memory */
  if ((imf->coef = (float *)malloc(sizeof(float) * imf->n_parts)) ==
      NULL)
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

  /* Write the masses in the correct attributes */
  imf->mass_min = imf->mass_limits[0];
  imf->mass_max = imf->mass_limits[imf->n_parts];

  /* Compute the coefficients */
  stellar_evolution_compute_initial_mass_function_coefficients(imf);

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

  /* Read constant terms */
  parser_get_param_float_array(params, "GEARLifetime:constant", 3, lt->constant);
  
}

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

__attribute__((always_inline)) INLINE static void stellar_evolution_init_supernovae_ii(
    struct supernovae_ii *snii, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct initial_mass_function *imf) {

  /* Read the minimal mass of a supernovae */
  snii->mass_min = parser_get_param_float(params, "GEARSupernovaeII:min_mass");

  /* Read the maximal mass of a supernovae */
  snii->mass_max = parser_get_param_float(params, "GEARSupernovaeII:max_mass");

  /* Get the IMF parameters */
  snii->exponent = stellar_evolution_get_imf_exponent(imf, snii->mass_min, snii->mass_max);
  snii->coef_exp = stellar_evolution_get_imf_coefficient(imf, snii->mass_min, snii->mass_max);
  snii->coef_exp /= snii->exponent;
}


/**
 * @brief Initialize the global properties of the stellar evolution scheme.
 *
 * @param sm The #stellar_model.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_props_init(
    struct stellar_model* sm, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct cosmology* cosmo) {

  /* Initialize the initial mass function */
  stellar_evolution_init_initial_mass_function(&sm->imf, phys_const, us, params);

  /* Initialize the lifetime model */
  stellar_evolution_init_lifetime(&sm->lifetime, phys_const, us, params);

  /* Initialize the supernovae Ia model */
  stellar_evolution_init_supernovae_ia(&sm->snia, phys_const, us, params, &sm->imf);
 
  /* Initialize the supernovae II model */
  stellar_evolution_init_supernovae_ii(&sm->snii, phys_const, us, params, &sm->imf);

}

#endif // SWIFT_STELLAR_EVOLUTION_IO_GEAR_H
