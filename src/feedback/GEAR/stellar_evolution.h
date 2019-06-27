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
#ifndef SWIFT_STELLAR_EVOLUTION_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_GEAR_H

#include "hdf5_functions.h"
#include "initial_mass_function.h"
#include "lifetime.h"
#include "random.h"
#include "stellar_evolution_struct.h"
#include "supernovae_ia.h"
#include "supernovae_ii.h"

#include <math.h>
#include <stddef.h>


/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function compute the SN rate and yields before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_evolve_spart(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct cosmology* cosmo, const struct unit_system* us,
    const integertime_t ti_begin,
    const double star_age_beg_step, const double dt) {

  /* Get the metallicity */
  const float metallicity = chemistry_spart_metal_mass_fraction(sp);

  /* Compute masses range */
  const float log_m_beg_step = lifetime_get_log_mass_from_lifetime(
    &sm->lifetime, log10(star_age_beg_step), metallicity);
  const float log_m_end_step = lifetime_get_log_mass_from_lifetime(
    &sm->lifetime, log10(star_age_beg_step + dt), metallicity);

  const float m_beg_step = pow(10, log_m_beg_step);
  const float m_end_step = pow(10, log_m_end_step);

  /* Check if the star can produce a supernovae */
  const int can_produce_snia = supernovae_ia_can_explode(&sm->snia, m_beg_step) ||
    supernovae_ia_can_explode(&sm->snia, m_beg_step);
  const int can_produce_snii = supernovae_ii_can_explode(&sm->snii, m_beg_step) ||
    supernovae_ii_can_explode(&sm->snii, m_end_step);

  if (can_produce_snia) {
    /* Compute rates */
    const float number_snia_f = supernovae_ia_get_number(
      &sm->snia, m_end_step, m_beg_step);

    /* Get the number of SN */
    const float rand_snia = random_unit_interval(sp->id, ti_begin,
						 random_number_stellar_feedback);

    /* Get the fraction part */
    const float frac_snia = number_snia_f - floor(number_snia_f);

    /* Get the integer number of SN */
    sp->feedback_data.number_snia = floor(number_snia_f) + (rand_snia < frac_snia)? 1 : 0;

  }

  if (can_produce_snii) {
    /* Compute rates */
    const float number_snii_f = supernovae_ii_get_number(
    &sm->snii, m_end_step, m_beg_step);

    /* Get the number of SN */
    const float rand_snii = random_unit_interval(sp->id, ti_begin,
					       random_number_stellar_feedback_2);

    /* Get the fraction part */
    const float frac_snii = number_snii_f - floor(number_snii_f);
    
    /* Get the integer number of SN */
    sp->feedback_data.number_snii = floor(number_snii_f) + (rand_snii < frac_snii)? 1 : 0;
  }

  if (sp->feedback_data.number_snia != 0 || sp->feedback_data.number_snii != 0) {
    /* Decrease star mass by amount of mass distributed to gas neighbours */
    // TODO
    float mass_ejected = 0.f;
    if (sp->mass < mass_ejected)
      error("Stars cannot have negative mass. (%g < %g)",
	    sp->mass, mass_ejected);
    sp->mass -= mass_ejected;
  }
}

/**
 * @brief Get the name of the element i.
 *
 * @param sm The #stellar_model.
 * @param i The element indice.
 */
__attribute__((always_inline)) INLINE static const char* stellar_evolution_get_element_name(
    const struct stellar_model *sm, int i) {

  return sm->elements_name + i * GEAR_LABELS_SIZE;
}

/**
 * @brief Read the name of all the elements present in the tables.
 *
 * @param sm The #stellar_model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_read_elements(
    struct stellar_model* sm, struct swift_params* params) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(params, "Data", &file_id, &group_id);

  /* Read the elements */
  io_read_string_array_attribute(
    group_id, "elts", sm->elements_name, CHEMISTRY_ELEMENT_COUNT,
    GEAR_LABELS_SIZE);

  /* Print the name of the elements */
  char txt[CHEMISTRY_ELEMENT_COUNT * (GEAR_LABELS_SIZE + 2)] = "";
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    if (i != 0) {
      strcat(txt, ", ");
    }
    strcat(txt, stellar_evolution_get_element_name(sm, i));
  }

  message("Chemistry elements: %s", txt);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);

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

  /* Read the list of elements */
  stellar_evolution_read_elements(sm, params);

  /* Initialize the initial mass function */
  initial_mass_function_init(&sm->imf, phys_const, us, params);

  /* Initialize the lifetime model */
  lifetime_init(&sm->lifetime, phys_const, us, params);

  /* Initialize the supernovae Ia model */
  supernovae_ia_init(&sm->snia, phys_const, us, params, sm);
 
  /* Initialize the supernovae II model */
  supernovae_ii_init(&sm->snii, phys_const, us, params, sm);


}

#endif // SWIFT_STELLAR_EVOLUTION_GEAR_H
