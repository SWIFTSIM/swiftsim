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
#ifndef SWIFT_SUPERNOVAE_II_GEAR_H
#define SWIFT_SUPERNOVAE_II_GEAR_H

#include "hdf5_functions.h"
#include "interpolation.h"
#include "stellar_evolution_struct.h"

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
    const struct supernovae_ii *snii, float m1, float m2) {
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

/**
 * @brief Get the SNII yields.
 *
 * @param snii The #supernovae_ii model.
 * @param m1 The lower mass.
 * @param m2 The upper mass.
 * @param yields The elements ejected.
 * @param mass_ejected The mass of non processsed elements.
 * @param mass_ejected_processed The mass of processed elements.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_get_supernovae_ii_yields(
    const struct supernovae_ii *snii, float m1, float m2,
    float *yields, float *mass_ejected, float *mass_ejected_processed) {

  error("TODO");
};


/**
 * @brief Read the SNII yields from the table.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_read_supernovae_ii_yields(
    struct supernovae_ii *snii, struct swift_params* params, const struct stellar_model *sm) {

  hid_t file_id, group_id;

  hsize_t previous_count = 0;

  /* Open IMF group */
  h5_open_group(params, "Data/SNII", &file_id, &group_id);

  /* Do all the elements */
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {

    /* Get the element name */
    const char *name = stellar_evolution_get_element_name(sm, i);

    
    /* Now let's get the number of elements */
    /* Open attribute */
    const hid_t h_dataset = H5Dopen(group_id, name, H5P_DEFAULT);
    if (h_dataset < 0) error("Error while opening attribute '%s'", name);

    /* Get the number of elements */
    hsize_t count = io_get_number_element_in_dataset(h_dataset);

    /* Check that all the arrays have the same size */
    if (i != 0 && count != previous_count) {
      error("The code is not able to deal with yields arrays of different size");
    }
    previous_count = count;

    /* Close the attribute */
    H5Dclose(h_dataset);

    /* Allocate the memory */
    float *data = (float*) malloc(sizeof(float) * count);
    if (data == NULL)
      error("Failed to allocate the SNII yields for %s.", name);

    /* Read the dataset */
    io_read_array_dataset(group_id, name, FLOAT,
			  data, count);

    /* Initialize the interpolation */
    interpolate_1d_init(&snii->yields[i], log10(snii->mass_min), log10(snii->mass_max), data, count, /* Integrate? */0);

    /* Cleanup the memory */
    free(data);
  }

  /* Read the mass ejected */

  /* Allocate the memory */
  float *mass_ejected = (float*) malloc(sizeof(float) * previous_count);
  if (mass_ejected == NULL)
    error("Failed to allocate the SNII yields for the mass ejected.");

  /* Read the dataset */
  io_read_array_dataset(group_id, "Ej", FLOAT,
			mass_ejected, previous_count);

  /* Initialize the interpolation */
  interpolate_1d_init(&snii->ejected_mass_processed, log10(snii->mass_min), log10(snii->mass_max), mass_ejected, previous_count, /* Integrate? */0);

  /* Read the mass ejected of non processed gas */

  /* Allocate the memory */
  float *mass_ejected_non_process = (float*) malloc(sizeof(float) * previous_count);
  if (mass_ejected_non_process == NULL)
    error("Failed to allocate the SNII yields for the (non processed) mass ejected.");

  /* Read the dataset */
  io_read_array_dataset(group_id, "Ejnp", FLOAT,
			mass_ejected_non_process, previous_count);

  /* Initialize the interpolation */
  interpolate_1d_init(&snii->ejected_mass, log10(snii->mass_min), log10(snii->mass_max), mass_ejected_non_process, previous_count, /* Integrate? */0);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);  
};


/**
 * @brief Reads the supernovae II parameters from parameters file.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_read_supernovae_ii_from_params(
    struct supernovae_ii *snii, struct swift_params* params) {

  /* Read the minimal mass of a supernovae */
  snii->mass_min = parser_get_opt_param_float(params, "GEARSupernovaeII:min_mass", snii->mass_min);

  /* Read the maximal mass of a supernovae */
  snii->mass_max = parser_get_opt_param_float(params, "GEARSupernovaeII:max_mass", snii->mass_max);

}

/**
 * @brief Reads the supernovae II parameters from tables.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_read_supernovae_ii_from_tables(
    struct supernovae_ii *snii, struct swift_params* params) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(params, "Data/SNII", &file_id, &group_id);

  /* Read the minimal mass of a supernovae */
  io_read_attribute(group_id, "Mmin", FLOAT, &snii->mass_min);

  /* Read the maximal mass of a supernovae */
  io_read_attribute(group_id, "Mmax", FLOAT, &snii->mass_max);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
}


/**
 * @brief Initialize the #supernovae_ii structure.
 *
 * @param snii The #supernovae_ii model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_init_supernovae_ii(
    struct supernovae_ii *snii, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct stellar_model *sm) {

  /* Read the parameters from the tables */
  stellar_evolution_read_supernovae_ii_from_tables(snii, params);
  
  /* Read the parameters from the params file */
  stellar_evolution_read_supernovae_ii_from_tables(snii, params);

  /* Read the supernovae yields */
  stellar_evolution_read_supernovae_ii_yields(snii, params, sm);
  
  /* Apply the unit changes */
  snii->mass_min *= phys_const->const_solar_mass;
  snii->mass_max *= phys_const->const_solar_mass;

  /* Apply the unit changes to the data */
  for(int i = 0; i < CHEMISTRY_ELEMENT_COUNT; i++) {
    interpolate_1d_change_units(&snii->yields[i], phys_const->const_solar_mass,
				phys_const->const_solar_mass);
    // TODO multiply by the mass ejected
  }

  /* Get the IMF parameters */
  snii->exponent = stellar_evolution_get_imf_exponent(&sm->imf, snii->mass_min, snii->mass_max);
  snii->coef_exp = stellar_evolution_get_imf_coefficient(&sm->imf, snii->mass_min, snii->mass_max);
  snii->coef_exp /= snii->exponent;
}


#endif // SWIFT_SUPERNOVAE_II_GEAR_H
