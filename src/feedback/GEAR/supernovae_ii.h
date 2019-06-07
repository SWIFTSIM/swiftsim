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

__attribute__((always_inline)) INLINE static float *stellar_evolution_get_supernovae_ii_yields(void) {
  return (float*)NULL;
};


/**
 * @brief Reads the supernovae II parameters from parameters file.
 *
 * @param snii The #supernovae_ii model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The simulation parameters.
 * @param imf The #initial_mass_function model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_read_supernovae_ii_from_params(
    struct supernovae_ii *snii, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

  /* Read the minimal mass of a supernovae */
  snii->mass_min = parser_get_opt_param_float(params, "GEARSupernovaeII:min_mass", snii->mass_min);

  /* Read the maximal mass of a supernovae */
  snii->mass_max = parser_get_opt_param_float(params, "GEARSupernovaeII:max_mass", snii->mass_max);

}

/**
 * @brief Reads the supernovae II parameters from tables.
 *
 * @param snii The #supernovae_ii model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The simulation parameters.
 * @param imf The #initial_mass_function model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_read_supernovae_ii_from_tables(
    struct supernovae_ii *snii, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

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
 * @param imf The #initial_mass_function model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_init_supernovae_ii(
    struct supernovae_ii *snii, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct initial_mass_function *imf) {

  /* Read the parameters from the tables */
  stellar_evolution_read_supernovae_ii_from_tables(snii, phys_const, us, params);
  
  /* Read the parameters from the params file */
  stellar_evolution_read_supernovae_ii_from_tables(snii, phys_const, us, params);
  
  /* Apply the unit changes */
  snii->mass_min *= phys_const->const_solar_mass;
  snii->mass_max *= phys_const->const_solar_mass;

  /* Get the IMF parameters */
  snii->exponent = stellar_evolution_get_imf_exponent(imf, snii->mass_min, snii->mass_max);
  snii->coef_exp = stellar_evolution_get_imf_coefficient(imf, snii->mass_min, snii->mass_max);
  snii->coef_exp /= snii->exponent;
}


#endif // SWIFT_SUPERNOVAE_II_GEAR_H
