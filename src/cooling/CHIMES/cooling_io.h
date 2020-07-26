/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_CHIMES_IO_H
#define SWIFT_COOLING_CHIMES_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cooling.h"
#include "engine.h"
#include "io_properties.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param cooling the parameters of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, hid_t h_grp_columns,
    const struct cooling_function_data* cooling) {

  io_write_attribute_s(h_grp, "Cooling Model", "CHIMES");

  /* Add the species names to the named columns */
  hsize_t dims[1] = {CHIMES_NETWORK_SIZE};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, CHIMES_NAME_STR_LENGTH);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "SpeciesFractions", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           cooling->chimes_species_names_reduced[0]);
  H5Dclose(dset);

  H5Tclose(type);
  H5Sclose(space);
}
#endif

INLINE static void convert_part_T(const struct engine* e, const struct part* p,
                                  const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_temperature(e->physical_constants, e->hydro_properties,
                                   e->internal_units, e->cosmology,
                                   e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_T(const struct engine* e,
                                      const struct part* p,
                                      const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_subgrid_temperature(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_rho(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_subgrid_density(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

INLINE static void convert_part_chimes_abundances(const struct engine* e,
                                                  const struct part* p,
                                                  const struct xpart* xp,
                                                  float* ret) {
  /* Create dummy part and xpart
   * structures that we can use to
   * temporarily update the CHIMES
   * abundance array for the snapshot. */
  struct part dummy_p = *p;
  struct xpart dummy_xp = *xp;

  /* Update abundance array for
   * particles on the EOS. */
  cooling_set_subgrid_properties(e->physical_constants, e->internal_units,
                                 e->cosmology, e->hydro_properties,
                                 e->entropy_floor, e->cooling_func, 
				 e->dustevo, &dummy_p, &dummy_xp);

  int i;
  for (i = 0; i < CHIMES_NETWORK_SIZE; i++)
    ret[i] = (float)dummy_xp.cooling_data.chimes_abundances[i];
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @param cooling The #cooling_function_data
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part* parts, const struct xpart* xparts, struct io_props* list,
    const struct cooling_function_data* cooling) {

  list[0] = io_make_output_field_convert_part(
      "SpeciesFractions", FLOAT, CHIMES_NETWORK_SIZE, UNIT_CONV_NO_UNITS, 0.f,
      parts, xparts, convert_part_chimes_abundances,
      "Species fractions array for all ions and molecules in the CHIMES "
      "network. The fraction of species i is defined in terms "
      "of its number density relative to hydrogen, i.e. n_i / n_H_tot.");

  list[1] = io_make_output_field_convert_part(
      "Temperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts, xparts,
      convert_part_T, "Temperature of the particles");

  list[2] = io_make_output_field_convert_part(
      "SubgridTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts,
      xparts, convert_part_sub_T,
      "The subgrid temperatures if the particles are within deltaT of the "
      "entropy floor the subgrid temperature is calculated assuming a "
      "pressure equilibrium on the entropy floor, if the particles are "
      "above deltaT of the entropy floor the subgrid temperature is "
      "identical to the SPH temperature.");

  list[3] = io_make_output_field_convert_part(
      "SubgridPhysicalDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, parts,
      xparts, convert_part_sub_rho,
      "The subgrid physical density if the particles are within deltaT of the "
      "entropy floor the subgrid density is calculated assuming a pressure "
      "equilibrium on the entropy floor, if the particles are above deltaT "
      "of the entropy floor the subgrid density is identical to the "
      "physical SPH density.");

  return 4;
}

#endif /* SWIFT_COOLING_CHIMES_IO_H */
