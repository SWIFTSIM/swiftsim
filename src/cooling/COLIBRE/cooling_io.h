/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_COLIBRE_IO_H
#define SWIFT_COOLING_COLIBRE_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cooling.h"
#include "io_properties.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of cooling to the file.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param cooling The #cooling_function_data
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, hid_t h_grp_columns,
    const struct cooling_function_data* cooling) {

  io_write_attribute_s(h_grp, "Cooling Model", "COLIBRE");
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

INLINE static void convert_part_sub_HI_frac(const struct engine* e,
                                            const struct part* p,
                                            const struct xpart* xp,
                                            float* ret) {

  ret[0] = cooling_get_subgrid_HI_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_HII_frac(const struct engine* e,
                                             const struct part* p,
                                             const struct xpart* xp,
                                             float* ret) {

  ret[0] = cooling_get_subgrid_HII_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_H2_frac(const struct engine* e,
                                            const struct part* p,
                                            const struct xpart* xp,
                                            float* ret) {

  ret[0] = cooling_get_subgrid_H2_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 * @param cooling The #cooling_function_data
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part* parts, const struct xpart* xparts, struct io_props* list,
    const struct cooling_function_data* cooling) {

  list[0] = io_make_output_field_convert_part(
      "Temperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts, xparts,
      convert_part_T, "Temperatures of the gas particles");

  list[1] = io_make_output_field_convert_part(
      "SubgridTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts,
      xparts, convert_part_sub_T,
      "The subgrid temperatures if the particles are within deltaT of the "
      "entropy floor the subgrid temperature is calculated assuming a "
      "pressure equilibrium on the entropy floor, if the particles are "
      "above deltaT of the entropy floor the subgrid temperature is "
      "identical to the SPH temperature.");

  list[2] = io_make_output_field_convert_part(
      "SubgridPhysicalDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, parts,
      xparts, convert_part_sub_rho,
      "The subgrid physical density if the particles are within deltaT of the "
      "entropy floor the subgrid density is calculated assuming a pressure "
      "equilibrium on the entropy floor, if the particles are above deltaT "
      "of the entropy floor the subgrid density is identical to the "
      "physical SPH density.");

  list[3] = io_make_output_field_convert_part(
      "HydrogenNeutralFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      xparts, convert_part_sub_HI_frac,
      "Fractions of neutral hydrogen atoms, nHI/nH, assuming equilibrium "
      "tables. If the particles are within deltaT of the entropy floor the "
      "fractions are calculated using the subgrid quantities, i.e. assuming a "
      "pressure equilibrium on the entropy floor. If the particles are "
      "above deltaT of the entropy floor, the normal hydro quantities are "
      "used.");

  list[4] = io_make_output_field_convert_part(
      "HydrogenIonizedFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      xparts, convert_part_sub_HII_frac,
      "Fractions of ionized hydrogen atoms, nHII/nH, assuming equilibrium "
      "tables. If the particles are within deltaT of the entropy floor the "
      "fractions are calculated using the subgrid quantities, i.e. assuming a "
      "pressure equilibrium on the entropy floor. If the particles are "
      "above deltaT of the entropy floor, the normal hydro quantities are "
      "used.");

  list[5] = io_make_output_field_convert_part(
      "HydrogenMolecularFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      xparts, convert_part_sub_H2_frac,
      "Fractions of hydrogen molecules, nH2/nH, assuming equilibrium tables. "
      "If the particles are within deltaT of the entropy floor the "
      "fractions are calculated using the subgrid quantities, i.e. assuming a "
      "pressure equilibrium on the entropy floor. If the particles are "
      "above deltaT of the entropy floor, the normal hydro quantities are "
      "used.");

  return 6;
}

#endif /* SWIFT_COOLING_COLIBRE_IO_H */
