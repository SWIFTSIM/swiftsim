/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PRESSURE_FLOOR_GEAR_H
#define SWIFT_PRESSURE_FLOOR_GEAR_H

#include "adiabatic_index.h"
#include "cosmology.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "parser.h"
#include "units.h"

/**
 * @file src/pressure_floor/GEAR/pressure_floor.h
 * @brief Pressure floor used in the GEAR model
 */

/**
 * @brief Properties of the pressure floor in the GEAR model.
 */
struct pressure_floor_properties {

  /*! Jeans factor power 2/3 */
  float n_jeans_2_3;

  /*! Gravitational constant */
  double G;
  
};

/**
 * @brief Compute the physical pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param cosmo The #cosmology.
 * @param pressure The pressure without any pressure floor.
 */
static INLINE float pressure_floor_get_physical_pressure(
    const struct part *p, const struct cosmology *cosmo,
    const float pressure) {

  /* Physical density in internal units */
  const float rho = hydro_get_physical_density(p, cosmo);

  /* Compute pressure floor */
  float floor = 4.0 * M_1_PI * p->h * p->h * rho * pressure_floor_props.n_jeans_2_3;
  // TODO add sigma (will be done once the SF is merged)
  floor *= rho * hydro_one_over_gamma;

  return fmax(pressure, floor);
}

/**
 * @brief Compute the comoving pressure floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param pressure The pressure without any pressure floor.
 */
static INLINE float pressure_floor_get_comoving_pressure(
    const struct part *p, const float pressure) {

  /* Physical density in internal units */
  const float rho = hydro_get_comoving_density(p);

  /* Compute pressure floor */
  float floor = 4.0 * M_1_PI * p->h * p->h * rho * pressure_floor_props.n_jeans_2_3;
  // TODO add sigma (will be done once the SF is merged)
  floor *= rho * hydro_one_over_gamma;


  return fmax(pressure, floor);
}

/**
 * @brief Initialise the pressure floor by reading the parameters and converting
 * to internal units.
 *
 * The input temperatures and number densities are converted to pressure and
 * density assuming a neutral gas of primoridal abundance.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_const The physical constants.
 * @param hydro_props The propoerties of the hydro scheme.
 * @param props The pressure floor properties to fill.
 */
static INLINE void pressure_floor_init(struct pressure_floor_properties *props,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct hydro_props *hydro_props,
                                      struct swift_params *params) {

  /* Save the gravitational constant */
  props->G = phys_const->const_newton_G;

  /* Read the jeans factor */
  props->n_jeans_2_3 = parser_get_param_float(
      params, "GEARPressureFloor:Jeans_factor");

  props->n_jeans_2_3 = pow(props->n_jeans_2_3, 2./3.);

}

/**
 * @brief Print the properties of the pressure floor to stdout.
 *
 * @param props The pressure floor properties.
 */
static INLINE void pressure_floor_print(
    const struct pressure_floor_properties *props) {

  message("Pressure floor is 'GEAR' with:");
  message("Jeans factor: %g", pow(props->n_jeans_2_3, 3./2.));
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of pressure floor to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void pressure_floor_print_snapshot(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Pressure floor", "GEAR");
}
#endif

/**
 * @brief Write an pressure floor struct to the given FILE as a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void pressure_floor_struct_dump(
    const struct pressure_floor_properties *props, FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct pressure_floor_properties),
                       1, stream, "pressure floor", "pressure floor properties");
}

/**
 * @brief Restore a pressure floor struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void pressure_floor_struct_restore(
    struct pressure_floor_properties *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct pressure_floor_properties), 1,
                      stream, NULL, "pressure floor properties");
}

#endif /* SWIFT_PRESSURE_FLOOR_GEAR_H */
