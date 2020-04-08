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
#ifndef SWIFT_ENTROPY_FLOOR_COLIBRE_H
#define SWIFT_ENTROPY_FLOOR_COLIBRE_H

#include "adiabatic_index.h"
#include "cosmology.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "parser.h"
#include "units.h"

/**
 * @file src/entropy_floor/COLIBRE/entropy_floor.h
 * @brief Entropy floor used in the COLIBRE model
 */

/**
 * @brief Properties of the entropy floor in the COLIBRE model.
 */
struct entropy_floor_properties {

  /*! Density normalisation for the Jeans floor in Hydrogen atoms per cubic cm
   */
  float Jeans_density_norm_H_p_cm3;

  /*! Density normalisation for the Jeans floor in internal units */
  float Jeans_density_norm;

  /*! Inverse of the density normalisation for the Jeans floor in internal units
   */
  float Jeans_density_norm_inv;

  /*! Slope of the Jeans floor power-law */
  float Jeans_gamma_effective;

  /*! Temperature of the Jeans floor at the density normalisation in Kelvin */
  float Jeans_temperature_norm_K;

  /*! Temperature of the Jeans floor at the density normalisation in internal
   * units */
  float Jeans_temperature_norm;

  /*! Pressure of the Jeans floor at the density thresh. in internal units */
  float Jeans_pressure_norm;

  /*! Density normalisation for the Cool floor in Hydrogen atoms per cubic cm */
  float Cool_density_norm_H_p_cm3;

  /*! Density normalisation for the Cool floor in internal units */
  float Cool_density_norm;

  /*! Inverse of the density normalisation for the Cool floor in internal units
   */
  float Cool_density_norm_inv;

  /*! Slope of the Cool floor power-law */
  float Cool_gamma_effective;

  /*! Temperature of the Cool floor at the density normalisation in Kelvin */
  float Cool_temperature_norm_K;

  /*! Temperature of the Cool floor at the density normalisation. in internal
   * units */
  float Cool_temperature_norm;

  /*! Pressure of the Cool floor at the density normalisation in internal units
   */
  float Cool_pressure_norm;
};

/**
 * @brief Compute the entropy floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
static INLINE float entropy_floor(
    const struct part *p, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props) {

  /* Physical density in internal units */
  const float rho_phys = hydro_get_physical_density(p, cosmo);

  /* Physical pressure */
  float pressure = 0.f;

  /* Apply the Jeans limit. */
  const float pressure_Jeans = props->Jeans_pressure_norm *
                               powf(rho_phys * props->Jeans_density_norm_inv,
                                    props->Jeans_gamma_effective);

  pressure = max(pressure, pressure_Jeans);

  /* Only apply the Cool limit if
   * Cool_temperature_norm > 0. */
  if (props->Cool_temperature_norm > 0.) {
    const float pressure_Cool = props->Cool_pressure_norm *
                                powf(rho_phys * props->Cool_density_norm_inv,
                                     props->Cool_gamma_effective);

    pressure = max(pressure, pressure_Cool);
  }

  /* Convert to an entropy.
   * (Recall that the entropy is the same in co-moving and physical frames) */
  return gas_entropy_from_pressure(rho_phys, pressure);
}

/**
 * @brief Compute the temperature from the entropy floor for a given #part
 *
 * Calculate the EoS temperature, the particle is not updated.
 * This is the temperature exactly corresponding to the imposed EoS shape.
 * It only matches the entropy returned by the entropy_floor() function
 * for a neutral gas with primoridal abundance.
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
static INLINE float entropy_floor_temperature(
    const struct part *p, const struct cosmology *cosmo,
    const struct entropy_floor_properties *props) {

  /* Physical density in internal units */
  const float rho_phys = hydro_get_physical_density(p, cosmo);

  /* Physical */
  float temperature = 0.f;

  /* Apply the Jeans limit. */
  const float jeans_slope = props->Jeans_gamma_effective - 1.f;

  const float temperature_Jeans =
      props->Jeans_temperature_norm *
      pow(rho_phys * props->Jeans_density_norm_inv, jeans_slope);

  temperature = max(temperature, temperature_Jeans);

  /* Only apply the Cool limit if
   * Cool_temperature_norm > 0. */
  if (props->Cool_temperature_norm > 0.) {
    const float cool_slope = props->Cool_gamma_effective - 1.f;

    const float temperature_Cool =
        props->Cool_temperature_norm *
        pow(rho_phys * props->Cool_density_norm_inv, cool_slope);

    temperature = max(temperature, temperature_Cool);
  }

  return temperature;
}

/**
 * @brief Initialise the entropy floor by reading the parameters and converting
 * to internal units.
 *
 * The input temperatures and number densities are converted to entropy and
 * density assuming a neutral gas of primoridal abundance.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_const The physical constants.
 * @param hydro_props The propoerties of the hydro scheme.
 * @param props The entropy floor properties to fill.
 */
static INLINE void entropy_floor_init(struct entropy_floor_properties *props,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct hydro_props *hydro_props,
                                      struct swift_params *params) {

  /* Read the parameters in the units they are set */
  props->Jeans_density_norm_H_p_cm3 = parser_get_param_float(
      params, "COLIBREEntropyFloor:Jeans_density_norm_H_p_cm3");
  props->Jeans_temperature_norm_K = parser_get_param_float(
      params, "COLIBREEntropyFloor:Jeans_temperature_norm_K");
  props->Jeans_gamma_effective = parser_get_param_float(
      params, "COLIBREEntropyFloor:Jeans_gamma_effective");

  props->Cool_density_norm_H_p_cm3 = parser_get_param_float(
      params, "COLIBREEntropyFloor:Cool_density_norm_H_p_cm3");
  props->Cool_temperature_norm_K = parser_get_param_float(
      params, "COLIBREEntropyFloor:Cool_temperature_norm_K");
  props->Cool_gamma_effective = parser_get_param_float(
      params, "COLIBREEntropyFloor:Cool_gamma_effective");

  /* Initial Hydrogen abundance (mass fraction) */
  const double X_H = hydro_props->hydrogen_mass_fraction;

  /* Now convert to internal units assuming primodial Hydrogen abundance */
  props->Jeans_temperature_norm =
      props->Jeans_temperature_norm_K /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  props->Jeans_density_norm =
      props->Jeans_density_norm_H_p_cm3 /
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY) *
      phys_const->const_proton_mass / X_H;

  props->Cool_temperature_norm =
      props->Cool_temperature_norm_K /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  props->Cool_density_norm =
      props->Cool_density_norm_H_p_cm3 /
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY) *
      phys_const->const_proton_mass / X_H;

  /* We assume neutral gas */
  const float mean_molecular_weight = hydro_props->mu_neutral;

  /* Get the common terms */
  props->Jeans_density_norm_inv = 1.f / props->Jeans_density_norm;
  props->Cool_density_norm_inv = 1.f / props->Cool_density_norm;

  /* P_norm = (k_B * T) / (m_p * mu) * rho_norm */
  props->Jeans_pressure_norm =
      ((phys_const->const_boltzmann_k * props->Jeans_temperature_norm) /
       (phys_const->const_proton_mass * mean_molecular_weight)) *
      props->Jeans_density_norm;

  props->Cool_pressure_norm =
      ((phys_const->const_boltzmann_k * props->Cool_temperature_norm) /
       (phys_const->const_proton_mass * mean_molecular_weight)) *
      props->Cool_density_norm;
}

/**
 * @brief Print the properties of the entropy floor to stdout.
 *
 * @param props The entropy floor properties.
 */
static INLINE void entropy_floor_print(
    const struct entropy_floor_properties *props) {

  message("Entropy floor is 'COLIBRE' with:");
  message("Jeans limiter with slope n=%.3f at rho=%e (%e H/cm^3) and T=%.1f K",
          props->Jeans_gamma_effective, props->Jeans_density_norm,
          props->Jeans_density_norm_H_p_cm3, props->Jeans_temperature_norm);

  if (props->Cool_temperature_norm > 0.)
    message(
        " Cool limiter with slope n=%.3f at rho=%e (%e H/cm^3) and T=%.1f K",
        props->Cool_gamma_effective, props->Cool_density_norm,
        props->Cool_density_norm_H_p_cm3, props->Cool_temperature_norm);
  else
    message(" Cool limiter disabled");
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of entropy floor to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void entropy_floor_write_flavour(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Entropy floor", "COLIBRE");
}
#endif

/**
 * @brief Write an entropy floor struct to the given FILE as a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void entropy_floor_struct_dump(
    const struct entropy_floor_properties *props, FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct entropy_floor_properties),
                       1, stream, "entropy floor", "entropy floor properties");
}

/**
 * @brief Restore a entropy floor struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
static INLINE void entropy_floor_struct_restore(
    struct entropy_floor_properties *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct entropy_floor_properties), 1,
                      stream, NULL, "entropy floor properties");
}

#endif /* SWIFT_ENTROPY_FLOOR_COLIBRE_H */
