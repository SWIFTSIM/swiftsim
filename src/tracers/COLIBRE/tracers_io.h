/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TRACERS_COLIBRE_IO_H
#define SWIFT_TRACERS_COLIBRE_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "io_properties.h"
#include "tracers.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of tracers to the file.
 *
 * @param h_grp The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void tracers_write_flavour(
    hid_t h_grp) {

  io_write_attribute_s(h_grp, "Tracers", "COLIBRE");
}
#endif

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int tracers_write_particles(
    const struct part* parts, const struct xpart* xparts, struct io_props* list,
    const int with_cosmology) {

  list[0] = io_make_output_field(
      "MaximalTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, xparts,
      tracers_data.maximum_temperature,
      "Maximal temperatures ever reached by the particles");

  if (with_cosmology) {
    list[1] = io_make_output_field(
        "MaximalTemperatureScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        xparts, tracers_data.maximum_temperature_scale_factor,
        "Scale-factors at which the maximal temperature was reached");

  } else {

    list[1] = io_make_output_field(
        "MaximalTemperatureTimes", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
        tracers_data.maximum_temperature_time,
        "Times at which the maximal temperature was reached");
  }

  list[2] = io_make_output_field(
      "StellarMomentaReceived", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f, xparts,
      tracers_data.momentum_received,
      "Momentum received from stellar winds in physical coordinates");

  list[3] = io_make_output_field(
      "HIIregionsEndTime", FLOAT, 1, UNIT_CONV_TIME, 0.f, xparts,
      tracers_data.HIIregion_timer_gas,
      "Time until particle is in HII region");

  list[4] = io_make_output_field(
      "HIIregionsStarID", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      tracers_data.HIIregion_starid,
      "ID of star particle responsible for this HII region");

  list[5] = io_make_output_field(
      "HydrogenNeutralFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      tracers_data.nHI_over_nH,
      "Fraction of neutral hydrogen atoms, nHI/nH");

  list[6] = io_make_output_field(
      "HydrogenIonizedFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      tracers_data.nHII_over_nH,
      "Fraction of ionized hydrogen atoms, nHII/nH");

  list[7] = io_make_output_field(
      "HydrogenMolecularFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
      tracers_data.nH2_over_nH,
      "Fraction of hydrogen molecules, nH2/nH");

  list[8] = io_make_output_field(
      "SubgridDensity", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, xparts,
      tracers_data.subgrid_dens,
      "Subgrid density");

  list[9] = io_make_output_field(
      "SubgridTemperature", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, xparts,
      tracers_data.subgrid_temp,
      "Subgrid temperature");

  return 10;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart* sparts, struct io_props* list,
    const int with_cosmology) {

  list[0] = io_make_output_field(
      "MaximalTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, sparts,
      tracers_data.maximum_temperature,
      "Maximal temperatures ever reached by the particles before they got "
      "converted to stars");

  if (with_cosmology) {
    list[1] = io_make_output_field(
        "MaximalTemperatureScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        sparts, tracers_data.maximum_temperature_scale_factor,
        "Scale-factors at which the maximal temperature was reached");

  } else {

    list[1] = io_make_output_field(
        "MaximalTemperatureTimes", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.maximum_temperature_time,
        "Times at which the maximal temperature was reached");
  }

  list[2] = io_make_output_field(
      "StellarMomentaReceived", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f, sparts,
      tracers_data.momentum_received,
      "Momentum received by the gas particles from stellar winds before it was "
      "converted to a star in physical coordinates");

  return 3;
}

#endif /* SWIFT_TRACERS_COLIBRE_IO_H */
