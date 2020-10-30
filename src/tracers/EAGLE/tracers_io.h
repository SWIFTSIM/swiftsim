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
#ifndef SWIFT_TRACERS_EAGLE_IO_H
#define SWIFT_TRACERS_EAGLE_IO_H

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

  io_write_attribute_s(h_grp, "Tracers", "EAGLE");
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
        "MaximalTemperatureTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, xparts,
        tracers_data.maximum_temperature_time,
        "Times at which the maximal temperature was reached");
  }

  list[2] =
      io_make_output_field(
        "HeatedBySNIIFeedback", INT, 1, UNIT_CONV_NO_UNITS,
        0.f, xparts, tracers_data.hit_by_SNII_feedback,
        "Number of times that the particles have been directly hit by a SNII "
        "feedback event at some point in the past.");

  list[3] =
      io_make_output_field(
        "HeatedByAGNFeedback", INT, 1, UNIT_CONV_NO_UNITS,
        0.f, xparts, tracers_data.hit_by_AGN_feedback,
        "Number of times that the particles have been directly hit by an AGN "
        "feedback event at some point in the past.");

  list[4] = io_make_output_field("EnergiesReceivedFromAGNFeedback", FLOAT, 1,
                                 UNIT_CONV_ENERGY, 0.f, xparts,
                                 tracers_data.AGN_feedback_energy,
                                 "Total amount of thermal energy from AGN "
                                 "feedback events received by the particles.");

  list[5] = io_make_output_field("EnergiesReceivedFromSNIIFeedback", FLOAT, 1,
                                 UNIT_CONV_ENERGY, 0.f, xparts,
                                 tracers_data.SNII_feedback_energy,
                                 "Total amount of thermal energy from SNII "
                                 "feedback events received by the particles.");

  list[6] = io_make_output_field(
    "AverageSNIIBirthDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, xparts,
    tracers_data.SNII_birth_density,
    "Average birth densities of the stars whose SNII feedback events have "
    "heated each particle. Set to -FLT_MAX for particles that have never "
    "been hit by SNII feedback.");

  list[7] = io_make_output_field(
    "MaximumStarFormingDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, xparts,
    tracers_data.maximum_physical_density,
    "Maximum physical densities of each gas particle while having a non-zero "
    "star formation rate.");


  return 8;
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
        "MaximalTemperatureTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.maximum_temperature_time,
        "Times at which the maximal temperature was reached");
  }

  list[2] =
      io_make_output_field(
        "HeatedBySNIIFeedback", INT, 1, UNIT_CONV_NO_UNITS,
        0.f, sparts, tracers_data.hit_by_SNII_feedback,
        "Number of times that the particles have been directly hit by a SNII "
        "feedback event at some point in the past.");

  list[3] =
      io_make_output_field(
        "HeatedByAGNFeedback", INT, 1, UNIT_CONV_NO_UNITS,
        0.f, sparts, tracers_data.hit_by_AGN_feedback,
        "Number of times that the particles have been directly hit by an AGN "
        "feedback event at some point in the past.");

  list[4] = io_make_output_field(
      "EnergiesReceivedFromAGNFeedback", FLOAT, 1, UNIT_CONV_ENERGY, 0.f,
      sparts, tracers_data.AGN_feedback_energy,
      "Total amount of thermal energy from AGN feedback events received by the "
      "particles when the particle was still a gas particle.");

  list[5] = io_make_output_field("EnergiesReceivedFromSNIIFeedback", FLOAT, 1,
                                 UNIT_CONV_ENERGY, 0.f, sparts,
                                 tracers_data.SNII_feedback_energy,
                                 "Total amount of thermal energy from SNII "
                                 "feedback events received by the particles.");

  list[6] = io_make_output_field(
    "AverageSNIIBirthDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
    tracers_data.SNII_birth_density,
    "Average birth densities of the stars whose SNII feedback events have "
    "heated each particle. Set to -FLT_MAX for particles that have never "
    "been hit by SNII feedback.");

  list[7] = io_make_output_field(
    "MaximumStarFormingDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
    tracers_data.maximum_physical_density,
    "Maximum physical densities of each progenitor gas particle while having a "
    "non-zero star formation rate.");

  return 8;
}

#endif /* SWIFT_TRACERS_EAGLE_IO_H */
