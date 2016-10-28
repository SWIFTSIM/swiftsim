/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_UNITS_H
#define SWIFT_UNITS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "parser.h"

/**
 * @brief The unit system used internally.
 *
 * This structure contains the conversion factors to the 7 cgs base units to the
 * internal units. It is used everytime a conversion is performed or an i/o
 * function is called.
 **/
struct UnitSystem {

  /*! Conversion factor from grams to internal mass units */
  double UnitMass_in_cgs;

  /*! Conversion factor from centimeters to internal length unit */
  double UnitLength_in_cgs;

  /*! Conversion factor from seconds to internal time units */
  double UnitTime_in_cgs;

  /*! Conversion factor from Ampere to internal current units */
  double UnitCurrent_in_cgs;

  /*! Conversion factor from Kelvins to internal temperature units. */
  double UnitTemperature_in_cgs;
};

/**
 * @brief The base units used in the cgs (and internal) system. All units are
 * derived from those.
 */
enum BaseUnits {
  UNIT_MASS = 0,
  UNIT_LENGTH = 1,
  UNIT_TIME = 2,
  UNIT_CURRENT = 3,
  UNIT_TEMPERATURE = 4
};

/**
 * @brief  The different conversion factors supported by default
 */
enum UnitConversionFactor {
  UNIT_CONV_NO_UNITS,
  UNIT_CONV_MASS,
  UNIT_CONV_LENGTH,
  UNIT_CONV_TIME,
  UNIT_CONV_DENSITY,
  UNIT_CONV_SPEED,
  UNIT_CONV_ACCELERATION,
  UNIT_CONV_FORCE,
  UNIT_CONV_ENERGY,
  UNIT_CONV_ENERGY_PER_UNIT_MASS,
  UNIT_CONV_ENTROPY,
  UNIT_CONV_ENTROPY_PER_UNIT_MASS,
  UNIT_CONV_POWER,
  UNIT_CONV_PRESSURE,
  UNIT_CONV_FREQUENCY,
  UNIT_CONV_ELECTRIC_CHARGE,
  UNIT_CONV_ELECTRIC_VOLTAGE,
  UNIT_CONV_ELECTRIC_CAPACITANCE,
  UNIT_CONV_ELECTRIC_RESISTANCE,
  UNIT_CONV_ELECTRIC_CONDUCTANCE,
  UNIT_CONV_MAGNETIC_FLUX,
  UNIT_CONV_MAGNETIC_FIELD,
  UNIT_CONV_MAGNETIC_INDUCTANCE,
  UNIT_CONV_TEMPERATURE,
  UNIT_CONV_VOLUME,
  UNIT_CONV_INV_VOLUME
};

void units_init_cgs(struct UnitSystem*);
void units_init(struct UnitSystem*, const struct swift_params*,
                const char* category);
void units_init_default(struct UnitSystem* us,
                        const struct swift_params* params, const char* category,
                        const struct UnitSystem* def);

int units_are_equal(const struct UnitSystem* a, const struct UnitSystem* b);

/* Base units */
double units_get_base_unit(const struct UnitSystem*, enum BaseUnits);
const char* units_get_base_unit_internal_symbol(enum BaseUnits);
const char* units_get_base_unit_cgs_symbol(enum BaseUnits);

/* Cosmology factors */
float units_general_h_factor(const struct UnitSystem* us,
                             const float baseUnitsExponants[5]);
float units_h_factor(const struct UnitSystem* us,
                     enum UnitConversionFactor unit);
float units_general_a_factor(const struct UnitSystem* us,
                             const float baseUnitsExponants[5]);
float units_a_factor(const struct UnitSystem* us,
                     enum UnitConversionFactor unit);

/* Conversion to CGS */
double units_general_cgs_conversion_factor(const struct UnitSystem* us,
                                           const float baseUnitsExponants[5]);
double units_cgs_conversion_factor(const struct UnitSystem* us,
                                   enum UnitConversionFactor unit);
void units_general_cgs_conversion_string(char* buffer,
                                         const struct UnitSystem* us,
                                         const float baseUnitsExponants[5]);
void units_cgs_conversion_string(char* buffer, const struct UnitSystem* us,
                                 enum UnitConversionFactor unit);

/* Conversion between systems */
double units_general_conversion_factor(const struct UnitSystem* from,
                                       const struct UnitSystem* to,
                                       const float baseUnitsExponants[5]);
double units_conversion_factor(const struct UnitSystem* from,
                               const struct UnitSystem* to,
                               enum UnitConversionFactor unit);

#endif /* SWIFT_UNITS_H */
