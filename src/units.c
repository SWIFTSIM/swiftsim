/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "units.h"

/* Includes. */
#include "adiabatic_index.h"
#include "const.h"
#include "error.h"

/**
 * @brief Initialises the UnitSystem structure with CGS system
 *
 * @param us The UnitSystem to initialize
 */
void units_init_cgs(struct UnitSystem* us) {

  us->UnitMass_in_cgs = 1.;
  us->UnitLength_in_cgs = 1.;
  us->UnitTime_in_cgs = 1.;
  us->UnitCurrent_in_cgs = 1.;
  us->UnitTemperature_in_cgs = 1.;
}

/**
 * @brief Initialises the UnitSystem structure with the constants given in
 * rhe parameter file.
 *
 * @param us The UnitSystem to initialize.
 * @param params The parsed parameter file.
 * @param category The section of the parameter file to read from.
 */
void units_init(struct UnitSystem* us, const struct swift_params* params,
                const char* category) {

  char buffer[200];
  sprintf(buffer, "%s:UnitMass_in_cgs", category);
  us->UnitMass_in_cgs = parser_get_param_double(params, buffer);
  sprintf(buffer, "%s:UnitLength_in_cgs", category);
  us->UnitLength_in_cgs = parser_get_param_double(params, buffer);
  sprintf(buffer, "%s:UnitVelocity_in_cgs", category);
  const double unitVelocity = parser_get_param_double(params, buffer);
  us->UnitTime_in_cgs = us->UnitLength_in_cgs / unitVelocity;
  sprintf(buffer, "%s:UnitCurrent_in_cgs", category);
  us->UnitCurrent_in_cgs = parser_get_param_double(params, buffer);
  sprintf(buffer, "%s:UnitTemp_in_cgs", category);
  us->UnitTemperature_in_cgs = parser_get_param_double(params, buffer);
}

/**
 * @brief Initialises the UnitSystem structure with the constants given in
 * rhe parameter file. Uses a default if the values are not present in the file.
 *
 * @param us The UnitSystem to initialize.
 * @param params The parsed parameter file.
 * @param category The section of the parameter file to read from.
 * @param def The default unit system to copy from if required.
 */
void units_init_default(struct UnitSystem* us,
                        const struct swift_params* params, const char* category,
                        const struct UnitSystem* def) {

  if (!def) error("Default UnitSystem not allocated");

  char buffer[200];
  sprintf(buffer, "%s:UnitMass_in_cgs", category);
  us->UnitMass_in_cgs =
      parser_get_opt_param_double(params, buffer, def->UnitMass_in_cgs);
  sprintf(buffer, "%s:UnitLength_in_cgs", category);
  us->UnitLength_in_cgs =
      parser_get_opt_param_double(params, buffer, def->UnitLength_in_cgs);
  sprintf(buffer, "%s:UnitVelocity_in_cgs", category);
  const double defaultVelocity = def->UnitLength_in_cgs / def->UnitTime_in_cgs;
  const double unitVelocity =
      parser_get_opt_param_double(params, buffer, defaultVelocity);
  us->UnitTime_in_cgs = us->UnitLength_in_cgs / unitVelocity;
  sprintf(buffer, "%s:UnitCurrent_in_cgs", category);
  us->UnitCurrent_in_cgs =
      parser_get_opt_param_double(params, buffer, def->UnitCurrent_in_cgs);
  sprintf(buffer, "%s:UnitTemp_in_cgs", category);
  us->UnitTemperature_in_cgs =
      parser_get_opt_param_double(params, buffer, def->UnitTemperature_in_cgs);
}

/**
 * @brief Returns the base unit conversion factor for a given unit system
 * @param us The UnitSystem used
 * @param baseUnit The base unit
 */
double units_get_base_unit(const struct UnitSystem* us,
                           enum BaseUnits baseUnit) {
  switch (baseUnit) {
    case UNIT_MASS:
      return us->UnitMass_in_cgs;
    case UNIT_LENGTH:
      return us->UnitLength_in_cgs;
    case UNIT_TIME:
      return us->UnitTime_in_cgs;
    case UNIT_CURRENT:
      return us->UnitCurrent_in_cgs;
    case UNIT_TEMPERATURE:
      return us->UnitTemperature_in_cgs;
    default:
      error("Invalid base Unit");
  }
  return 0.0;
}

/**
 * @brief Returns the base unit symbol used internally
 * @param baseUnit The base unit
 */
const char* units_get_base_unit_internal_symbol(enum BaseUnits baseUnit) {
  switch (baseUnit) {
    case UNIT_MASS:
      return "U_M";
    case UNIT_LENGTH:
      return "U_L";
    case UNIT_TIME:
      return "U_t";
    case UNIT_CURRENT:
      return "U_I";
    case UNIT_TEMPERATURE:
      return "U_T";
    default:
      error("Invalid base Unit");
  }
  return "";
}

/**
 * @brief Returns the base unit symbol in the cgs system
 * @param baseUnit The base unit
 */
const char* units_get_base_unit_cgs_symbol(enum BaseUnits baseUnit) {
  switch (baseUnit) {
    case UNIT_MASS:
      return "g";
    case UNIT_LENGTH:
      return "cm";
    case UNIT_TIME:
      return "s";
    case UNIT_CURRENT:
      return "A";
    case UNIT_TEMPERATURE:
      return "K";
    default:
      error("Invalid base Unit");
  }
  return "";
}

void units_get_base_unit_exponants_array(float baseUnitsExp[5],
                                         enum UnitConversionFactor unit) {
  switch (unit) {
    case UNIT_CONV_NO_UNITS:
      break;

    case UNIT_CONV_MASS:
      baseUnitsExp[UNIT_MASS] = 1.f;
      break;

    case UNIT_CONV_LENGTH:
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      break;

    case UNIT_CONV_TIME:
      baseUnitsExp[UNIT_TIME] = 1.f;
      break;

    case UNIT_CONV_FREQUENCY:
      baseUnitsExp[UNIT_TIME] = -1.f;
      break;

    case UNIT_CONV_DENSITY:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = -3.f;
      break;

    case UNIT_CONV_SPEED:
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      baseUnitsExp[UNIT_TIME] = -1.f;
      break;

    case UNIT_CONV_ACCELERATION:
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_FORCE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENERGY:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENERGY_PER_UNIT_MASS:
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENTROPY:
      baseUnitsExp[UNIT_MASS] = 1.f - hydro_gamma;
      baseUnitsExp[UNIT_LENGTH] = 3.f * hydro_gamma - 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ENTROPY_PER_UNIT_MASS:
      baseUnitsExp[UNIT_MASS] = -hydro_gamma;
      baseUnitsExp[UNIT_LENGTH] = 3.f * hydro_gamma - 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_POWER:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -3.f;
      break;

    case UNIT_CONV_PRESSURE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = -1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      break;

    case UNIT_CONV_ELECTRIC_CHARGE:
      baseUnitsExp[UNIT_TIME] = 1.f;
      baseUnitsExp[UNIT_CURRENT] = 1.f;
      break;

    case UNIT_CONV_ELECTRIC_VOLTAGE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -3.f;
      baseUnitsExp[UNIT_CURRENT] = -1.f;
      break;

    case UNIT_CONV_ELECTRIC_CAPACITANCE:
      baseUnitsExp[UNIT_MASS] = -1.f;
      baseUnitsExp[UNIT_LENGTH] = -2.f;
      baseUnitsExp[UNIT_TIME] = 4;
      baseUnitsExp[UNIT_CURRENT] = 2.f;
      break;

    case UNIT_CONV_ELECTRIC_RESISTANCE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -3.f;
      baseUnitsExp[UNIT_CURRENT] = -2.f;
      break;

    case UNIT_CONV_ELECTRIC_CONDUCTANCE:
      baseUnitsExp[UNIT_MASS] = -1.f;
      baseUnitsExp[UNIT_LENGTH] = -2.f;
      baseUnitsExp[UNIT_TIME] = 3.f;
      baseUnitsExp[UNIT_CURRENT] = 2.f;
      break;

    case UNIT_CONV_MAGNETIC_FLUX:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      baseUnitsExp[UNIT_CURRENT] = -1.f;
      break;

    case UNIT_CONV_MAGNETIC_FIELD:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      baseUnitsExp[UNIT_CURRENT] = -1.f;
      break;

    case UNIT_CONV_MAGNETIC_INDUCTANCE:
      baseUnitsExp[UNIT_MASS] = 1.f;
      baseUnitsExp[UNIT_LENGTH] = 2.f;
      baseUnitsExp[UNIT_TIME] = -2.f;
      baseUnitsExp[UNIT_CURRENT] = -2.f;
      break;

    case UNIT_CONV_TEMPERATURE:
      baseUnitsExp[UNIT_TEMPERATURE] = 1.f;
      break;

    case UNIT_CONV_VOLUME:
      baseUnitsExp[UNIT_LENGTH] = 3.f;
      break;

    case UNIT_CONV_INV_VOLUME:
      baseUnitsExp[UNIT_LENGTH] = -3.f;
      break;

    default:
      error("Invalid choice of pre-defined units");
      break;
  }
}

/**
 * @brief Returns the conversion factor for a given unit in the chosen unit
 * system
 * @param us The system of units in use
 * @param unit The unit to convert
 */
double units_cgs_conversion_factor(const struct UnitSystem* us,
                                   enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  units_get_base_unit_exponants_array(baseUnitsExp, unit);

  return units_general_cgs_conversion_factor(us, baseUnitsExp);
}

/**
 * @brief Returns the h factor exponentiation for a given unit
 * @param us The system of units in use
 * @param unit The unit to convert
 */
float units_h_factor(const struct UnitSystem* us,
                     enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  units_get_base_unit_exponants_array(baseUnitsExp, unit);

  return units_general_h_factor(us, baseUnitsExp);
}

/**
 * @brief Returns the scaling factor exponentiation for a given unit
 * @param us The system of units in use
 * @param unit The unit to convert
 */
float units_a_factor(const struct UnitSystem* us,
                     enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  units_get_base_unit_exponants_array(baseUnitsExp, unit);

  return units_general_a_factor(us, baseUnitsExp);
}

/**
 * @brief Returns a string containing the exponents of the base units making up
 * the conversion factors
 */
void units_cgs_conversion_string(char* buffer, const struct UnitSystem* us,
                                 enum UnitConversionFactor unit) {
  float baseUnitsExp[5] = {0.f};

  units_get_base_unit_exponants_array(baseUnitsExp, unit);

  units_general_cgs_conversion_string(buffer, us, baseUnitsExp);
}

/**
 * @brief Returns the conversion factor for a given unit (expressed in terms of
 * the 5 fundamental units) in the chosen unit system
 * @param us The unit system used
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
double units_general_cgs_conversion_factor(const struct UnitSystem* us,
                                           const float baseUnitsExponants[5]) {
  double factor = 1.;

  for (int i = 0; i < 5; ++i)
    if (baseUnitsExponants[i] != 0)
      factor *= pow(units_get_base_unit(us, (enum BaseUnits)i),
                    baseUnitsExponants[i]);
  return factor;
}

/**
 * @brief Returns the h factor exponentiation for a given unit (expressed in
 * terms of the 5 fundamental units)
 * @param us The unit system used
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
float units_general_h_factor(const struct UnitSystem* us,
                             const float baseUnitsExponants[5]) {
  float factor_exp = 0.f;

  factor_exp += -baseUnitsExponants[UNIT_MASS];
  factor_exp += -baseUnitsExponants[UNIT_LENGTH];
  factor_exp += -baseUnitsExponants[UNIT_TIME];

  return factor_exp;
}

/**
 * @brief Returns the scaling factor exponentiation for a given unit (expressed
 * in terms of the 5 fundamental units)
 * @param us The unit system used
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
float units_general_a_factor(const struct UnitSystem* us,
                             const float baseUnitsExponants[5]) {
  float factor_exp = 0.f;

  factor_exp += baseUnitsExponants[UNIT_LENGTH];

  return factor_exp;
}

/**
 * @brief Returns a string containing the exponents of the base units making up
 * the conversion factors (expressed in terms of the 5 fundamental units)
 * @param buffer The buffer in which to write (The buffer must be long enough,
 * 140 chars at most)
 * @param us The UnitsSystem in use.
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
void units_general_cgs_conversion_string(char* buffer,
                                         const struct UnitSystem* us,
                                         const float baseUnitsExponants[5]) {
  char temp[14];
  const double a_exp = units_general_a_factor(us, baseUnitsExponants);
  const double h_exp = units_general_h_factor(us, baseUnitsExponants);

  /* Check whether we are unitless or not */
  char isAllNonZero = 1;
  for (int i = 0; i < 5; ++i)
    if (baseUnitsExponants[i] != 0.) isAllNonZero = 0;

  if (isAllNonZero) {
    sprintf(buffer, "[ - ] ");
    return;
  }

  /* Add a-factor */
  if (a_exp == 0)
    sprintf(buffer, " ");
  else if (a_exp == 1)
    sprintf(buffer, "a ");
  else if (remainder(a_exp, 1.) == 0)
    sprintf(buffer, "a^%d ", (int)a_exp);
  else
    sprintf(buffer, "a^%7.4f ", a_exp);

  /* Add h-factor */
  if (h_exp == 0)
    sprintf(temp, " ");
  else if (h_exp == 1)
    sprintf(temp, "h ");
  else if (remainder(h_exp, 1.) == 0)
    sprintf(temp, "h^%d ", (int)h_exp);
  else
    sprintf(temp, "h^%7.4f ", h_exp);
  strncat(buffer, temp, 12);

  /* Add conversion units */
  for (int i = 0; i < 5; ++i)
    if (baseUnitsExponants[i] != 0) {
      if (baseUnitsExponants[i] == 0.)
        sprintf(temp, " ");
      else if (baseUnitsExponants[i] == 1.)
        sprintf(temp, "%s ",
                units_get_base_unit_internal_symbol((enum BaseUnits)i));
      else if (remainder(baseUnitsExponants[i], 1.) == 0)
        sprintf(temp, "%s^%d ",
                units_get_base_unit_internal_symbol((enum BaseUnits)i),
                (int)baseUnitsExponants[i]);
      else
        sprintf(temp, "%s^%7.4f ",
                units_get_base_unit_internal_symbol((enum BaseUnits)i),
                baseUnitsExponants[i]);
      strncat(buffer, temp, 12);
    }

  /* Add CGS units */
  strncat(buffer, " [ ", 3);

  for (int i = 0; i < 5; ++i) {
    if (baseUnitsExponants[i] != 0) {
      if (baseUnitsExponants[i] == 0.)
        continue;
      else if (baseUnitsExponants[i] == 1.)
        sprintf(temp, "%s ", units_get_base_unit_cgs_symbol((enum BaseUnits)i));
      else if (remainder(baseUnitsExponants[i], 1.) == 0)
        sprintf(temp, "%s^%d ",
                units_get_base_unit_cgs_symbol((enum BaseUnits)i),
                (int)baseUnitsExponants[i]);
      else
        sprintf(temp, "%s^%7.4f ",
                units_get_base_unit_cgs_symbol((enum BaseUnits)i),
                baseUnitsExponants[i]);
      strncat(buffer, temp, 12);
    }
  }

  strncat(buffer, "]", 2);
}

/**
 * @brief Are the two unit systems equal ?
 *
 * @param a The First #UnitSystem
 * @param b The second #UnitSystem
 * @return 1 if the systems are the same, 0 otherwise
 */
int units_are_equal(const struct UnitSystem* a, const struct UnitSystem* b) {

  if (a->UnitMass_in_cgs != b->UnitMass_in_cgs) return 0;
  if (a->UnitLength_in_cgs != b->UnitLength_in_cgs) return 0;
  if (a->UnitTime_in_cgs != b->UnitTime_in_cgs) return 0;
  if (a->UnitCurrent_in_cgs != b->UnitCurrent_in_cgs) return 0;
  if (a->UnitTemperature_in_cgs != b->UnitTemperature_in_cgs) return 0;

  return 1;
}

/**
 * @brief Return the unit conversion factor between two systems
 *
 * @param from The #UnitSystem we are converting from
 * @param to The #UnitSystem we are converting to
 * @param baseUnitsExponants The exponent of each base units required to form
 * the desired quantity. See conversionFactor() for a working example
 */
double units_general_conversion_factor(const struct UnitSystem* from,
                                       const struct UnitSystem* to,
                                       const float baseUnitsExponants[5]) {

  const double from_cgs =
      units_general_cgs_conversion_factor(from, baseUnitsExponants);
  const double to_cgs =
      units_general_cgs_conversion_factor(to, baseUnitsExponants);

  return from_cgs / to_cgs;
}

/**
 * @brief Return the unit conversion factor between two systems
 *
 * @param from The #UnitSystem we are converting from
 * @param to The #UnitSystem we are converting to
 * @param unit The unit we are converting
 *
 * @return The conversion factor
 */
double units_conversion_factor(const struct UnitSystem* from,
                               const struct UnitSystem* to,
                               enum UnitConversionFactor unit) {

  float baseUnitsExp[5] = {0.f};

  units_get_base_unit_exponants_array(baseUnitsExp, unit);

  return units_general_conversion_factor(from, to, baseUnitsExp);
}
