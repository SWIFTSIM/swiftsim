/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H


/**
 * @brief Model for the initial mass function.
 *
 * Describe a model such as Kroupa 2001:
 *
 * f(m) = coef[i] * pow(m, exp[i])
 */
struct initial_mass_function {

  /*! Mass limits between IMF parts (n_parts + 1 elements) */
  float *mass_limits;

  /*! Exponent of each IMF parts (n_parts elements) */
  float *exp;

  /*! Coefficient of each IMF parts (n_parts elements) */
  float *coef;

  /*! Number of parts in the function */
  int n_parts;

  /*! Minimal mass contained in mass_limits, copied for more clarity*/
  float mass_min;

  /*! Maximal mass contained in mass_limits, copied for more clarity */
  float mass_max;

};

/**
 * @brief Model for the stellar lifetime.
 */
struct lifetime {

  /*! Coefficients for the log10(m)^2 term */
  float quadratic[3];

  /*! Coefficients for the log10(m) term */
  float linear[3];

  /*! Coefficients for the constant term */
  float constant[3];
};

/**
 * @brief Model for SNIa.
 */
struct supernovae_ia {
  /*! Yields TODO more comment */
  struct interpolation *yields;

};

/**
 * @brief Model for SNII.
 */
struct supernovae_ii {

  /*! Yields TODO more comment */
  struct interpolation *yields;
  
};

/**
 * @brief The complete stellar model.
 */
struct stellar_model {
  // TODO elements
  
  /*! The initial mass function */
  struct initial_mass_function imf;

  /*! The stellar lifetime */
  struct lifetime lifetime;

  /*! The supernovae type Ia */
  struct supernovae_ia snia;

  /*! The supernovae type II */
  struct supernovae_ii snii;

};



#endif // SWIFT_STELLAR_EVOLUTION_STRUCT_GEAR_H
