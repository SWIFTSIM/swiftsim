/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_SNIA_DTD_STRUCT_H
#define SWIFT_SNIA_DTD_STRUCT_H

/**
 * @file src/snia_dtd.h
 * @brief Branches between the different SNIa delay time distributions recipies.
 */

/* Config parameters. */
#include "../config.h"

/* Import the right SNIa delay time distribution (DTD) struct definition */
#if defined(SNIA_DTD_EXP)

/**
 * @brief SNIa struct for DTD properties
 */
struct dtd {

  /*! SNIa time scale */
  double SNIa_timescale_Gyr;

  /*! Inverse SNIa time scale*/
  double SNIa_timescale_Gyr_inv;

  /*! Efficiency of the SNIa model */
  double SNIa_efficiency;
};

#elif defined(SNIA_DTD_POWER)

/**
 * @brief SNIa struct for DTD properties
 */
struct dtd {

  /*! Efficiency of the SNIa model */
  double SNIa_efficiency;

  /*! Power of power law */
  double power;

  /*! Normalization time */
  double normalization_timescale_Gyr;

  /*! Delay time */
  double delay_time_Gyr;

  /*! DTD normalization */
  double norm;
};

#elif defined(SNIA_DTD_POWER_BETA1)

/**
 * @brief SNIa struct for DTD properties
 */
struct dtd {

  /*! Efficiency of the SNIa model */
  double SNIa_efficiency;

  /*! Normalization time */
  double normalization_timescale_Gyr;

  /*! Delay time */
  double delay_time_Gyr;

  /*! DTD normalization */
  double norm;
};

#elif defined(SNIA_DTD_GAUSSIAN)

/**
 * @brief SNIa struct for DTD properties
 */
struct dtd {

  /*! Efficiency of the SNIa model constant part */
  double SNIa_efficiency_const;

  /*! Efficiency of the SNIa model constant part */
  double SNIa_efficiency_gauss;

  /*! Normalization time */
  double normalization_timescale_Gyr;

  /*! Characteristic time */
  double characteristic_time_Gyr;

  /*! deviation of the characteristic time */
  double std_characteristic_time_Gyr;

  /*! inverse of deviation of the characteristic time */
  double std_characteristic_time_Gyr_inv;

  /*! Delay time */
  double delay_time_Gyr;

  /*! DTD normalization */
  double norm_const;
};

#elif defined(SNIA_DTD_CONST)

/**
 * @brief SNIa struct for DTD properties
 */
struct dtd {

  /*! SNIa time scale */
  double normalization_timescale_Gyr;

  /*! normalization of the SNIa constant DTD */
  double norm;

  /*! Efficiency of the SNIa model */
  double SNIa_efficiency;
};

#elif defined(SNIA_DTD_BROKEN_POWER_LAW)

/**
 * @brief SNIa struct for DTD properties
 */
struct dtd {

  /*! Efficiency of the SNIa model */
  double SNIa_efficiency;

  /*! Power short time power law */
  double power_short_time;

  /*! Power long time power law */
  double power_long_time;

  /*! Power law break time */
  double break_time_Gyr;

  /*! Normalization time scale */
  double normalization_timescale_Gyr;

  /*! Delay time */
  double delay_time_Gyr;

  /*! DTD normalization short time */
  double norm_short;

  /*! DTD normalization long time */
  double norm_long;
};

#else
#error "Invalid choice of the SNIa delay time distribution (DTD)"
#endif

#endif /* SWIFT_SNIA_DTD_STRUCT_H */
 
