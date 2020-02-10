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
#ifndef SWIFT_TRACERS_STRUCT_COLIBRE_H
#define SWIFT_TRACERS_STRUCT_COLIBRE_H

/**
 * @brief Properties of the tracers stored in the extended particle data.
 */
struct tracers_xpart_data {

  /*! Maximum temperature achieved by this particle */
  float maximum_temperature;

  float momentum_received;

  /* subgrid temperature */
  float subgrid_temp;

  /* subgrid density */
  float subgrid_dens;

#if !defined(COOLING_CHIMES) && !defined(COOLING_CHIMES_HYBRID)
  /*! Hydrogen fractions */
  float nHI_over_nH;
  float nHII_over_nH;
  float nH2_over_nH;
#endif

  /*! HII region timer */
  float HIIregion_timer_gas;

  /*! Id of the star particle responsible for HII region */
  long long HIIregion_starid;

  /*! Anonymous union for the cosmological non-cosmological runs distinction */
  union {

    /*! Scale-factor at which the maximal temperature was reached */
    float maximum_temperature_scale_factor;

    /*! Time at which the maximal temperature was reached */
    float maximum_temperature_time;
  };
};

#endif /* SWIFT_TRACERS_STRUCT_COLIBRE_H */
