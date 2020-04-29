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

  /*! Cumulative momentum received from stellar winds (physical internal units)
   */
  float momentum_received;

  /*! Subgrid temperature */
  float subgrid_temp;

  /*! Subgrid density (internal units, physical frame) */
  float subgrid_dens;

#if !defined(COOLING_CHIMES)
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

  union {

    /*! Scale-factor at which the particle was last kicked */
    float last_momentum_kick_scale_factor;

    /*! Time at which the particle was last kicked */
    float last_momentum_kick_time;
  };

  union {

    /*! Scale-factor at which the particle last received energy from type II SNe
     */
    float last_SNII_injection_scale_factor;

    /*! Time at which the particle last received energy from type II SNe */
    float last_SNII_injection_time;
  };

  union {

    /*! Scale-factor at which the particle last received energy from type Ia SNe
     */
    float last_SNIa_injection_scale_factor;

    /*! Time at which the particle last received energy from type Ia SNe */
    float last_SNIa_injection_time;
  };

  union {

    /*! Scale-factor at which the particle last received energy from AGN */
    float last_AGN_injection_scale_factor;

    /*! Time at which the particle last received energy from AGN */
    float last_AGN_injection_time;
  };

  /*! Density of the gas at the last SNII feedback event
   * (physical internal units) */
  float density_at_last_SNII_feedback_event;

  /*! Total amount of AGN feedback energy received by this particle
   * (physical internal units) */
  float AGN_feedback_energy;
};

#endif /* SWIFT_TRACERS_STRUCT_COLIBRE_H */
