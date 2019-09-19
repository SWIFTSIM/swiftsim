/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_STRUCT_CHIMES_H
#define SWIFT_COOLING_STRUCT_CHIMES_H

/**
 * @file src/cooling/CHIMES/cooling_struct.h
 * @brief Infrastructure for CHIMES cooling. 
 */

/* Local includes. */
#include "cooling/CHIMES/chimes/chimes_vars.h" 
#include "cooling/CHIMES/chimes/chimes_proto.h"

/* Maximum size of CHIMES abundance arrays. */ 
#define CHIMES_NETWORK_SIZE 157 

/**
 * @brief Properties of the cooling function. 
 * This includes the globalVaraibles structure 
 * that holds the parameters used to control 
 * the behaviour of the CHIMES module. 
 */
struct cooling_function_data {

  /* CHIMES global variables. */ 
  struct globalVariables ChimesGlobalVars; 

  /* Flags to control UV field and 
   * shielding options. */ 
  int UV_field_flag; 
  int Shielding_flag; 

  /* Stores normalisation of user-supplied 
   * radiation fields. */ 
  ChimesFloat *isotropic_photon_density; 

  /* Flags to control eqm mode and 
   * thermal evolution. */ 
  int ChemistryEqmMode; 
  int ThermEvolOn; 
  
  /* Cosmic ray ionisation rate of HI. */ 
  double cosmic_ray_rate; 

  /* Temperature of the CMB at present day */
  double T_CMB_0;

  /* Parameters to compute S and Ca 
   * from Si. */ 
  float S_over_Si_ratio_in_solar; 
  float Ca_over_Si_ratio_in_solar; 
  
  float Si_solar_mass_fraction; 
  float S_solar_mass_fraction; 
  float Ca_solar_mass_fraction; 

  /* Solar metallicity */ 
  float Zsol; 
};

/**
 * @brief Properties of the cooling stored in the particle data
 */
struct cooling_xpart_data {
  /* CHIMES abundance array */ 
  ChimesFloat chimes_abundances[CHIMES_NETWORK_SIZE]; 

  /* Cumulative energy radiated by the particle */
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_CHIMES_H */
