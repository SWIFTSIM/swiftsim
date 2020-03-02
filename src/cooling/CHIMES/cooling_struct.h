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
#include "colibre_tables.h"

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
  int use_redshift_dependent_UVB;

  /* User parameter to scale the
   * normalisation of the radiation
   * field by a constant factor. */
  ChimesFloat rad_field_norm_factor;

  /* Factor to re-scale shielding
   * length. */
  double shielding_length_factor;

  /* Maximum shielding length, in code,
   * units. If negative, do not impose
   * a maximum. */
  double max_shielding_length;

  /* Parameters used for the
   * COLIBRE ISRF. */
  double N_H0;

  /* Flags to control eqm mode and
   * thermal evolution. */
  int ChemistryEqmMode;
  int ThermEvolOn;

  /* Flag to determine how we set
   * the initial CHIMES abundances. */
  int init_abundance_mode;

  /* For init_abundance_mode == 1, all
   * elements are initially set to one
   * ionisation state, determined by
   * the InitIonState parameter. */
  int InitIonState;

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

  /* Flag to implement dust depletion
   * as in COLIBRE. */
  int colibre_metal_depletion;

  /* Dust depletion factors in the solar neighbourhood.
   * Taken from Ploeckinger et al. (in prep). */
  double f_dust0_C;
  double f_dust0_O;
  double f_dust0_Mg;
  double f_dust0_Si;
  double f_dust0_Ca;
  double f_dust0_Fe;

  /* Distance from EOS to use thermal equilibrium
   * temperature for subgrid props, and to evolve
   * the chemistry in equilibrium. */
  float dlogT_EOS;

  /* Flag to use the subgrid density and
   * temperature from the COLIBRE cooling tables. */
  int use_colibre_subgrid_EOS;

  /* Flag to set particles to chemical equilibrium 
   * if they have been heated by feedback. */ 
  int set_FB_particles_to_eqm; 

  /* Threshold to switch between rapid and slow cooling regimes. */
  double rapid_cooling_threshold;

  /* Ionization fraction of gas particles tagged as HII regions */
  float HIIregion_fion;

  /* Temperature of gas particles tagged as HII regions */
  float HIIregion_temp;

  /* Colibre cooling table */
  struct colibre_cooling_tables colibre_table;
};

/**
 * @brief Properties of the cooling stored in the particle data
 */
struct cooling_xpart_data {
  /* CHIMES abundance array */
  double chimes_abundances[CHIMES_NETWORK_SIZE];

  /* Cumulative energy radiated by the particle */
  float radiated_energy;

  /* Flag to indicate when the particle has been heated by feedback. */ 
  int heated_by_FB; 
};

#endif /* SWIFT_COOLING_STRUCT_CHIMES_H */
