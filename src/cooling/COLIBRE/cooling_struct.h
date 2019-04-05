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
#ifndef SWIFT_COOLING_STRUCT_COLIBRE_H
#define SWIFT_COOLING_STRUCT_COLIBRE_H

#define eagle_table_path_name_length 500

/**
 * @brief struct containing cooling tables
 */
struct cooling_tables {

  /* array of all cooling processes (temperature) */
  float *Tcooling;

  /* array of all cooling processes (internal energy) */
  float *Ucooling;

  /* array of all heating processes (temperature) */
  float *Theating;

  /* array of all heating processes (internal energy) */
  float *Uheating;

  /* array of all electron abundances (temperature) */
  float *Telectron_fraction;

  /* array of all electron abundances (internal energy) */
  float *Uelectron_fraction;

  /* array to get T from U */
  float *T_from_U;

  /* array to get U from T */
  float *U_from_T;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling tables */
  struct cooling_tables table;

  /*! Redshift bins */
  float *Redshifts;

  /*! Hydrogen number density bins */
  float *nH;

  /*! Temperature bins */
  float *Temp;

  /*! Metallicity bins */
  float *Metallicity;

  /*! Internal energy bins */
  float *Therm;

  /*! Abundance ratios for each metallicity bin and for each included element */
  float *LogAbundances;
  float *Abundances;
  float *Abundances_inv;

  /*! Atomic masses for all included elements */
  float *atomicmass;

  /*! Mass fractions of all included elements */
  float *LogMassFractions;
  float *MassFractions;

  /*! Index for solar metallicity in the metallicity dimension */
  int indxZsol;

  /*! Solar metallicity */
  float *Zsol;

  /*! Filepath to the directory containing the HDF5 cooling tables */
  char cooling_table_path[eagle_table_path_name_length];

  /*! Redshit of H reionization */
  float H_reion_z;
  /*! Ca over Si abundance divided by the solar ratio for these elements */
  float Ca_over_Si_ratio_in_solar;

  /*! S over Si abundance divided by the solar ratio for these elements */
  float S_over_Si_ratio_in_solar;

  /*! Redshift of He reionization */
  float He_reion_z_centre;

  /*! Spread of the He reionization */
  float He_reion_z_sigma;

  /*! He reionization energy in CGS units */
  float He_reion_heat_cgs;

  /*! Internal energy conversion from internal units to CGS (for quick access)
   */
  double internal_energy_to_cgs;

  /*! Internal energy conversion from CGS to internal units (for quick access)
   */
  double internal_energy_from_cgs;

  /*! Number density conversion from internal units to CGS (for quick access) */
  double number_density_to_cgs;

  /*! Inverse of proton mass in cgs (for quick access) */
  double inv_proton_mass_cgs;

  /*! Temperatur of the CMB at present day (for quick access) */
  double T_CMB_0;

  /*! Compton rate in cgs units */
  double compton_rate_cgs;

  /*! Are we doing Newton-Raphson iterations? */
  int newton_flag;
};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {

  /*! Cumulative energy radiated by the particle */
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_COLIBRE_H */
