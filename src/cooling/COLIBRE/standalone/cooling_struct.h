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
#define table_path_name_length 500

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
 * @brief Properties of the particle / resolution element
 */
struct particle_data {
   float *metal_mass_fraction;
   float *metal_mass;
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

  /*! Constant factor for Compton cooling */
  double Compton_constant_cgs;

  /*! Filepath to the directory containing the HDF5 cooling tables */
  char cooling_table_path[table_path_name_length];

};

