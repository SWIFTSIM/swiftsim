/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COLIBRE_COOL_TABLES_H
#define SWIFT_COLIBRE_COOL_TABLES_H

/**
 * @file src/cooling/COLIBRE/cooling.h
 * @brief COLIBRE cooling function
 */

/* Config parameters. */
#include "../config.h"

#include "cooling_struct.h"

/*! Number of different bins along the temperature axis of the tables */
#define colibre_cooling_N_temperature 86

/*! Number of different bins along the redshift axis of the tables */
#define colibre_cooling_N_redshifts 10

/*! Number of different bins along the density axis of the tables */
#define colibre_cooling_N_density 71

/*! Number of different bins along the metallicity axis of the tables */
#define colibre_cooling_N_metallicity 11

/*! Number of different bins along the internal energy axis of the tables */
#define colibre_cooling_N_internalenergy 191

/*! Number of different cooling channels in the tables */
#define colibre_cooling_N_cooltypes 22

/*! Number of different heating channels in the tables */
#define colibre_cooling_N_heattypes 24

/*! Number of different electron fractions (each element - other atoms
 *  + tot prim + tot metal + tot)  in the tables */
#define colibre_cooling_N_electrontypes 14

/*! Number of different elements in the tables */
#define colibre_cooling_N_elementtypes 12

enum {
  element_H,
  element_He,
  element_C,
  element_N,
  element_O,
  element_Ne,
  element_Mg,
  element_Si,
  element_S,
  element_Ca,
  element_Fe,
  element_OA
};

enum {
  cooltype_H2 = 12,
  cooltype_molecules,
  cooltype_HD,
  cooltype_NetFFH,
  cooltype_NetFFM,
  cooltype_eeBrems,
  cooltype_Compton,
  cooltype_Dust
};

enum {
  heattype_H2 = 12,
  heattype_COdiss,
  heattype_CosmicRay,
  heattype_UTA,
  heattype_line,
  heattype_Hlin,
  heattype_ChaT,
  heattype_HFF,
  heattype_Compton,
  heattype_Dust
};

void get_cooling_redshifts(struct cooling_function_data *cooling);
void read_cooling_header(struct cooling_function_data *cooling);
void read_cooling_tables(struct cooling_function_data *cooling);

#endif
