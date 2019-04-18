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

#define element_H 0
#define element_He 1
#define element_C 2
#define element_N 3
#define element_O 4
#define element_Ne 5
#define element_Mg 6
#define element_Si 7
#define element_S 8
#define element_Ca 9
#define element_Fe 10
#define element_OA 11

#define cooltype_H2 12
#define cooltype_molecules 13
#define cooltype_HD 14
#define cooltype_NetFFH 15
#define cooltype_NetFFM 16
#define cooltype_eeBrems 17
#define cooltype_Compton 18
#define cooltype_Dust 19

#define heattype_H2 12
#define heattype_COdiss 13
#define heattype_CosmicRay 14
#define heattype_UTA 15
#define heattype_line 16
#define heattype_Hlin 17
#define heattype_ChaT 18
#define heattype_HFF 19
#define heattype_Compton 20
#define heattype_Dust 21

void get_cooling_redshifts(struct cooling_function_data *cooling);

void read_cooling_header(struct cooling_function_data *cooling);
void read_cooling_tables(struct cooling_function_data *cooling);

#endif
