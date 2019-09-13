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

  /* Temperature of the CMB at present day */
  double T_CMB_0;

};

/**
 * @brief Properties of the cooling stored in the particle data
 */
struct cooling_xpart_data {};

#endif /* SWIFT_COOLING_STRUCT_CHIMES_H */
