/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_DEBUG_H
#define SWIFT_RT_IO_DEBUG_H

#include "io_properties.h"

/**
 * @file src/rt/debug/rt_io.h
 * @brief Main header file for the debug radiative transfer scheme IO routines.
 */

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  list[0] =
      io_make_output_field("RTStarIact", INT, 1, UNIT_CONV_NO_UNITS, 0, parts,
                           rt_data.iact_stars_inject,
                           "number of interactions between this hydro particle"
                           " and any star particle during injection step");
  list[1] = io_make_output_field("RTTotalCalls", INT, 1, UNIT_CONV_NO_UNITS, 0,
                                 parts, rt_data.calls_tot,
                                 "total number of calls to this "
                                 "particle during the run");
  list[2] = io_make_output_field("RTCallsThisStep", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, parts, rt_data.calls_per_step,
                                 "number of calls "
                                 "to this particle during one time step");
  list[3] =
      io_make_output_field("RTCallsSelfInjection", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_self_inject,
                           "number of calls to this particle during one time "
                           "step in injection self task");
  list[4] =
      io_make_output_field("RTCallsPairInjection", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_pair_inject,
                           "number of calls to this particle during one time "
                           "step in injection pair task");
  list[5] =
      io_make_output_field("RTPhotonsUpdated", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.photon_number_updated,
                           "=1 if photon number has been updated in this step");
  list[6] =
      io_make_output_field("RTCallsIactGradient", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_iact_gradient,
                           "number of calls to this particle during the"
                           "gradient interaction loop");
  list[7] =
      io_make_output_field("RTCallsIactGradientSym", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_iact_gradient_sym,
                           "number of symmetric calls to this particle during the"
                           "gradient interaction loop");
  list[8] =
      io_make_output_field("RTCallsIactGradientNonSym", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_iact_gradient_nonsym,
                           "number of nonsymmetric calls to this particle during the"
                           "gradient interaction loop");
  list[9] =
      io_make_output_field("RTCallsIactTransport", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_iact_transport,
                           "number of calls to this particle during the"
                           "transport interaction loop");
  list[10] =
      io_make_output_field("RTCallsIactTransportSym", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_iact_transport_sym,
                           "number of symmetric calls to this particle during the"
                           "transport interaction loop");
  list[11] =
      io_make_output_field("RTCallsIactTransportNonSym", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_iact_transport_nonsym,
                           "number of non-symmetric calls to this particle during the"
                           "transport interaction loop");

  list[12] =
      io_make_output_field("RTNeighsIactGrad", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.neigh_iact_grad,
                           "neighbours interacted with during gradient loop");
  list[13] =
      io_make_output_field("RTNeighCellsIactGrad", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.neigh_cell_iact_grad,
                           "cells of neighbours interacted with during gradient loop");
  list[14] =
      io_make_output_field("RTNrNeighIactGrad", INT, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.neigh_iact_grad_free,
                           "number of neighbours interacted with during gradient loop");
  list[15] =
      io_make_output_field("RTNeighsIactTransport", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.neigh_iact_transp,
                           "neighbours interacted with during gradient loop");

  list[16] =
      io_make_output_field("RTNeighCellsIactTransport", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.neigh_cell_iact_transp,
                           "cells of neighbours interacted with during transport loop");

  list[17] =
      io_make_output_field("RTNrNeighIactTransport", INT, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.neigh_iact_transp_free,
                           "number of neighbours interacted with during transport loop");
  list[18] =
      io_make_output_field("RTThisCellGrad", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.this_cell_grad,
                           "the cell this particle is in, from gradient tasks");
  list[19] =
      io_make_output_field("RTThisCellTransport", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.this_cell_transport,
                           "the cell this particle is in, from transport tasks");


  list[20] =
      io_make_output_field("RTHydroNeighsIactGrad", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_neigh_iact_grad,
                           "neighbours interacted with during gradient loop");
  list[21] =
      io_make_output_field("RTHydroNeighCellsIactGrad", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_neigh_cell_iact_grad,
                           "cells of neighbours interacted with during gradient loop");
  list[22] =
      io_make_output_field("RTHydroNrNeighIactGrad", INT, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_neigh_iact_grad_free,
                           "number of neighbours interacted with during gradient loop");
  list[23] =
      io_make_output_field("RTHydroNeighsIactTransport", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_neigh_iact_transp,
                           "neighbours interacted with during gradient loop");

  list[24] =
      io_make_output_field("RTHydroNeighCellsIactTransport", LONGLONG, 400, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_neigh_cell_iact_transp,
                           "cells of neighbours interacted with during transport loop");

  list[25] =
      io_make_output_field("RTHydroNrNeighIactTransport", INT, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_neigh_iact_transp_free,
                           "number of neighbours interacted with during transport loop");
  list[26] =
      io_make_output_field("RTHydroThisCellGrad", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_this_cell_grad,
                           "the cell this particle is in, according to gradient tasks");
  list[27] =
      io_make_output_field("RTHydroThisCellTransport", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0, parts, rt_data.hydro_this_cell_transport, 
                           "the cell this particle is in, according to force tasks");

  list[28] =
      io_make_output_field("RTHydroCallsIactGradient", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_hydro_iact_gradient,
                           "number of calls to this particle during the"
                           "gradient interaction loop");
  list[29] =
      io_make_output_field("RTHydroCallsIactGradientSym", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_hydro_iact_gradient_sym,
                           "number of symmetric calls to this particle during the"
                           "gradient interaction loop");
  list[30] =
      io_make_output_field("RTHydroCallsIactGradientNonSym", INT, 1, UNIT_CONV_NO_UNITS, 0,
                           parts, rt_data.calls_hydro_iact_gradient_nonsym,
                           "number of nonsymmetric calls to this particle during the"
                           "gradient interaction loop");
  list[31] =
      io_make_output_field("RTHydroCallsIactForce", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_hydro_iact_force,
                           "number of calls to this particle during the"
                           "force interaction loop");
  list[32] =
      io_make_output_field("RTHydroCallsIactForceSym", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_hydro_iact_force_sym,
                           "number of symmetric calls to this particle during the"
                           "force interaction loop");
  list[33] =
      io_make_output_field("RTHydroCallsIactForceNonSym", INT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.calls_hydro_iact_force_nonsym,
                           "number of non-symmetric calls to this particle during the"
                           "force interaction loop");
  list[34] =
      io_make_output_field("RTSmlGrad", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.h_grad,
                           "h used in RT gradient interaction loop");
  list[35] =
      io_make_output_field("RTSmlTransport", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.h_transport,
                           "h used in RT transport interaction loop");
  list[36] =
      io_make_output_field("RTHydroSmlGrad", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.h_hydro_grad,
                           "h used in hydro gradient interaction loop");
  list[37] =
      io_make_output_field("RTHydroSmlForce", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           0, parts, rt_data.h_force,
                           "h used in hydro force interaction loop");
  return 38;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {

  list[0] = io_make_output_field("RTHydroIact", INT, 1, UNIT_CONV_NO_UNITS, 0,
                                 sparts, rt_data.iact_hydro_inject,
                                 "number of interactions between this hydro "
                                 "particle and any star particle");
  list[1] = io_make_output_field("RTTotalCalls", INT, 1, UNIT_CONV_NO_UNITS, 0,
                                 sparts, rt_data.calls_tot,
                                 "total number of calls "
                                 "to this particle during the run");
  list[2] = io_make_output_field("RTCallsThisStep", INT, 1, UNIT_CONV_NO_UNITS,
                                 0, sparts, rt_data.calls_per_step,
                                 "number of calls to "
                                 "this particle during one time step");
  list[3] =
      io_make_output_field("RTCallsSelfInjection", INT, 1, UNIT_CONV_NO_UNITS,
                           0, sparts, rt_data.calls_self_inject,
                           "number of calls to this particle during one time "
                           "step in injection self task");
  list[4] =
      io_make_output_field("RTCallsPairInjection", INT, 1, UNIT_CONV_NO_UNITS,
                           0, sparts, rt_data.calls_pair_inject,
                           "number of calls to this particle during one time "
                           "step in injection pair task");
  list[5] = io_make_output_field(
      "RTEmissionRateSet", INT, 1, UNIT_CONV_NO_UNITS, 0, sparts,
      rt_data.emission_rate_set, "Stellar photon emission rates set?");

  return 6;
}
#endif /* SWIFT_RT_IO_DEBUG_H */
