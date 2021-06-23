/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_GEAR_RT_RIEMANN_GLF_H
#define SWIFT_GEAR_RT_RIEMANN_GLF_H

#include "rt_getters.h"
#include "rt_parameters.h"


/**
 * @file src/rt/GEAR/rt_riemann_GLF.h
 * @brief  The Global Lax-Friedrich (GLF) Riemann solver for the moments of the
 * radiative transfer equation following Rosdahl et al 2013
 * */

/**
 * @brief Solve the Riemann problem for the RT equations and return the
 * flux at the interface.
 *
 * @param UL left state
 * @param UR right state
 * @param flux_half the resulting flux at the interface
 */
__attribute__((always_inline)) INLINE static void rt_riemann_solve_for_flux(
    const float UL[4], const float UR[4], float flux_half[4][3]) {


  float fluxL[4][3];
  rt_get_hyperbolic_flux(UL, fluxL);
  float fluxR[4][3];
  rt_get_hyperbolic_flux(UR, fluxR);

  const float c_red = rt_params.reduced_speed_of_light;
  const float cdU[4] = {c_red * (UR[0] - UL[0]), c_red * (UR[1] - UL[1]), c_red * (UR[2] - UL[2]), c_red * (UR[3] - UL[3])};

  flux_half[0][0] = 0.5f * (fluxL[0][0] + fluxR[0][0] - cdU[0]);
  flux_half[0][1] = 0.5f * (fluxL[0][1] + fluxR[0][1] - cdU[0]);
  flux_half[0][2] = 0.5f * (fluxL[0][2] + fluxR[0][2] - cdU[0]);

  flux_half[1][0] = 0.5f * (fluxL[1][0] + fluxR[1][0] - cdU[1]);
  flux_half[1][1] = 0.5f * (fluxL[1][1] + fluxR[1][1] - cdU[1]);
  flux_half[1][2] = 0.5f * (fluxL[1][2] + fluxR[1][2] - cdU[1]);

  flux_half[2][0] = 0.5f * (fluxL[2][0] + fluxR[2][0] - cdU[2]);
  flux_half[2][1] = 0.5f * (fluxL[2][1] + fluxR[2][1] - cdU[2]);
  flux_half[2][2] = 0.5f * (fluxL[2][2] + fluxR[2][2] - cdU[2]);

  flux_half[3][0] = 0.5f * (fluxL[3][0] + fluxR[3][0] - cdU[3]);
  flux_half[3][1] = 0.5f * (fluxL[3][1] + fluxR[3][1] - cdU[3]);
  flux_half[3][2] = 0.5f * (fluxL[3][2] + fluxR[3][2] - cdU[3]);
}

#endif /* SWIFT_GEAR_RT_RIEMANN_GLF_H */