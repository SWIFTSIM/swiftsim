/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_DEBUG_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_DEBUG_H

/**
 * @file PressureEnergy/hydro_debug.h
 * @brief Pressure-Energy implementation of SPH (Debugging routines)
 *
 * The thermal variable is the energy (u) and the pressure is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows equations (16), (17) and (18) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const struct part* p, const struct xpart* xp) {
  printf(
      "x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e],v_full=[%.3e,%.3e,%.3e] \n a=[%.3e,%.3e,%.3e],\n "
      "h=%.3e, wcount=%.3f, wcount_dh=%.3e, time_bin=%d\n"
      "p_bar=%.3e, p_bar_dh=%3.e, u=%3.e\n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], xp->v_full[0],
      xp->v_full[1], xp->v_full[2], p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->h, p->density.wcount, p->density.wcount_dh, p->time_bin,
      p->pressure_bar, p->pressure_bar_dh, p->u);
}

#endif /* SWIFT_PRESSURE_ENERGY_HYDRO_DEBUG_H */
