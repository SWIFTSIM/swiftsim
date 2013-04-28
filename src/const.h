/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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


/* Hydrodynamical constants. */
#define const_gamma             (5.0f/3.0f)
#define const_viscosity_alpha   0.8f

/* Time integration constants. */
#define const_cfl               0.3f
#define const_ln_max_h_change   log(1.26f)    /* Particle can't change volume by more than a factor of 2=1.26^3 over one time step */

/* Neighbour search constants. */
#define const_eta_kernel        1.1272f       /* Corresponds to 48 ngbs with the cubic spline kernel */
#define const_delta_nwneigh     1.f
#define CUBIC_SPLINE_KERNEL
