/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
 *               
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

#ifndef SWIFT_MINIMAL_HYDRO_DEFAULT_H
#define SWIFT_MINIMAL_HYDRO_DEFAULT_H

/**
 * @file Minimal/hydro_defaults.h
 * @brief Minimal conservative implementation of SPH . (default parameters)
 *        
 *        This file defines a number of things that are used in
 *        hydro_properties.c as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 */

/* Viscosity parameters -- FIXED -- MUST BE DEFINED AT COMPILE-TIME */

/* Cosmology default beta=3.0. Planetary default beta=4.0
 * Alpha can be set in the parameter file.
 * Beta is defined as in e.g. Price (2010) Eqn (103) */
#define const_viscosity_beta 3.0f

/* The viscosity that the particles are reset to after being hit by a
 * feedback event. This should be set to the same value as the
 * hydro_props_default_viscosity_alpha in fixed schemes, and likely
 * to hydro_props_default_viscosity_alpha_max in variable schemes. */
#define hydro_props_default_viscosity_alpha_feedback_reset 0.8f


/* Viscosity paramaters -- Defaults; can be changed at run-time */

/* The "initial" hydro viscosity, or the fixed value for non-variable
 * schemes. This usually takes the value 0.8. */
#define hydro_props_default_viscosity_alpha 0.8f

/* Minimal value for the viscosity alpha in variable schemes. In 
 * non-variable schemes this must be defined but is not used. */
#define hydro_props_default_viscosity_alpha_min 0.8f

/* Maximal value for the viscosity alpha in variable schemes. In 
 * non-variable schemes this must be defined but is not used. */
#define hydro_props_default_viscosity_alpha_max 0.8f

/* Decay length for the viscosity scheme. This is scheme dependent. In
 * non-variable schemes this must be defined but is not used. This also
 * sets the decay length for the diffusion. */
#define hydro_props_default_viscosity_length 0.25f


/* Diffusion parameters -- Defaults; can be changed at run-time */

/* The "initial" diffusion, or the fixed value for non-variable
 * schemes. This usually takes the value 0.0. */
#define hydro_props_default_diffusion_alpha 0.0f

/* Beta coefficient for the diffusion. This controls how fast the
 * diffusion coefficient peaks, and how high it can get. Chosen to be
 * very small in schemes where little diffusion is needed, 0.2-1.0 in
 * schemes (e.g. density-energy) where diffusion is needed to solve
 * the contact discontinuity problem. */
#define hydro_props_default_diffusion_beta 0.0f

/* Maximal value for the diffusion alpha in variable schemes. In 
 * non-variable schemes this must be defined but is not used. */
#define hydro_props_default_diffusion_alpha_max 0.0f

/* Minimal value for the diffusion alpha in variable schemes. In 
 * non-variable schemes this must be defined but is not used. */
#define hydro_props_default_diffusion_alpha_min 0.0f

#endif /* SWIFT_MINIMAL_HYDRO_DEFAULT_H */