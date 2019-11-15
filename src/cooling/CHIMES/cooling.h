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
#ifndef SWIFT_COOLING_CHIMES_H
#define SWIFT_COOLING_CHIMES_H

/**
 * @file src/cooling/CHIMES/cooling.h
 * @brief Empty infrastructure for CHIMES cooling. 
 */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "cooling_struct.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct space;

void cooling_init_backend(struct swift_params* parameter_file,
			  const struct unit_system* us,
			  const struct phys_const* phys_const,
                          const struct hydro_props *hydro_props,
			  struct cooling_function_data* cooling);  

void cooling_print_backend(const struct cooling_function_data *cooling); 

void cooling_first_init_part(const struct phys_const* restrict phys_const,
			     const struct unit_system* restrict us,
			     const struct cosmology* restrict cosmo,
			     const struct cooling_function_data* data, 
			     struct part* restrict p,
			     struct xpart* restrict xp); 

void chimes_set_init_eqm(const struct phys_const* restrict phys_const,
			 const struct unit_system* restrict us,
			 const struct cosmology* restrict cosmo,
			 const struct hydro_props *hydro_properties,
			 const struct entropy_floor_properties *floor_props,
			 const struct cooling_function_data* data, 
			 struct part* restrict p,
			 struct xpart* restrict xp); 

void cooling_update(const struct cosmology* cosmo,
		    struct cooling_function_data* cooling,
		    struct space* s); 

void chimes_update_gas_vars(const double u_cgs,
			    const struct phys_const *phys_const,
			    const struct unit_system *us,
			    const struct cosmology *cosmo,
			    const struct hydro_props *hydro_properties,
			    const struct entropy_floor_properties *floor_props,
			    const struct cooling_function_data *cooling,
			    struct part *restrict p, struct xpart *restrict xp,
			    struct gasVariables *ChimesGasVars, 
			    const float dt_cgs); 

void chimes_update_element_abundances(const struct phys_const *phys_const,
				      const struct unit_system *us,
				      const struct cosmology *cosmo,
				      const struct cooling_function_data *cooling,
				      struct part *restrict p, 
				      struct xpart *restrict xp,
				      struct gasVariables *ChimesGasVars, 
				      const int mode);

void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm, const double time); 

float cooling_get_radiated_energy(const struct xpart* restrict xp); 

double chimes_mu(const struct cooling_function_data *cooling,
		 const struct part *restrict p, 
		 const struct xpart *restrict xp); 

float cooling_get_temperature(const struct phys_const* restrict phys_const,
			      const struct hydro_props* restrict hydro_props,
			      const struct unit_system* restrict us,
			      const struct cosmology* restrict cosmo,
			      const struct cooling_function_data* restrict cooling,
			      const struct part* restrict p, 
			      const struct xpart* restrict xp);

double calculate_colibre_N_ref(const struct phys_const *phys_const,
			       const struct unit_system *us,
			       const struct cosmology *cosmo,
			       const struct cooling_function_data *cooling,
			       struct part *restrict p, struct xpart* restrict xp, 
			       const double mu); 

void cooling_struct_dump(const struct cooling_function_data* cooling, 
			 FILE* stream); \

void cooling_struct_restore(struct cooling_function_data* cooling,
			    FILE* stream, const struct cosmology* cosmo); 

void cooling_convert_quantities(struct part *restrict p, 
				struct xpart *restrict xp,
				const struct cosmology *cosmo, 
				const struct hydro_props *hydro_props, 
				const struct phys_const *phys_const, 
				const struct unit_system* us, 
				const struct entropy_floor_properties *floor_props, 
				const struct cooling_function_data* cooling); 

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended data of the particle.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props, const struct part* restrict p,
    const struct xpart* restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
static INLINE void cooling_clean(struct cooling_function_data* cooling) {}

#endif /* SWIFT_COOLING_CHIMES_H */
