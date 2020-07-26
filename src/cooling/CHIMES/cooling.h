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
#include "dust.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct space;

void cooling_init_backend(struct swift_params *parameter_file,
			    const struct unit_system *us,
			    const struct phys_const *phys_const,
			    const struct hydro_props *hydro_props,
			    struct cooling_function_data *cooling,
			    struct dustevo_props *dp);

void cooling_print_backend(const struct cooling_function_data *cooling);

void cooling_update(const struct cosmology *cosmo,
                    struct cooling_function_data *cooling, struct space *s);

void chimes_update_gas_vars(const double u_cgs,
                            const struct phys_const *phys_const,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct hydro_props *hydro_properties,
                            const struct entropy_floor_properties *floor_props,
                            const struct cooling_function_data *cooling,
			    const struct dustevo_props *dp,
                            struct part *restrict p, struct xpart *restrict xp,
                            struct gasVariables *ChimesGasVars,
                            const float dt_cgs);

void chimes_update_element_abundances(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct cooling_function_data *cooling,
    struct part *restrict p, struct xpart *restrict xp,
    struct gasVariables *ChimesGasVars, const int mode);

void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
		       const struct dustevo_props *dp,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm, const double time);

float cooling_get_radiated_energy(const struct xpart *restrict xp);

void cooling_split_part(struct part *p, struct xpart *xp, double n);

double chimes_mu(const struct cooling_function_data *cooling,
                 const struct part *restrict p,
                 const struct xpart *restrict xp);

float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp);

double calculate_colibre_N_ref(const struct phys_const *phys_const,
                               const struct unit_system *us,
                               const struct cosmology *cosmo,
                               const struct cooling_function_data *cooling,
                               struct part *restrict p,
                               struct xpart *restrict xp, const double mu);

void cooling_struct_dump(const struct cooling_function_data *cooling,
			 struct dustevo_props *dp,
                         FILE *stream);

void cooling_struct_restore(struct cooling_function_data *cooling,
			    struct dustevo_props* dp, FILE *stream,
                            const struct cosmology *cosmo);

void cooling_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling,
    const struct dustevo_props *dp);

void cooling_set_subgrid_properties(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, 
    const struct dustevo_props *dp, struct part *p,
    struct xpart *xp);

float cooling_get_subgrid_temperature(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

float cooling_get_subgrid_density(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

void cooling_set_HIIregion_chimes_abundances(
    struct gasVariables *ChimesGasVars,
    const struct cooling_function_data *cooling);

void cooling_set_FB_particle_chimes_abundances(
    struct gasVariables *ChimesGasVars,
    const struct cooling_function_data *cooling);

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
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props, const struct part *restrict p,
    const struct xpart *restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * @param cooling the cooling data structure.
 */
static INLINE void cooling_clean(struct cooling_function_data *cooling) {}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * Initialises the cooling_data.heated_by_FB flag to zero.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo The current cosmological model.
 * @param data The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *data, struct part *restrict p,
    struct xpart *restrict xp) {
  xp->cooling_data.heated_by_FB = 0;
}

/**
 * @brief Updates cooling properties of particle hit by feedback.
 *
 * Sets the cooling_data.heated_by_FB flag to 1.
 *
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
cooling_update_feedback_particle(struct xpart *restrict xp) {
  xp->cooling_data.heated_by_FB = 1;
}
#endif /* SWIFT_COOLING_CHIMES_H */
