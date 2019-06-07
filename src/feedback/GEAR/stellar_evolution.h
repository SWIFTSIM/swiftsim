/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STELLAR_EVOLUTION_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_GEAR_H

#include "stellar_evolution_struct.h"
#include "lifetime.h"
#include "initial_mass_function.h"
#include "supernovae_ia.h"
#include "supernovae_ii.h"

#include <math.h>
#include <stddef.h>


/**
 * @brief Initialize the global properties of the stellar evolution scheme.
 *
 * @param sm The #stellar_model.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void stellar_evolution_props_init(
    struct stellar_model* sm, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct cosmology* cosmo) {

  /* Initialize the initial mass function */
  stellar_evolution_init_initial_mass_function(&sm->imf, phys_const, us, params);

  /* Initialize the lifetime model */
  stellar_evolution_init_lifetime(&sm->lifetime, phys_const, us, params);

  /* Initialize the supernovae Ia model */
  stellar_evolution_init_supernovae_ia(&sm->snia, phys_const, us, params, &sm->imf);
 
  /* Initialize the supernovae II model */
  stellar_evolution_init_supernovae_ii(&sm->snii, phys_const, us, params, &sm->imf);


}

#endif // SWIFT_STELLAR_EVOLUTION_GEAR_H
