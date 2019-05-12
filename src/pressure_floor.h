/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PRESSURE_FLOOR_H
#define SWIFT_PRESSURE_FLOOR_H

/**
 * @file src/pressure_floor.h
 * @brief Branches between the different pressure floor models
 */

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "common_io.h"
#include "cosmology.h"
#include "error.h"
#include "inline.h"


extern struct pressure_floor_properties pressure_floor_props;

/* Pre-declarations to avoid cyclic inclusions */
static INLINE float hydro_get_physical_density(const struct part *restrict p,
                                               const struct cosmology *cosmo);
static INLINE float hydro_get_comoving_density(const struct part *restrict p);

static INLINE float pressure_floor_get_physical_pressure(
    const struct part *p, const struct cosmology *cosmo,
    const float pressure);

static INLINE float pressure_floor_get_comoving_pressure(
    const struct part *p, const float pressure);

/* Check if pressure floor is implemented in hydro */
#ifndef PRESSURE_FLOOR_NONE
#ifdef GADGET2_SPH
/* Implemented */
#else
#error Pressure floor not implemented with this hydro scheme
#endif

#endif
/* Import the right pressure floor definition */
#if defined(PRESSURE_FLOOR_NONE)
#include "./pressure_floor/none/pressure_floor.h"
#elif defined(PRESSURE_FLOOR_GEAR)
#include "./pressure_floor/GEAR/pressure_floor.h"
#endif

#endif /* SWIFT_PRESSURE_FLOOR_H */
