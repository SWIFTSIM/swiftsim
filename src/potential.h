/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016  Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_POTENTIAL_H
#define SWIFT_POTENTIAL_H

/**
 * @file src/potential.h
 * @brief Branches between the different external gravitational potentials.
 */

/* Config parameters. */
#include "../config.h"

/* Import the right external potential definition */
#if defined(EXTERNAL_POTENTIAL_NONE)
#include "./potential/none/potential.h"
#elif defined(EXTERNAL_POTENTIAL_POINTMASS)
#include "./potential/point_mass/potential.h"
#elif defined(EXTERNAL_POTENTIAL_ISOTHERMAL)
#include "./potential/isothermal/potential.h"
#elif defined(EXTERNAL_POTENTIAL_DISC_PATCH)
#include "./potential/disc_patch/potential.h"
#else
#error "Invalid choice of external potential"
#endif

/* Now, some generic functions, defined in the source file */
void potential_init(const struct swift_params* parameter_file,
                    const struct phys_const* phys_const,
                    const struct UnitSystem* us, const struct space* s,
                    struct external_potential* potential);

void potential_print(const struct external_potential* potential);

#endif /* SWIFT_POTENTIAL_H */
