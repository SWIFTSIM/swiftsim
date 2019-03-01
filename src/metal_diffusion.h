/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Camila Correa (correa@strw.leidenuniv.nl)
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
#ifndef SWIFT_METAL_DIFFUSION_H
#define SWIFT_METAL_DIFFUSION_H

/**
 * @file src/.?./diffusion.h
 * @brief Branches between the two diffusion functions:
 * 1) Diffusion of pasive scalar diffusion (useful for 2D hydrodynamic tests),
 * 2) Diffusion of element abundance
 */

/* Config parameters. */
#include "../config.h"

/* Import the right diffusion definition */
#if defined(DIFFUSION_NONE)
#include "./metal_diffusion/none/diffusion.h"
#include "./metal_diffusion/none/diffusion_iact.h"
#elif defined(SCALAR_DIFFUSION)
#include "./metal_diffusion/DiffusionOfScalar/diffusion.h"
#include "./metal_diffusion/DiffusionOfScalar/diffusion_iact.h"
#elif defined(METAL_DIFFUSION)
#include "./metal_diffusion/DiffusionOfElementAbundance/diffusion.h"
#include "./metal_diffusion/DiffusionOfElementAbundance/diffusion_iact.h"
#else
#error "Invalid choice of diffusion function."
#endif

/* Common functions */
void diffusion_init(struct swift_params* parameter_file,
                    struct diffusion_global_data* data);

void diffusion_print(const struct diffusion_global_data* data);

/* Dump/restore. */
void diffusion_struct_dump(const struct diffusion_global_data* diffusion,
                           FILE* stream);
void diffusion_struct_restore(const struct diffusion_global_data* diffusion,
                              FILE* stream);


#endif
