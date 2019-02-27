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
/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "diffusion.h"

/**
 * @brief Initialises the passive scalar arrays.
 *
 * Calls diffusion_init_backend for the chosen diffusion function.
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.  */
void diffusion_init(struct swift_params* parameter_file,
                    struct diffusion_global_data* data) {
    
    diffusion_init_backend(parameter_file, data);
}
