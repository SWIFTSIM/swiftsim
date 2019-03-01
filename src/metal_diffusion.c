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
#include "metal_diffusion.h"

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

/**
 * @brief Prints the properties of the diffusion model to stdout.
 *
 * Calls diffusion_print_backend for the chosen diffusion model.
 *
 * @brief The #diffusion_global_data containing information about the current
 * model.
 */
void diffusion_print(const struct diffusion_global_data* data) {
    diffusion_print_backend(data);
}

/**
 * @brief Write a diffusion struct to the given FILE as a stream of bytes.
 *
 * @param diffusion the struct
 * @param stream the file stream
 */
void diffusion_struct_dump(const struct diffusion_global_data* diffusion,
                           FILE* stream) {
    restart_write_blocks((void*)diffusion, sizeof(struct diffusion_global_data),
                         1, stream, "diffusion", "diffusion function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param diffusion the struct
 * @param stream the file stream
 */
void diffusion_struct_restore(const struct diffusion_global_data* diffusion,
                              FILE* stream) {
    restart_read_blocks((void*)diffusion, sizeof(struct diffusion_global_data), 1,
                        stream, NULL, "diffusion function");
}
