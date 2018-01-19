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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "chemistry.h"

/**
 * @brief Initialises the chemistry properties.
 *
 * Calls chemistry_init_backend for the chosen chemistry function.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
void chemistry_init(const struct swift_params* parameter_file,
                    const struct unit_system* us,
                    const struct phys_const* phys_const,
                    struct chemistry_data* data) {

  chemistry_init_backend(parameter_file, us, phys_const, data);
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * Calls chemistry_print_backend for the chosen chemistry model.
 *
 * @brief The #chemistry_data containing information about the current model.
 */
void chemistry_print(const struct chemistry_data* data) {
  chemistry_print_backend(data);
}
