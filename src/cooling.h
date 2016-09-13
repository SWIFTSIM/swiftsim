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
#ifndef SWIFT_COOLING_H
#define SWIFT_COOLING_H

/**
 * @file src/cooling.h
 * @brief Branches between the different cooling functions.
 */

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"

#include "cooling/cooling_models.h"

/* Import the right cooling definition */
//#include "./cooling/none/cooling.h"

#include "./cooling/const_du/cooling.h"
#include "./cooling/const_lambda/cooling.h"
//#elif defined(COOLING_GRACKLE)
//#include "./cooling/grackle/cooling.h"
//#else
//#error "Invalid choice of cooling function."
//#endif

/* Common functions */
void cooling_init(const struct swift_params* parameter_file,
                  const struct UnitSystem* us,
                  const struct phys_const* phys_const,
                  cooling_function_data_handle* cooling_handle_ptr);

void cooling_print(const cooling_function_data_handle cooling_handle);

#endif /* SWIFT_COOLING_H */
