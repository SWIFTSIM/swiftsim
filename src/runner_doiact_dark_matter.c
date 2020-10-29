/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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

/* Local headers. */
#include "active.h"
#include "dark_matter.h"
#include "dark_matter_iact.h"
#include "dark_matter_logger.h"
#include "cell.h"
#include "engine.h"
#include "runner.h"
#include "space_getsid.h"
#include "timers.h"

/* Import the dark matter density loop functions. */
#define FUNCTION dark_matter_density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_functions_dark_matter.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

