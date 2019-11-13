/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#include "config.h"

/* This file's header */
#include "event_logger.h"

#ifdef SWIFT_DEBUG_CHECKS
struct feedback_history_debug log_SNIa_debug;
struct feedback_history_debug log_SNII_debug;
#endif

/* Declare the feedback structures */

#ifdef __APPLE__
/*
 * The clang compiler and linker on OSX incorrectly optimize
 * out the logger global objects before the final linking stage, which
 * leads to a compilation errors.
 * The fake initialisations below forces the compiler to keep the
 * instance and pass it to the linker stage.
 */
struct feedback_history_SNIa log_SNIa = {.core.fp = NULL};
struct feedback_history_SNII log_SNII = {.core.fp = NULL};
struct feedback_history_r_processes log_r_processes = {.core.fp = NULL};

#else  /* i.e. not __APPLE__ */
struct feedback_history_SNIa log_SNIa;
struct feedback_history_SNII log_SNII;
struct feedback_history_r_processes log_r_processes;
#endif /* __APPLE__ */
