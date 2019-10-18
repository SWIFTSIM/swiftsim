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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_H

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "feedback_logger_struct.h"
#include "feedback_properties.h"
#include "engine.h"

/**
 * @brief Initialize the times of the feedback logger
 *
 * @param e the engine on this node
 */
INLINE static void feedback_logger_init(const struct engine e) {}

INLINE static void feedback_logger_init_log_file(const struct engine e) {}

INLINE static void feedback_logger_log_data(const struct engine e) {}

INLINE static void feedback_logger_open_files(void) {}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication between all the nodes regarding the feedback logger
 *
 * @param nodeID the nodeID
 * @param SNII the SNII feedback structure
 * @param SNIa the SNIa feedback structure
 * @param r_processes the r-processes feedback structure
 * @param logger_time the current time
 * @param delta_logger_time the delta time stepping of the feedback logger
 * @return the new time variable used to log the feedback data
 */
INLINE static double feedback_logger_MPI(const struct engine e) {}
#endif




#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_H */
