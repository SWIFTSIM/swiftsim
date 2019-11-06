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
#ifndef SWIFT_NONE_EVENT_LOGGER_H
#define SWIFT_NONE_EVENT_LOGGER_H

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "engine.h"

/**
 * @brief Initialize the times of the feedback logger
 *
 * @param e the engine on this node
 */
INLINE static void event_logger_init(const struct engine *restrict e) {}

/**
 * @brief Initialize the log files
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_init_log_file(
    const struct engine *restrict e) {}

/**
 * @brief Write the data to the log file if we have the correct time, do this
 * for all the log files
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_log_data(const struct engine *restrict e) {}

/**
 * @brief Open the files of the feedback loggers
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_open_files(const struct engine *restrict e) {
}

/**
 * @brief Update the core values in the logger for the time step
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_time_step(const struct engine *restrict e) {}

/**
 * @brief Close the files for the feedback logger
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_close(const struct engine *restrict e) {}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication between all the nodes regarding the feedback
 * logger
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_MPI_Reduce(const struct engine *restrict e) {
}
#endif

/**
 * @brief dump the feedback logger info to the restart file
 *
 * @param stream, the data stream
 */
INLINE static void event_logger_struct_dump(FILE *stream) {}

/**
 * @brief restore the feedback logger info to the restart file
 *
 * @param stream, the data stream
 */
INLINE static void event_logger_struct_restore(FILE *stream) {}

#endif /* SWIFT_NONE_EVENT_LOGGER_H */
