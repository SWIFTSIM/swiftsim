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

#include "engine.h"
#include "feedback_logger_SNII.h"
#include "feedback_logger_SNIa.h"
#include "feedback_logger_r_processes.h"
#include "feedback_logger_struct.h"
#include "feedback_properties.h"

/**
 * @brief Initialize the times of the feedback logger
 *
 * @param e the engine on this node
 */
INLINE static void feedback_logger_init(const struct engine *restrict e) {
  feedback_logger_SNII_init(e);
  feedback_logger_SNIa_init(e);
  feedback_logger_r_processes_init(e);

#ifdef SWIFT_DEBUG_CHECKS
  feedback_logger_SNIa_init_debug(e);
  feedback_logger_SNII_init_debug(e);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief Initialize the log files
 *
 * @param e the engine we are running on
 */
INLINE static void feedback_logger_init_log_file(
    const struct engine *restrict e) {
  feedback_logger_SNII_init_log_file(e);
  feedback_logger_SNIa_init_log_file(e);
  feedback_logger_r_processes_init_log_file(e);

#ifdef SWIFT_DEBUG_CHECKS
  feedback_logger_SNIa_init_log_file_debug(e);
  feedback_logger_SNII_init_log_file_debug(e);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief Write the data to the log file if we have the correct time, do this
 * for all the log files
 *
 * @param e the engine we are running on
 */
INLINE static void feedback_logger_log_data(const struct engine *restrict e) {
  feedback_logger_SNII_log_data(e);
  feedback_logger_SNIa_log_data(e);
  feedback_logger_r_processes_log_data(e);
}

/**
 * @brief Open the files of the feedback loggers
 *
 * @param e the engine we are running on
 */
INLINE static void feedback_logger_open_files(const struct engine *restrict e) {
  log_SNII.core.fp = fopen("SNII.txt", "w");
  log_SNIa.core.fp = fopen("SNIa.txt", "w");
  log_r_processes.core.fp = fopen("r_processes.txt", "w");

#ifdef SWIFT_DEBUG_CHECKS
  /* Open SNIa debugging file */
  char savename_SNIa[50];
  snprintf(savename_SNIa, 50, "SNIa_%d.txt", e->nodeID);
  log_SNIa_debug.fp = fopen(savename_SNIa, "w");

  /* Open SNII debugging file */
  char savename_SNII[50];
  snprintf(savename_SNII, 50, "SNII_%d.txt", e->nodeID);
  log_SNII_debug.fp = fopen(savename_SNII, "w");
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief Update the core values in the logger for the time step
 *
 * @param e the engine we are running on
 */
INLINE static void feedback_logger_time_step(const struct engine *restrict e) {
  feedback_logger_SNII_time_step(e);
  feedback_logger_SNIa_time_step(e);
  feedback_logger_r_processes_time_step(e);
}

/**
 * @brief Close the files for the feedback logger
 *
 * @param e the engine we are running on
 */
INLINE static void feedback_logger_close(const struct engine *restrict e) {
  feedback_logger_SNII_log_data_end(e);
  feedback_logger_SNIa_log_data_end(e);
  feedback_logger_r_processes_log_data_end(e);
}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication between all the nodes regarding the feedback
 * logger
 *
 * @param e the engine we are running on
 */
INLINE static void feedback_logger_MPI_Reduce(const struct engine *restrict e) {

  feedback_logger_SNII_MPI_Reduce(e);
  feedback_logger_SNIa_MPI_Reduce(e);
  feedback_logger_r_processes_MPI_Reduce(e);
}
#endif

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_H */
