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
#ifndef SWIFT_COLIBRE_EVENT_LOGGER_H
#define SWIFT_COLIBRE_EVENT_LOGGER_H

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local includes */
#include "engine.h"
#include "event_logger_SNII.h"
#include "event_logger_SNIa.h"
#include "event_logger_r_processes.h"

/**
 * @brief Initialize the times of the feedback logger
 *
 * @param e the engine on this node
 */
INLINE static void event_logger_init(const struct engine *e) {
  event_logger_SNII_init(e);
  event_logger_SNIa_init(e);
  event_logger_r_processes_init(e);

#ifdef SWIFT_DEBUG_CHECKS
  event_logger_SNIa_init_debug(e);
  event_logger_SNII_init_debug(e);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief Initialize the log files
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_init_log_file(const struct engine *e) {
  event_logger_SNII_init_log_file(e);
  event_logger_SNIa_init_log_file(e);
  event_logger_r_processes_init_log_file(e);

#ifdef SWIFT_DEBUG_CHECKS
  event_logger_SNIa_init_log_file_debug(e);
  event_logger_SNII_init_log_file_debug(e);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief Write the data to the log file if we have the correct time, do this
 * for all the log files
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_log_data(const struct engine *e) {
  event_logger_SNII_log_data(e);
  event_logger_SNIa_log_data(e);
  event_logger_r_processes_log_data(e);
}

/**
 * @brief Open the files of the feedback loggers
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_open_files(const struct engine *e,
                                           const char *mode) {
  log_SNII.core.fp = fopen("SNII.txt", mode);
  log_SNIa.core.fp = fopen("SNIa.txt", mode);
  log_r_processes.core.fp = fopen("r_processes.txt", mode);

#ifdef SWIFT_DEBUG_CHECKS
  /* Open SNIa debugging file */
  char savename_SNIa[50];
  snprintf(savename_SNIa, 50, "SNIa_%d.txt", e->nodeID);
  log_SNIa_debug.fp = fopen(savename_SNIa, mode);

  /* Open SNII debugging file */
  char savename_SNII[50];
  snprintf(savename_SNII, 50, "SNII_%d.txt", e->nodeID);
  log_SNII_debug.fp = fopen(savename_SNII, mode);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief Update the core values in the logger for the time step
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_time_step(const struct engine *e) {
  event_logger_SNII_time_step(e);
  event_logger_SNIa_time_step(e);
  event_logger_r_processes_time_step(e);
}

/**
 * @brief Close the files for the feedback logger
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_close(const struct engine *e) {
  event_logger_SNII_log_data_end(e);
  event_logger_SNIa_log_data_end(e);
  event_logger_r_processes_log_data_end(e);
}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication between all the nodes regarding the feedback
 * logger
 *
 * @param e the engine we are running on
 */
INLINE static void event_logger_MPI_Reduce(const struct engine *e) {

  event_logger_SNII_MPI_Reduce(e);
  event_logger_SNIa_MPI_Reduce(e);
  event_logger_r_processes_MPI_Reduce(e);
}
#endif

/**
 * @brief dump the feedback logger info to the restart file
 *
 * @param stream, the data stream
 */
INLINE static void event_logger_struct_dump(FILE *stream) {
  restart_write_blocks((void *)&log_SNII, sizeof(struct feedback_history_SNII),
                       1, stream, "logger SNII", "Event logger for SNII");
  restart_write_blocks((void *)&log_SNIa, sizeof(struct feedback_history_SNIa),
                       1, stream, "logger SNIa", "Event logger for SNIa");
  restart_write_blocks((void *)&log_r_processes,
                       sizeof(struct feedback_history_r_processes), 1, stream,
                       "logger r-proc", "Event logger for r-processes");
#ifdef SWIFT_DEBUG_CHECKS
  restart_write_blocks(
      (void *)&log_SNII_debug, sizeof(struct feedback_history_debug), 1, stream,
      "logger SNII debug", "Event logger for SNII during debugging");
  restart_write_blocks(
      (void *)&log_SNIa_debug, sizeof(struct feedback_history_debug), 1, stream,
      "logger SNIa debug", "Event logger for SNIa during debugging");
#endif /* SWIFT_DEBUG_CHECKS*/
}

/**
 * @brief restore the feedback logger info to the restart file
 *
 * @param stream, the data stream
 */
INLINE static void event_logger_struct_restore(FILE *stream) {
  restart_read_blocks((void *)&log_SNII, sizeof(struct feedback_history_SNII),
                      1, stream, NULL, "Event logger for SNII");
  restart_read_blocks((void *)&log_SNIa, sizeof(struct feedback_history_SNIa),
                      1, stream, NULL, "Event logger for SNIa");
  restart_read_blocks((void *)&log_r_processes,
                      sizeof(struct feedback_history_r_processes), 1, stream,
                      NULL, "Event logger for r-processes");
#ifdef SWIFT_DEBUG_CHECKS
  restart_read_blocks((void *)&log_SNII_debug,
                      sizeof(struct feedback_history_debug), 1, stream, NULL,
                      "Event logger for SNII during debugging");
  restart_read_blocks((void *)&log_SNIa_debug,
                      sizeof(struct feedback_history_debug), 1, stream, NULL,
                      "Event logger for SNIa during debugging");
#endif /* SWIFT_DEBUG_CHECKS*/
}

#endif /* SWIFT_COLIBRE_EVENT_LOGGER_H */
