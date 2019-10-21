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
#include "feedback_logger_struct.h"
#include "feedback_properties.h"
#include "feedback_logger_SNIa.h"
#include "feedback_logger_SNII.h"
#include "feedback_logger_r_processes.h"

/**
 * @brief Initialize the times of the feedback logger
 *
 * @param e the engine on this node
 */
INLINE static void feedback_logger_init(const struct engine *restrict e) {
  feedback_logger_SNII_init(e);
  feedback_logger_SNIa_init(e);
  feedback_logger_r_processes_init(e);
}

INLINE static void feedback_logger_init_log_file(const struct engine *restrict e) {
  feedback_logger_SNII_init_log_file(e); 
  feedback_logger_SNIa_init_log_file(e);
  feedback_logger_r_processes_init_log_file(e);
}

INLINE static void feedback_logger_log_data(const struct engine *restrict e) {
  feedback_logger_SNII_log_data(e);
  feedback_logger_SNIa_log_data(e);
  feedback_logger_r_processes_log_data(e);
}

INLINE static void feedback_logger_open_files(void) {
  log_SNII.core.fp = fopen("SNII.txt", "w");
  log_SNIa.core.fp = fopen("SNIa.txt", "w");
  log_r_processes.core.fp = fopen("r_processes.txt", "w");
}

INLINE static void feedback_logger_time_step(const struct engine *restrict e) {
  feedback_logger_SNII_time_step(e);
  feedback_logger_SNIa_time_step(e);
  feedback_logger_r_processes_time_step(e);
}


INLINE static void feedback_logger_close(const struct engine *restrict e) {
  fclose(log_SNII.core.fp);
  fclose(log_SNIa.core.fp);
  fclose(log_r_processes.core.fp);
}

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
INLINE static void feedback_logger_MPI(const struct engine *restrict e) {
  
  feedback_logger_SNII_MPI(e);
  feedback_logger_SNIa_MPI(e);
  feedback_logger_r_processes_MPI(e);

}
#endif

/**
 * @brief Log the event in case of debugging
 *
 * @param fp the file pointer to write the debugging information to
 * @param time the current simulation time
 * @param si the star particle
 * @param pj the gas particle
 * @param xpj the extra information of the gas particle
 * @param cosmo the cosmology struct
 * @param step the current simulation step
 */
INLINE static void feedback_logger_SNIa_log_event_debug(
    FILE *fp, const double time, const struct spart *restrict si,
    struct part *restrict pj, struct xpart *restrict xpj,
    const struct cosmology *restrict cosmo, const int step) {}

/**
 * @brief Initialize the SNIa logger debug file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_logger_SNIa_init_log_file_debug(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {}

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_H */
