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
#ifndef SWIFT_EVENT_LOGGER_STRUCT_H
#define SWIFT_EVENT_LOGGER_STRUCT_H

/* feedback history general struct for any quantity */
struct event_history_logger {

  /* Feedback history lock */
  swift_lock_type lock;

  /* Feedback history file pointer */
  FILE *fp;

  /* Feedback history logger time since the last log to the file */
  double logger_time_since_last_log;

  /* Feedback history delta logger time */
  double delta_logger_time;

  /*! Previous storage step*/
  int step_prev;

  /*! Previous storage time */
  float time_prev;

  /*! Previous storage redshift */
  float z_prev;

  /*! Previous storage scale factor */
  float a_prev;
};

/**
 * @brief Initialize the logger core
 *
 * @param e the engine we are running
 * @param fhl the feedback history logger that is the core
 * @param log_type the log type to initialize
 */
INLINE static void event_logger_core_init(
    const struct engine *restrict e,
    struct event_history_logger *restrict fhl) {

  /* Initialize the lock*/
  lock_init(&fhl->lock);

  /* initialization of the step */
  fhl->step_prev = 0;

  /* Initialize the first time */
  fhl->time_prev = e->time;

  /* Initialize the redshift */
  fhl->z_prev = e->cosmology->z;

  /* Initialize the scale factor */
  fhl->a_prev = e->cosmology->a;

  /* Initialize the logger time */
  fhl->logger_time_since_last_log = 0.;
}

/**
 * @brief update the current time with the timestep in the logger core
 *
 * @param e the engine we are running
 * @param fhl the feedback history logger that is the core
 */
INLINE static void event_logger_core_time_step(
    const struct engine *restrict e,
    struct event_history_logger *restrict fhl) {

  fhl->logger_time_since_last_log += e->time_step;
}

/**
 * @brief Check if we need to log for this core
 *
 * @param e the engine we are running
 * @param fhl the feedback history logger that is the core
 */
INLINE static int event_logger_core_log(const struct engine *restrict e,
                                           struct event_history_logger *restrict
                                               fhl) {

  return (fhl->logger_time_since_last_log >= fhl->delta_logger_time);
}

/**
 * @brief Update the logger core values after logging
 *
 * @param e the engine we are running
 * @param fhl the feedback history logger that is the core
 */
INLINE static void event_logger_core_update(
    const struct engine *restrict e,
    struct event_history_logger *restrict fhl) {

  /* Update the core values */

  /* The step number */
  fhl->step_prev = e->step;

  /* Update the time */
  fhl->time_prev = e->time;

  /* Update the scale factor */
  fhl->a_prev = e->cosmology->a;

  /* Update the redshift */
  fhl->z_prev = e->cosmology->z;

  /* Update the logger time */
  fhl->logger_time_since_last_log -= fhl->delta_logger_time;
}

#endif /* SWIFT_EVENT_LOGGER_STRUCT_H */
