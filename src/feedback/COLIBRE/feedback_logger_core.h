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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_CORE_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_CORE_H

INLINE static void feedback_logger_core_init(const struct engine *restrict e, struct feedback_history_logger *restrict fhl, const int log_type) {
  
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
  fhl->logger_time = 0.;

  /* Make a constant for the physical constants */
  const struct phys_const *phys_const = e->physical_constants;

  /* Define the swift parameter file */
  struct swift_params *params = e->parameter_file;

  /* Define the delta logger time */
  double delta_logger_time_Myr;

  /* Make a difference between the 3 types of logger */
  switch (log_type) {
    case 1:
      /* Time interval to use to store the SNIa events */
      delta_logger_time_Myr = parser_get_param_double(
      params, "Feedback_logger:delta_time_SNII_Myr");
      break;
    case 2:
      /* Time interval to use to store the SNIa events */
      delta_logger_time_Myr = parser_get_param_double(
      params, "Feedback_logger:delta_time_SNIa_Myr");
      break;
    case 3:
      /* Time interval to use to store the SNIa events */
      delta_logger_time_Myr = parser_get_param_double(
      params, "Feedback_logger:delta_time_r_processes_Myr");
      break;
  }

  fhl->delta_logger_time = delta_logger_time_Myr * 1e6 * phys_const->const_year;
}

INLINE static void feedback_logger_core_time_step(const struct engine *restrict e, struct feedback_history_logger *restrict fhl) {
  
  fhl->logger_time += e->time_step;
}

INLINE static int feedback_logger_core_log(const struct engine *restrict e, struct feedback_history_logger *restrict fhl) {
  
  if (fhl->logger_time < fhl->delta_logger_time) return 0;

  return 1;
}

INLINE static void feedback_logger_core_update(const struct engine *restrict e, struct feedback_history_logger *restrict fhl) {
  
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
  fhl->logger_time -= fhl->delta_logger_time;

}


#endif
