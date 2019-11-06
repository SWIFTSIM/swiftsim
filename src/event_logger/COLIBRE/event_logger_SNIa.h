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
#ifndef SWIFT_COLIBRE_EVENT_LOGGER_SNIA_H
#define SWIFT_COLIBRE_EVENT_LOGGER_SNIA_H

/* Local includes */
#include "event_logger_core.h"
#include "event_logger_struct.h"
#include "feedback_properties.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* feedback history struct for SNIa */
struct feedback_history_SNIa {

  /* Load the core of logging functions */
  struct event_history_logger core;

  /*! Total new SNIa injected energy */
  double SNIa_energy;

  /*! Number of heating events */
  int events;
};

/**
 * @brief Initialize the SNIa logger file
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_init_log_file(const struct engine *e) {

  /* Load the structures of the internal units and the physical constants */
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;

  /* Use the File pointer */
  FILE *fp = log_SNIa.core.fp;

  /* Calculate the energy unit */
  const double E_unit = us->UnitMass_in_cgs * us->UnitLength_in_cgs *
                        us->UnitLength_in_cgs /
                        (us->UnitTime_in_cgs * us->UnitTime_in_cgs);

  /* Write some general text to the logger file */
  fprintf(fp, "# Stochastic SNIa logger file\n");
  fprintf(fp, "######################################################\n");
  fprintf(fp, "# The quantities are all given in internal physical units!\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# (0)  Simulation step\n");
  fprintf(fp, "# (1)  Previous simulation step\n");
  fprintf(fp,
          "# (2)  Time since Big Bang (cosmological run), Time since start of "
          "the simulation (non-cosmological run).\n");
  fprintf(fp,
          "# (3)  Previous time since Big Bang (cosmological run), Time since "
          "start of "
          "the simulation (non-cosmological run).\n");
  fprintf(fp, "#      Unit = %e seconds\n", us->UnitTime_in_cgs);
  fprintf(fp, "#      Unit = %e yr or %e Myr\n", 1.f / phys_const->const_year,
          1.f / phys_const->const_year / 1e6);
  fprintf(fp, "# (4)  Scale factor              (no unit)\n");
  fprintf(fp, "# (5)  Previous scale factor     (no unit)\n");
  fprintf(fp, "# (6)  Redshift                  (no unit)\n");
  fprintf(fp, "# (7)  Previous redshift         (no unit)\n");
  fprintf(fp, "# (8)  Injected energy of SNIa events\n");
  fprintf(fp, "#      Unit = %e erg\n", E_unit);
  fprintf(fp, "#      Unit = %e x 10^51 erg\n", E_unit / 1e51);
  fprintf(fp, "# (9)  Number of SNIa   (number, no unit)\n");
  fprintf(fp,
          "# (10) Number of SNIa per time (binned in time between current time "
          "and previous time).\n");
  fprintf(fp, "#      Unit = %e #/seconds\n", 1. / us->UnitTime_in_cgs);
  fprintf(fp, "#      Unit = %e #/yr or %e #/Myr\n", phys_const->const_year,
          phys_const->const_year * 1e6);
  fprintf(fp,
          "# (11) Number of SNIa per time per comoving volume (binned in time "
          "between current time and previous time).\n");
  fprintf(fp, "#      Unit = %e #/seconds/cm^3\n",
          1. / us->UnitTime_in_cgs / pow(us->UnitLength_in_cgs, 3));
  fprintf(fp, "#      Unit = %e #/yr/Mpc^3\n",
          phys_const->const_year * pow(phys_const->const_parsec * 1e6, 3));
  fprintf(fp, "# (12) Number of heating events   (no unit)\n");
  fprintf(fp, "#\n");
  fprintf(
      fp,
      "#  (0)      (1)         (2)              (3)           (4)           "
      " (5)         (6)          (7)          (8)         (9)          "
      "   (10)           (11)         (12) \n");
  fprintf(fp,
          "# step  prev. step      time          prev. time        a         "
          " prev a          z          prev z    injection E    Numb SNIa    "
          "   SNIa rate     SNIa rate/V   Number\n");
  fflush(fp);
}

/**
 * @brief Initialize the SNIa global struct
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_init(const struct engine *e) {

  /* Initialize the core variables */
  event_logger_core_init(e, &log_SNIa.core);

  /* Make a constant for the physical constants */
  const struct phys_const *phys_const = e->physical_constants;

  /* Define the swift parameter file */
  struct swift_params *params = e->parameter_file;

  /* Initialize the detla time core value */
  const double delta_logger_time_Myr =
      parser_get_param_double(params, "Event_logger:delta_time_SNIa_Myr");

  /* Convert the time to internal units */
  log_SNIa.core.delta_logger_time =
      delta_logger_time_Myr * 1e6 * phys_const->const_year;

  /* Initialize the energy to zero */
  log_SNIa.SNIa_energy = 0.;

  /* Initialize the number of heating events to zero */
  log_SNIa.events = 0;
}

/**
 * @brief Do a time step and update the SNIa global struct
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_time_step(const struct engine *e) {

  event_logger_core_time_step(e, &log_SNIa.core);
}

/**
 * @brief Function that writes to the logger file
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_log_data_general(const struct engine *e,
                                                      const double dt) {

  /* Get the core struct */
  struct event_history_logger *core = &log_SNIa.core;

  /* Get the feedback structure */
  const struct feedback_props *feedback_properties = e->feedback_props;

  /* Calculate the volume of the box */
  const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];

  /* Calculate Delta time */
  const double delta_time = dt;

  /* Get the total amount of SNIa energy */
  const double E_SNIa = log_SNIa.SNIa_energy;

  /* Get the Energy of a single SNIa */
  const double E_single_SNIa = feedback_properties->E_SNIa;

  /* Calculate the number of SNIas in the simulation */
  const double N_SNIa = E_SNIa / E_single_SNIa;

  /* Calculate the number of SNIa per time and per time per volume */
  const double N_SNIa_p_time = N_SNIa / delta_time;
  const double N_SNIa_p_time_p_volume = N_SNIa_p_time / volume;

  /* Get the number of heating events */
  const int N_heating_events = log_SNIa.events;

  /* Set constants of time */
  const double a = e->cosmology->a;
  const double z = e->cosmology->z;
  const int step = e->step;
  const double time = e->time;

  /* Print the data to the file */
  fprintf(core->fp,
          "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f  %12.7e  %12.7e  "
          "%12.7e  %12.7e %7d \n",
          step, core->step_prev, time, core->time_prev, a, core->a_prev, z,
          core->z_prev, E_SNIa, N_SNIa, N_SNIa_p_time, N_SNIa_p_time_p_volume,
          N_heating_events);
  fflush(core->fp);
}

/**
 * @brief Write data to the feedback logger file if we are on a write step
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_log_data(const struct engine *e) {

  /* Get the core struct */
  struct event_history_logger *core = &log_SNIa.core;

  /* Are we one a logger time step? */
  if (!event_logger_core_log(e, core)) return;

  /* We need to log */
  event_logger_SNIa_log_data_general(e, log_SNIa.core.delta_logger_time);

  /* Update the logger core */
  event_logger_core_update(e, core);

  /* Update this logger specific values */
  log_SNIa.SNIa_energy = 0.;
  log_SNIa.events = 0;
}

/**
 * @brief Write data to the feedback logger file on the last time step
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_log_data_end(const struct engine *e) {

  /* We need to log before closing */
  event_logger_SNIa_log_data_general(e,
                                     log_SNIa.core.logger_time_since_last_log);

  /* Lets close the file */
  fclose(log_SNIa.core.fp);

#ifdef SWIFT_DEBUG_CHECKS
  /* Close the debugging file */
  fclose(log_SNIa_debug.fp);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief log a SNIa event
 *
 * @param si the spart of the feedback event pair
 * @param pj the part of the feedback event pair
 * @param xpj the xpart of the part of the feedback event pair
 * @param cosmo the cosmology struct
 */
INLINE static void event_logger_SNIa_log_event(const struct spart *si,
                                               const struct part *pj,
                                               const struct xpart *xpj,
                                               const struct cosmology *cosmo) {

  /* Get the injected energy */
  const double mass_init = pj->mass;
  const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
  const double deltaE = delta_u * mass_init;

  /* Write to the log */
  if (lock_lock(&log_SNIa.core.lock) == 0) {
    log_SNIa.SNIa_energy += deltaE;
    log_SNIa.events += 1;
  }
  if (lock_unlock(&log_SNIa.core.lock) != 0) error("Failed to unlock the lock");
}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication for the SNIa logger
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_MPI_Reduce(const struct engine *e) {

  /* Are we one a logger time step? */
  if (!event_logger_core_log(e, &log_SNIa.core)) return;

  /* Define empty variables for the MPI communication */
  int number_events_received;
  double total_SNIa_energy;

  MPI_Reduce(&log_SNIa.events, &number_events_received, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&log_SNIa.SNIa_energy, &total_SNIa_energy, 1, MPI_DOUBLE, MPI_SUM,
             0, MPI_COMM_WORLD);

  if (e->nodeID != 0) {
    /* Get the core struct */
    struct event_history_logger *core = &log_SNIa.core;

    /* Update the core struct */
    event_logger_core_update(e, core);

    /* Update the SNIa variables */
    log_SNIa.SNIa_energy = 0.;
    log_SNIa.events = 0;
    return;
  }

  /* Update the variables for node 0 */
  log_SNIa.SNIa_energy = total_SNIa_energy;
  log_SNIa.events = number_events_received;
}
#endif

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Initialize the SNIa logger debug file
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_init_log_file_debug(
    const struct engine *e) {

  /* Load the structures of the internal units and the physical constants */
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;

  /* Calculate the energy unit */
  const double E_unit = us->UnitMass_in_cgs * us->UnitLength_in_cgs *
                        us->UnitLength_in_cgs /
                        (us->UnitTime_in_cgs * us->UnitTime_in_cgs);

  /* Use the File pointer */
  FILE *fp = log_SNIa_debug.fp;

  /* Write some general text to the logger file */
  fprintf(fp, "# Stochastic SNIa Logger file\n");
  fprintf(fp, "######################################################\n");
  fprintf(fp, "# The quantities are all given in internal physical units!\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# (0) Simulation step\n");
  fprintf(fp,
          "# (1) Time since Big Bang (cosmological run), Time since start of "
          "the simulation (non-cosmological run).\n");
  fprintf(fp, "#     Unit = %e seconds\n", us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e yr or %e Myr\n", 1.f / phys_const->const_year,
          1.f / phys_const->const_year / 1e6);
  fprintf(fp, "# (2) Scale factor     (no unit)\n");
  fprintf(fp, "# (3) Redshift         (no unit)\n");
  fprintf(fp, "# (4) ID star particle (no unit)\n");
  fprintf(fp, "# (5) ID gas particle  (no unit)\n");
  fprintf(fp, "# (6) Injected energy of SNIa events\n");
  fprintf(fp, "#     Unit = %e erg\n", E_unit);
  fprintf(fp, "#     Unit = %e 10^51 erg\n", E_unit / 1e51);
  fprintf(fp, "# (7) Number of SNIa   (number, no unit)\n");
  fprintf(fp, "#\n");
  fprintf(
      fp,
      "# (0)         (1)            (2)          (3)            (4)           "
      " (5)            (6)            (7)\n");
  fprintf(fp,
          "#            Time             a            z        ID star part.  "
          "ID gas part.   Injected Energy  Number of SNIa\n");
}

/**
 * @brief Log the event in case of debugging
 *
 * @param time the current simulation time
 * @param si the star particle
 * @param pj the gas particle
 * @param xpj the extra information of the gas particle
 * @param cosmo the cosmology struct
 * @param step the current simulation step
 */
INLINE static void event_logger_SNIa_log_event_debug(
    const double time, const struct spart *si, struct part *pj,
    struct xpart *xpj, const struct cosmology *cosmo, const int step) {

  if (lock_lock(&log_SNIa_debug.lock) == 0) {

    /* Use the File pointer */
    FILE *fp = log_SNIa_debug.fp;

    /* Get the times */
    const double a = cosmo->a;
    const double z = cosmo->z;

    /* Get the injected energy */
    const double mass_init = si->mass_init;
    const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
    const double deltaE = delta_u * mass_init;

    fprintf(fp, "%6d %16e %12.7f %12.7f %14llu %14llu %16e %16e\n", step, time,
            a, z, si->id, pj->id, deltaE, deltaE * 1.9884e2);
    fflush(fp);
  }
  if (lock_unlock(&log_SNIa_debug.lock) != 0)
    error("Failed to unlock the lock");
}

/**
 * @brief Initialize the SNIa debugging global struct
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_SNIa_init_debug(const struct engine *e) {
  /* Initialize the lock*/
  lock_init(&log_SNIa_debug.lock);
}
#endif /* SWIFT_DEBUG_CHECKS */

#endif /* SWIFT_COLIBRE_EVENT_LOGGER_SNIA_H */
