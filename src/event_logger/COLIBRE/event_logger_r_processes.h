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
#ifndef SWIFT_COLIBRE_EVENT_LOGGER_R_PROCESSES_H
#define SWIFT_COLIBRE_EVENT_LOGGER_R_PROCESSES_H

/* Local includes */
#include "event_logger_core.h"
#include "event_logger_struct.h"
#include "feedback_properties.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* feedback history struct for r-processes */
struct feedback_history_r_processes {

  /* Load the core of logging functions */
  struct event_history_logger core;

  /*! Total new r-processes mass in the simulation */
  double NSM_enrichment_mass;
  double CEJSN_enrichment_mass;
  double collapsar_enrichment_mass;

  /*! Number of r-processes events */
  int NSM_events;
  int CEJSN_events;
  int collapsar_events;
};

/**
 * @brief Initialize the r-processes logger file
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_init_log_file(
    const struct engine *e) {

  /* Load the structures of the internal units and the physical constants */
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;

  /* Use the File pointer */
  FILE *fp = log_r_processes.core.fp;

  /* Write some general text to the logger file */
  fprintf(fp, "# Stochastic r-processes logger file\n");
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
  fprintf(fp, "# (8)  Injected mass of r-processes\n");
  fprintf(fp, "#      Unit = %e g\n", us->UnitMass_in_cgs);
  fprintf(fp, "#      Unit = %e solar mass\n",
          1. / phys_const->const_solar_mass);
  fprintf(fp, "# (9)  Injected mass of r-processes per time \n");
  fprintf(fp, "#      Unit = %e g/s\n",
          us->UnitMass_in_cgs / us->UnitTime_in_cgs);
  fprintf(fp, "#      Unit = %e solar mass/yr\n",
          1. / phys_const->const_solar_mass * phys_const->const_year);
  fprintf(fp, "# (10) Injected mass of r-processes per time per volume\n");
  fprintf(fp, "#      Unit = %e g/s/cm^3\n",
          us->UnitMass_in_cgs / us->UnitTime_in_cgs /
              pow(us->UnitLength_in_cgs, 3));
  fprintf(fp, "#      Unit = %e solar mass/yr/Mpc^3\n",
          1. / phys_const->const_solar_mass * phys_const->const_year *
              pow(phys_const->const_parsec * 1e6, 3));
  fprintf(fp,
          "# (11) Number of r-processes (binned in time between current time "
          "and previous time)\n");
  fprintf(fp, "#      Unit = no unit\n");
  fprintf(fp,
          "# (12) Number of r-processes per time (binned in time between "
          "current time "
          "and previous time).\n");
  fprintf(fp, "#      Unit = %e #/seconds\n", 1. / us->UnitTime_in_cgs);
  fprintf(fp, "#      Unit = %e #/yr or %e #/Myr\n", phys_const->const_year,
          phys_const->const_year * 1e6);
  fprintf(fp,
          "# (13) Number of r-processes per time per comoving volume (binned "
          "in time "
          "between current time and previous time).\n");
  fprintf(fp, "#      Unit = %e #/seconds/cm^3\n",
          1. / us->UnitTime_in_cgs / pow(us->UnitLength_in_cgs, 3));
  fprintf(fp, "#      Unit = %e #/yr/Mpc^3\n",
          phys_const->const_year * pow(phys_const->const_parsec * 1e6, 3));
  fprintf(fp,
          "# (14) Number of NSM events (binned in time between current time "
          "and previous time)\n");
  fprintf(fp, "#      Unit = no unit\n");
  fprintf(fp,
          "# (15) Number of common-envelop jets SN events (binned in time "
          "between current time "
          "and previous time)\n");
  fprintf(fp, "#      Unit = no unit\n");
  fprintf(
      fp,
      "# (16) Number of collapsar events (binned in time between current time "
      "and previous time)\n");
  fprintf(fp, "#      Unit = no unit\n");

  fprintf(fp, "# (17) Injected mass by NSM events\n");
  fprintf(fp, "#      Unit = %e g\n", us->UnitMass_in_cgs);
  fprintf(fp, "#      Unit = %e solar mass\n",
          1. / phys_const->const_solar_mass);
  fprintf(fp, "# (18) Injected mass by CEJSN events\n");
  fprintf(fp, "#      Unit = %e g\n", us->UnitMass_in_cgs);
  fprintf(fp, "#      Unit = %e solar mass\n",
          1. / phys_const->const_solar_mass);
  fprintf(fp, "# (19) Injected mass by collapsar events\n");
  fprintf(fp, "#      Unit = %e g\n", us->UnitMass_in_cgs);
  fprintf(fp, "#      Unit = %e solar mass\n",
          1. / phys_const->const_solar_mass);
  fprintf(fp, "#\n");
  fprintf(
      fp,
      "#  (0)      (1)         (2)              (3)           (4)           "
      " (5)         (6)          (7)          (8)         (9)          "
      "   (10)           (11)        (12)         (13)            (14)   "
      " (15)    (16)    (17)          (18)          (19)\n");
  fprintf(fp,
          "# step  prev. step      time          prev. time        a         "
          " prev a          z          prev z     Inj. mass     Inj. mass rate"
          "  Inj.mass rate/V   N        N rate       N rate/V           N     "
          "  N       N   Inj. mass       Inj. mass       Inj. mass\n");
  fflush(fp);
}

/**
 * @brief Initialize the r-processes global struct
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_init(const struct engine *e) {

  /* Initialize the core variables */
  event_logger_core_init(e, &log_r_processes.core);

  /* Make a constant for the physical constants */
  const struct phys_const *phys_const = e->physical_constants;

  /* Define the swift parameter file */
  struct swift_params *params = e->parameter_file;

  /* Initialize the detla time core value */
  const double delta_logger_time_Myr = parser_get_param_double(
      params, "Event_logger:delta_time_r_processes_Myr");

  /* Convert the time to internal units */
  log_r_processes.core.delta_logger_time =
      delta_logger_time_Myr * 1e6 * phys_const->const_year;

  /* Initialize the energy to zero */
  log_r_processes.NSM_enrichment_mass = 0.;
  log_r_processes.CEJSN_enrichment_mass = 0.;
  log_r_processes.collapsar_enrichment_mass = 0.;

  /* Initialize the number of heating events to zero */
  log_r_processes.NSM_events = 0;
  log_r_processes.CEJSN_events = 0;
  log_r_processes.collapsar_events = 0;
}

/**
 * @brief Do a time step and update the r-processes global struct
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_time_step(const struct engine *e) {

  event_logger_core_time_step(e, &log_r_processes.core);
}

/**
 * @brief Function that writes to the logger file
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_log_data_general(
    const struct engine *e, const double dt) {

  /* Get the core struct */
  struct event_history_logger *core = &log_r_processes.core;

  /* Calculate the volume of the box */
  const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];

  /* Calculate the inverse of the volume of the box */
  const double volume_inv = 1. / volume;

  /* Calculate Delta time */
  const double delta_time = log_r_processes.core.delta_logger_time;

  /* Calculate the inverse of the delta time */
  const double delta_time_inv = 1. / delta_time;

  /* Get the total number of r-process events */
  const int N_r_processes = log_r_processes.NSM_events +
                            log_r_processes.CEJSN_events +
                            log_r_processes.collapsar_events;

  /* Get the number of NSM, CEJSN and collapsar events separately */
  const int N_NSM_events = log_r_processes.NSM_events;
  const int N_CEJSN_events = log_r_processes.CEJSN_events;
  const int N_collapsar_events = log_r_processes.collapsar_events;

  /* Get the total number of r-processes per time */
  const double N_r_processes_p_time = (double)N_r_processes * delta_time_inv;

  /* Get the total number of r-processes per time per volume */
  const double N_r_processes_p_time_p_volume =
      N_r_processes_p_time * volume_inv;

  /* Get the total enrichment mass */
  const double delta_mass = log_r_processes.NSM_enrichment_mass +
                            log_r_processes.CEJSN_enrichment_mass +
                            log_r_processes.collapsar_enrichment_mass;

  /* Get enrichment mass for each channel */
  const double NSM_delta_mass = log_r_processes.NSM_enrichment_mass;
  const double CEJSN_delta_mass = log_r_processes.CEJSN_enrichment_mass;
  const double collapsar_delta_mass = log_r_processes.collapsar_enrichment_mass;

  /* Get the enrichment mass per time */
  const double delta_mass_p_time = delta_mass * delta_time_inv;

  /* Get the enrichment mass per time per volume */
  const double delta_mass_p_time_p_volume = delta_mass_p_time * volume_inv;

  /* Set constants of time */
  const double a = e->cosmology->a;
  const double z = e->cosmology->z;
  const int step = e->step;
  const double time = e->time;

  /* Print the data to the file */
  fprintf(core->fp,
          "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f  %12.7e  %12.7e  "
          "%12.7e %7d      %12.7e %12.7e %7d %7d %7d  %12.7e  %12.7e  %12.7e\n",
          step, core->step_prev, time, core->time_prev, a, core->a_prev, z,
          core->z_prev, delta_mass, delta_mass_p_time,
          delta_mass_p_time_p_volume, N_r_processes, N_r_processes_p_time,
          N_r_processes_p_time_p_volume, N_NSM_events, N_CEJSN_events,
          N_collapsar_events, NSM_delta_mass, CEJSN_delta_mass,
          collapsar_delta_mass);
  fflush(core->fp);
}

/**
 * @brief Write data to the feedback logger file if we are on a write step
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_log_data(const struct engine *e) {

  /* Get the core struct */
  struct event_history_logger *core = &log_r_processes.core;

  /* Are we one a logger time step? */
  if (!event_logger_core_log(e, core)) return;

  /* We need to log */
  event_logger_r_processes_log_data_general(
      e, log_r_processes.core.delta_logger_time);

  /* Update the logger core */
  event_logger_core_update(e, core);

  /* Update the type specific variables */
  log_r_processes.NSM_events = 0;
  log_r_processes.CEJSN_events = 0;
  log_r_processes.collapsar_events = 0;
  log_r_processes.NSM_enrichment_mass = 0.;
  log_r_processes.CEJSN_enrichment_mass = 0.;
  log_r_processes.collapsar_enrichment_mass = 0.;
}

/**
 * @brief Write data to the feedback logger file on the last time step
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_log_data_end(
    const struct engine *e) {

  /* Write on the last time step */
  event_logger_r_processes_log_data_general(
      e, log_r_processes.core.logger_time_since_last_log);

  /* Close the logger file */
  fclose(log_r_processes.core.fp);
}

/**
 * @brief log a r-process event
 *
 * @param si the spart of the enrichment event
 * @param cosmo the cosmology struct
 * @param delta_mass the enrichment mass of the r-processes
 * @param num_events number of events per time step
 * @param flag, indication of which event: neutron stars (flag==0), rare SN
 * (flag==1), collapsars (flag==2)
 */
INLINE static void event_logger_r_processes_log_event(
    const struct spart *si, const struct cosmology *cosmo,
    const double delta_mass, const int num_events, const int flag) {

  if (lock_lock(&log_r_processes.core.lock) == 0) {

    /* Store the injected mass and add an event */

    if (flag == 0) {
      log_r_processes.NSM_events += num_events;
      log_r_processes.NSM_enrichment_mass += delta_mass;
    }
    if (flag == 1) {
      log_r_processes.CEJSN_events += num_events;
      log_r_processes.CEJSN_enrichment_mass += delta_mass;
    }
    if (flag == 2) {
      log_r_processes.collapsar_events += num_events;
      log_r_processes.collapsar_enrichment_mass += delta_mass;
    }
  }
  if (lock_unlock(&log_r_processes.core.lock) != 0)
    error("Failed to unlock the lock");
}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication for the r-processes logger
 *
 * @param e the engine we are running
 */
INLINE static void event_logger_r_processes_MPI_Reduce(const struct engine *e) {

  /* Are we one a logger time step? */
  if (!event_logger_core_log(e, &log_r_processes.core)) return;

  /* Define empty variables for the MPI communication */
  int number_events_received_NSM;
  int number_events_received_CEJSN;
  int number_events_received_coll;
  double total_mass_NSM;
  double total_mass_CEJSN;
  double total_mass_coll;

  MPI_Reduce(&log_r_processes.NSM_events, &number_events_received_NSM, 1, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_r_processes.CEJSN_events, &number_events_received_CEJSN, 1, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_r_processes.collapsar_events, &number_events_received_coll, 1,
             MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_r_processes.NSM_enrichment_mass, &total_mass_NSM, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_r_processes.CEJSN_enrichment_mass, &total_mass_CEJSN, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_r_processes.collapsar_enrichment_mass, &total_mass_coll, 1,
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (e->nodeID != 0) {
    /* Get the core struct */
    struct event_history_logger *core = &log_r_processes.core;

    /* Update the core struct */
    event_logger_core_update(e, core);

    /* Update the r-process variables */
    log_r_processes.NSM_enrichment_mass = 0.;
    log_r_processes.CEJSN_enrichment_mass = 0.;
    log_r_processes.collapsar_enrichment_mass = 0.;
    log_r_processes.NSM_events = 0;
    log_r_processes.CEJSN_events = 0;
    log_r_processes.collapsar_events = 0;
    return;
  }

  /* Update the variables for node 0 */
  log_r_processes.NSM_enrichment_mass = total_mass;
  log_r_processes.CEJSN_enrichment_mass = total_mass;
  log_r_processes.collapsar_enrichment_mass = total_mass;
  log_r_processes.NSM_events = number_events_received;
  log_r_processes.CEJSN_events = number_events_received;
  log_r_processes.collapsar_events = number_events_received;
}
#endif

#endif /* SWIFT_COLIBRE_EVENT_LOGGER_R_PROCESSES_H */
