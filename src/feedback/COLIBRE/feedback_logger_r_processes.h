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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_R_PROCESSES_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_R_PROCESSES_H

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Initialize the r-processes logger file
 *
 * @param e the engine we are running 
 */
INLINE static void feedback_logger_r_processes_init_log_file(const struct engine *restrict e) {

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
  fprintf(fp, "#\n");
  fprintf(
      fp,
      "#  (0)      (1)         (2)              (3)           (4)           "
      " (5)         (6)          (7)          (8)         (9)          "
      "   (10)           (11)        (12)         (13)\n");
  fprintf(fp,
          "# step  prev. step      time          prev. time        a         "
          " prev a          z          prev z     Inj. mass     Inj. mass rate"
          "  Inj.mass rate/V   N        N rate       N rate/V \n");
  fflush(fp);
}

INLINE static void feedback_logger_r_processes_init(const struct engine *restrict e) {
 
  /* Initialize the core variables */ 
  feedback_logger_core_init(e, &log_r_processes.core,3);

  /* Initialize the energy to zero */
  log_r_processes.enrichment_mass = 0.;

  /* Initialize the number of heating events to zero */
  log_r_processes.events = 0;
 
}

INLINE static void feedback_logger_r_processes_time_step(const struct engine *restrict e) {
  
  feedback_logger_core_time_step(e, &log_r_processes.core);
}

INLINE static void feedback_logger_r_processes_log_data(const struct engine *restrict e) {

  /* Are we one a logger time step? */
  if (!feedback_logger_core_log(e, &log_r_processes.core)) return;

  /* Get the core struct */
  struct feedback_history_logger *core = &log_r_processes.core;

  /* Calculate the volume of the box */
  const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];

  /* Calculate Delta time */
  const double delta_time = log_r_processes.core.delta_logger_time;

  const int N_r_processes = log_r_processes.events;

  const double N_r_processes_p_time = (double)N_r_processes / delta_time;

  const double N_r_processes_p_time_p_volume = N_r_processes_p_time / volume;

  const double delta_mass = log_r_processes.enrichment_mass;

  const double delta_mass_p_time = delta_mass / delta_time;

  const double delta_mass_p_time_p_volume = delta_mass_p_time / volume;

  /* Set constants of time */
  const double a = e->cosmology->a;
  const double z = e->cosmology->z;
  const int step = e->step; 
  const double time = e->time;

  /* Print the data to the file */
  fprintf(core->fp,
          "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f  %12.7e  %12.7e  "
          "%12.7e %7d      %12.7e %12.7e \n",
          step, core->step_prev, time, core->time_prev, a, core->a_prev, z,
          core->z_prev, delta_mass, delta_mass_p_time,
          delta_mass_p_time_p_volume, N_r_processes, N_r_processes_p_time,
          N_r_processes_p_time_p_volume);
  fflush(core->fp);

  feedback_logger_core_update(e,core);

  log_r_processes.events = 0;
  log_r_processes.enrichment_mass = 0.;
}

INLINE static void feedback_logger_r_processes_log_event(
    const struct spart *restrict si, const struct part *restrict pj,
    const struct xpart *restrict xpj, const struct cosmology *restrict cosmo, const double delta_mass) {
  
  if (lock_lock(&log_r_processes.core.lock) == 0) {
    /* Get the injected energy */

    log_r_processes.enrichment_mass += delta_mass;
    log_r_processes.events += 1;
  }
  if (lock_unlock(&log_r_processes.core.lock) != 0) error("Failed to unlock the lock");
}

#ifdef WITH_MPI
INLINE static void feedback_logger_r_processes_MPI(const struct engine *restrict e) {

  /* Are we one a logger time step? */
  if (!feedback_logger_core_log(e, &log_r_processes.core)) return;

  /* Define empty variables for the MPI communication */
  int number_events_received;
  double total_mass;

  MPI_Reduce(&log_r_processes.events, &number_events_received, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_r_processes.enrichment_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (e->nodeID != 0) {
    /* Get the core struct */
    struct feedback_history_logger *core = &log_r_processes.core;   

    /* Update the core struct */
    feedback_logger_core_update(e,core);

    /* Update the SNIa variables */
    log_r_processes.enrichment_mass = 0.;
    log_r_processes.events = 0;
    return;
  }

  /* Update the variables for node 0 */
  log_r_processes.enrichment_mass = total_mass;
  log_r_processes.events = number_events_received;

}
#endif

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_R_PROCESSES_H */
