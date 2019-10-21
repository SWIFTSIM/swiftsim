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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_SNIA_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_SNIA_H

#include "feedback_logger_core.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Initialize the SNIa logger file
 *
 * @param e the engine we are running 
 */
INLINE static void feedback_logger_SNIa_init_log_file(const struct engine *restrict e) {

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
}

INLINE static void feedback_logger_SNIa_init(const struct engine *restrict e) {
 
  /* Initialize the core variables */ 
  feedback_logger_core_init(e, &log_SNIa.core,2);

  /* Initialize the energy to zero */
  log_SNIa.SNIa_energy = 0.;

  /* Initialize the number of heating events to zero */
  log_SNIa.events = 0;
 
}

INLINE static void feedback_logger_SNIa_time_step(const struct engine *restrict e) {
  
  feedback_logger_core_time_step(e, &log_SNIa.core);
}

INLINE static void feedback_logger_SNIa_log_data(const struct engine *restrict e) {
  
  /* Are we one a logger time step? */
  if (!feedback_logger_core_log(e, &log_SNIa.core)) return;

  /* We need to log */
  
  /* Get the feedback structure */
  const struct feedback_props *feedback_properties = e->feedback_props; 

  /* Get the core struct */
  struct feedback_history_logger *core = &log_SNIa.core;

  /* Calculate the volume of the box */
  const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];

  /* Calculate Delta time */
  const double delta_time = log_SNIa.core.delta_logger_time;

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
  
  feedback_logger_core_update(e,core);

}

#ifdef WITH_MPI
INLINE static void feedback_logger_SNIa_MPI(const struct engine *restrict e) {

  /* Are we one a logger time step? */
  if (!feedback_logger_core_log(e, &log_SNIa.core)) return;

  /* Define empty variables for the MPI communication */
  int number_events_received;
  double total_SNIa_energy;

  MPI_Reduce(&log_SNIa.events, &number_events_received, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&log_SNIa.SNIa_energy, &total_SNIa_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (e->nodeID != 0) {
    /* Get the core struct */
    struct feedback_history_logger *core = &log_SNIa.core;   

    /* Update the core struct */
    feedback_logger_core_update(e,core);

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

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_SNIA_H */
