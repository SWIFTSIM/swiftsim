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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_SNII_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_SNII_H

#include "feedback_logger_struct.h"
#include "feedback_logger_core.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Initialize the SNII logger file
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_init_log_file(
    const struct engine *restrict e) {

  /* Load the structures of the internal units and the physical constants */
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;

  /* Use the File pointer */
  FILE *fp = log_SNII.core.fp;

  /* Calculate the energy unit */
  const double E_unit = us->UnitMass_in_cgs * us->UnitLength_in_cgs *
                        us->UnitLength_in_cgs /
                        (us->UnitTime_in_cgs * us->UnitTime_in_cgs);

  /* Write some general text to the logger file */
  fprintf(fp, "# Stochastic SNII logger file\n");
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
          " prev a          z          prev z    injection E    Numb SNII    "
          "   SNII rate     SNII rate/V   Number\n");
  fflush(fp);
}

/**
 * @brief Initialize the SNII global struct
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_init(const struct engine *restrict e) {

  /* Initialize the core variables */
  feedback_logger_core_init(e, &log_SNII.core, 1);

  /* Initialize the energy to zero */
  log_SNII.SNII_energy = 0.;

  /* Initialize the number of SNII to zero */
  log_SNII.N_SNII = 0.;

  /* Initialize the number of heating events to zero */
  log_SNII.events = 0;
}

/**
 * @brief Do a time step and update the SNIa global struct
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_time_step(
    const struct engine *restrict e) {

  feedback_logger_core_time_step(e, &log_SNII.core);
}

/**
 * @brief Write data to the feedback logger file if we are on a write step
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_log_data_general(
    const struct engine *restrict e, const double dt) {

  /* Get the core struct */
  struct feedback_history_logger *core = &log_SNII.core;

  /* Get the feedback structure */
  const struct feedback_props *feedback_properties = e->feedback_props;

  /* Calculate the volume of the box */
  const double volume = e->s->dim[0] * e->s->dim[1] * e->s->dim[2];

  /* Calculate Delta time */
  const double delta_time = dt;

  /* Get the total amount of SNIa energy */
  const double E_SNII = log_SNII.SNII_energy;

  /* Get the Energy of a single SNIa */
  const double E_single_SNII = feedback_properties->E_SNII;

  /* Calculate the number of SNIas in the simulation */
  const double N_SNII = E_SNII / E_single_SNII;

  /* Calculate the number of SNIa per time and per time per volume */
  const double N_SNII_p_time = N_SNII / delta_time;
  const double N_SNII_p_time_p_volume = N_SNII_p_time / volume;

  /* Get the number of heating events */
  const int N_heating_events = log_SNII.events;

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
          core->z_prev, E_SNII, N_SNII, N_SNII_p_time, N_SNII_p_time_p_volume,
          N_heating_events);
  fflush(core->fp);
}


/**
 * @brief Write data to the feedback logger file if we are on a write step
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_log_data(
    const struct engine *restrict e) {

  /* Get the core struct */
  struct feedback_history_logger *core = &log_SNII.core;

  if (!feedback_logger_core_log(e, core)) return;

  /* We need to log */
  feedback_logger_SNII_log_data_general(e, log_SNII.core.delta_logger_time);

  /* Update the logger core */
  feedback_logger_core_update(e, core);

  /* Update the specific logger values of this logger */
  log_SNII.SNII_energy = 0.;
  log_SNII.events = 0;
  log_SNII.N_SNII = 0.;
}

/**
 * @brief Write data to the feedback logger file on the last time step
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_log_data_end(
    const struct engine *restrict e) {

  /* End of simulation so we need to log */
  feedback_logger_SNII_log_data_general(e, log_SNII.core.logger_time);  

  /* Close the logger file */
  fclose(log_SNII.core.fp);

#ifdef SWIFT_DEBUG_CHECKS
  /* Close the debugging file */
  fclose(log_SNII_debug.fp);
#endif /* SWIFT_DEBUG_CHECKS */
}

/**
 * @brief log a SNII event
 *
 * @param si the spart of the feedback event pair
 * @param pj the part of the feedback event pair
 * @param xpj the xpart of the part of the feedback event pair
 * @param cosmo the cosmology struct
 * @param f_E the energy fraction of the SNII event
 */
INLINE static void feedback_logger_SNII_log_event(
    const struct spart *restrict si, const struct part *restrict pj,
    const struct xpart *restrict xpj, const struct cosmology *restrict cosmo,
    const double f_E) {

  if (lock_lock(&log_SNII.core.lock) == 0) {

    /* Get the injected energy */
    const double mass_init = pj->mass;
    const double delta_u = si->feedback_data.to_distribute.SNII_delta_u;
    const double deltaE = delta_u * mass_init;

    /* Update the total SNII energy */
    log_SNII.SNII_energy += deltaE;
    log_SNII.events += 1;

    /* For the number of SNIIs first divide by the energy fraction, rest is done
     * while writing the data */
    log_SNII.N_SNII += deltaE / f_E;
  }
  if (lock_unlock(&log_SNII.core.lock) != 0) error("Failed to unlock the lock");
}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication for the SNII logger
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_MPI(const struct engine *restrict e) {

  /* Are we one a logger time step? */
  if (!feedback_logger_core_log(e, &log_SNII.core)) return;

  /* Define empty variables for the MPI communication */
  int number_events_received;
  double doubles_received[2];
  const double logger_doubles_send[2] = {log_SNII.SNII_energy, log_SNII.N_SNII};

  MPI_Reduce(&log_SNII.events, &number_events_received, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&logger_doubles_send, &doubles_received, 2, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);

  if (e->nodeID != 0) {
    /* Get the core struct */
    struct feedback_history_logger *core = &log_SNII.core;

    /* Update the core struct */
    feedback_logger_core_update(e, core);

    /* Update the SNIa variables */
    log_SNII.SNII_energy = 0.;
    log_SNII.N_SNII = 0.;
    log_SNII.events = 0;
    return;
  }

  /* Update the variables for node 0 */
  log_SNII.SNII_energy = doubles_received[0];
  log_SNII.N_SNII = doubles_received[1];
  log_SNII.events = number_events_received;
}
#endif

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Initialize the SNII logger debug file
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_init_log_file_debug(
    const struct engine *restrict e) {

  /* Load the structures of the internal units and the physical constants */
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;

  /* Calculate the energy unit */
  const double E_unit = us->UnitMass_in_cgs * us->UnitLength_in_cgs *
                        us->UnitLength_in_cgs /
                        (us->UnitTime_in_cgs * us->UnitTime_in_cgs);

  /* Use the File pointer */
  FILE *fp = log_SNII_debug.fp;

  /* Write some general text to the logger file */
  fprintf(fp, "# Stochastic SNII Debugging Logger file\n");
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
  fprintf(fp, "# (7) Age of the star particle (Myr)\n");
  fprintf(fp, "#\n");
  fprintf(
      fp,
      "# (0)         (1)            (2)          (3)            (4)           "
      " (5)            (6)            (7)\n");
  fprintf(fp,
          "#            Time             a            z        ID star part.  "
          "ID gas part.   Injected Energy  Age of star \n");
  fflush(fp);
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
INLINE static void feedback_logger_SNII_log_event_debug(
    const double time, const struct spart *restrict si,
    struct part *restrict pj, struct xpart *restrict xpj,
    const struct cosmology *restrict cosmo, const int step) {

  if (lock_lock(&log_SNII_debug.lock) == 0) {

    /* Use the File pointer */
    FILE *fp = log_SNII_debug.fp;

    /* Get the times */
    const double a = cosmo->a;
    const double z = cosmo->z;

    /* Get the injected energy */
    const double mass_init = si->mass_init;
    const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
    const double deltaE = delta_u * mass_init;

    /* The age of the star particle */
    const float age_star = si->feedback_data.to_distribute.SNII_star_age_Myr;

    fprintf(fp, "%6d %16e %12.7f %12.7f %14llu %14llu %16e %.4f\n", step, time,
            a, z, si->id, pj->id, deltaE, age_star);
    fflush(fp);
  }
  if (lock_unlock(&log_SNII_debug.lock) != 0)
    error("Failed to unlock the lock");
}

/**
 * @brief Initialize the SNII debugging global struct
 *
 * @param e the engine we are running
 */
INLINE static void feedback_logger_SNII_init_debug(
    const struct engine *restrict e) {
  /* Initialize the lock*/
  lock_init(&log_SNII_debug.lock);
}
#endif /* SWIFT_DEBUG_CHECKS */

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_SNII_H */
