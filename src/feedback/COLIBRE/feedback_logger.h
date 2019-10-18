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

#include "feedback_logger_struct.h"
#include "feedback_properties.h"
#include "engine.h"

extern swift_lock_type lock_SNIa;
extern swift_lock_type lock_SNII;
extern swift_lock_type lock_r_processes;

/**
 * @brief Initialize the SNIa logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_logger_SNIa_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

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

/**
 * @brief Initialize the SNIa logger debug file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_logger_SNIa_init_log_file_debug(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

  /* Calculate the energy unit */
  const double E_unit = us->UnitMass_in_cgs * us->UnitLength_in_cgs *
                        us->UnitLength_in_cgs /
                        (us->UnitTime_in_cgs * us->UnitTime_in_cgs);

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
 * @brief Initialize the times of the feedback logger
 *
 * @param fha the feedback history accumulator struct
 * @param time the current simulation time
 * @param a the current scale factor
 * @param z the current redshift
 */
INLINE static void feedback_logger_init(
    struct feedback_history_accumulator *fha, const double time, const double a,
    const double z) {

  /* Initialize the start times of the simulation */
  fha->step_prev = 0;
  fha->time_prev = time;
  fha->z_prev = z;
  fha->a_prev = a;

  /* Initialize the feedback history logger locks */
  lock_init(&lock_SNIa);
  lock_init(&lock_SNII);
  lock_init(&lock_r_processes);
}

/**
 * @brief log all the collected data to the logger file
 *
 * @param feedback_properties, the feedback structure
 * @param fp the file pointer to write to.
 * @param SNIa the external variable that stores all the feedback information
 * @param fha the feedback history accumulator struct
 * @param step the current simulation step
 * @param time the current simulation time
 * @param a the current scale factor
 * @param z the current redshift
 * @param volume the simulation volume
 */
INLINE static void feedback_logger_SNIa_log_data(
    const struct feedback_props *restrict feedback_properties, FILE *fp,
    struct feedback_history_SNIa *restrict SNIa,
    const struct feedback_history_accumulator *fha, const int step,
    const double time, const double a, const double z, const double volume) {

  /* Calculate Delta time */
  const double delta_time = time - fha->time_prev;

  /* Get the total amount of SNIa energy */
  const double E_SNIa = SNIa->SNIa_energy;

  /* Get the Energy of a single SNIa */
  const double E_single_SNIa = feedback_properties->E_SNIa;

  /* Calculate the number of SNIas in the simulation */
  const double N_SNIa = E_SNIa / E_single_SNIa;

  /* Calculate the number of SNIa per time and per time per volume */
  const double N_SNIa_p_time = N_SNIa / delta_time;
  const double N_SNIa_p_time_p_volume = N_SNIa_p_time / volume;

  /* Get the number of heating events */
  const int N_heating_events = SNIa->heating;

  /* Print the data to the file */
  fprintf(fp,
          "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f  %12.7e  %12.7e  "
          "%12.7e  %12.7e %7d \n",
          step, fha->step_prev, time, fha->time_prev, a, fha->a_prev, z,
          fha->z_prev, E_SNIa, N_SNIa, N_SNIa_p_time, N_SNIa_p_time_p_volume,
          N_heating_events);
}

INLINE static void feedback_logger_update_times(
    struct feedback_history_accumulator *fha, const int step, const double time,
    const double a, const double z) {

  /* Set the times to the new values */
  fha->step_prev = step;
  fha->time_prev = time;
  fha->a_prev = a;
  fha->z_prev = z;
}

/**
 * @brief log a SNIa event
 *
 * @param SNIa the external variable that stores all the feedback information
 * @param time the current simulation time
 * @param si the star particle
 * @param pj the gas particle
 * @param xpj the extra information of the gas particle
 * @param cosmo the cosmology struct
 * @param step the current simulation step
 */
INLINE static void feedback_logger_SNIa_log_event(
    struct feedback_history_SNIa *restrict SNIa, const double time,
    const struct spart *restrict si, const struct part *restrict pj,
    const struct xpart *restrict xpj, const struct cosmology *restrict cosmo,
    const int step) {

  if (lock_lock(&lock_SNIa) == 0) {

    /* Get the injected energy */
    const double mass_init =
        pj->mass;  // I first had this but this doesn't seem right: const double
                   // mass_init = si->mass_init;
    const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
    const double deltaE = delta_u * mass_init;

    SNIa->SNIa_energy += deltaE;
    SNIa->heating += 1;
    message("Event!!!! BAM!! E=%e", SNIa->SNIa_energy);
  }
  if (lock_unlock(&lock_SNIa) != 0) error("Failed to unlock the lock");
}

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
    const struct cosmology *restrict cosmo, const int step) {

  if (lock_lock(&lock_SNIa) == 0) {

    /* Get the times */
    const double a = cosmo->a;
    const double z = cosmo->z;

    /* Get the injected energy */
    const double mass_init = si->mass_init;
    const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
    const double deltaE = delta_u * mass_init;

    fprintf(fp, "%6d %16e %12.7f %12.7f %14llu %14llu %16e %16e\n", step, time,
            a, z, si->id, pj->id, deltaE, deltaE * 1.9884e2);
  }
  if (lock_unlock(&lock_SNIa) != 0) error("Failed to unlock the lock");
}

/**
 * @brief function to reinitialize the SNIa logger.
 *
 * @param SNIa the external variable that stores all the feedback information
 */
INLINE static void feedback_logger_SNIa_clear(
    struct feedback_history_SNIa *restrict SNIa) {

  /* Set all the variables to zero */
  if (lock_lock(&lock_SNIa) == 0) {
    SNIa->SNIa_energy = 0;
    SNIa->heating = 0;
  }
  if (lock_unlock(&lock_SNIa) != 0) error("Failed to unlock the lock");
}

/**
 * @brief Initialize the SNII logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_logger_SNII_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

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
}

/**
 * @brief log all the collected data to the SNII logger file
 *
 * @param feedback_properties, the feedback structure
 * @param fp the file pointer to write to.
 * @param SNII the external variable that stores all the feedback information
 * @param fha the feedback history accumulator struct
 * @param step the current simulation step
 * @param time the current simulation time
 * @param a the current scale factor
 * @param z the current redshift
 * @param volume the simulation volume
 */
INLINE static void feedback_logger_SNII_log_data(
    const struct feedback_props *restrict feedback_properties, FILE *fp,
    struct feedback_history_SNII *restrict SNII,
    const struct feedback_history_accumulator *fha, const int step,
    const double time, const double a, const double z, const double volume) {

  /* Calculate Delta time */
  const double delta_time = time - fha->time_prev;

  /* Get the total amount of SNII energy */
  const double E_SNII = SNII->SNII_energy;

  /* Get the Energy of a single SNII */
  const double E_single_SNII = feedback_properties->E_SNII;

  /* Calculate the number of SNIIs in the simulation */
  const double N_SNII = SNII->N_SNII / E_single_SNII;

  /* Calculate the number of SNIa per time and per time per volume */
  const double N_SNII_p_time = N_SNII / delta_time;
  const double N_SNII_p_time_p_volume = N_SNII_p_time / volume;

  /* Get the number of heating events */
  const int N_heating_events = SNII->heating;

  /* Print the data to the file */
  fprintf(fp,
          "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f  %12.7e  %12.7e  "
          "%12.7e  %12.7e %7d \n",
          step, fha->step_prev, time, fha->time_prev, a, fha->a_prev, z,
          fha->z_prev, E_SNII, N_SNII, N_SNII_p_time, N_SNII_p_time_p_volume,
          N_heating_events);
}

/**
 * @brief function to reinitialize the SNII logger.
 *
 * @param SNII the external variable that stores all the feedback information
 */
INLINE static void feedback_logger_SNII_clear(
    struct feedback_history_SNII *restrict SNII) {

  /* Set all the variables to zero */
  if (lock_lock(&lock_SNII) == 0) {
    SNII->SNII_energy = 0.f;
    SNII->heating = 0;
    SNII->N_SNII = 0.f;
  }
  if (lock_unlock(&lock_SNII) != 0) error("Failed to unlock the lock");
}

/**
 * @brief log a SNII event
 *
 * @param SNII the external variable that stores all the feedback information
 * @param time the current simulation time
 * @param si the star particle
 * @param pj the gas particle
 * @param xpj the extra information of the gas particle
 * @param cosmo the cosmology struct
 * @param step the current simulation step
 */
INLINE static void feedback_logger_SNII_log_event(
    struct feedback_history_SNII *restrict SNII, const double time,
    const struct spart *restrict si, const struct part *restrict pj,
    const struct xpart *restrict xpj, const struct cosmology *restrict cosmo,
    const int step, const float f_E) {

  if (lock_lock(&lock_SNII) == 0) {

    /* Get the injected energy */
    const double mass_init =
        pj->mass;  // I first had this but this doesn't seem right: const double
                   // mass_init = si->mass_init;
    const double delta_u = si->feedback_data.to_distribute.SNII_delta_u;
    const double deltaE = delta_u * mass_init;

    /* Update the total SNII energy */
    SNII->SNII_energy += deltaE;
    SNII->heating += 1;

    /* For the number of SNIIs first divide by the energy fraction, rest is done
     * while writing the data */
    SNII->N_SNII += deltaE / f_E;
  }
  if (lock_unlock(&lock_SNII) != 0) error("Failed to unlock the lock");
}

/**
 * @brief Initialize the r-processes logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_logger_r_processes_init_log_file(
    FILE *fp, const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

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
}

/**
 * @brief log all the collected data to the r-processes logger file
 *
 * @param feedback_properties, the feedback structure
 * @param fp the file pointer to write to.
 * @param r_process the external variable that stores all the feedback
 * information
 * @param fha the feedback history accumulator struct
 * @param step the current simulation step
 * @param time the current simulation time
 * @param a the current scale factor
 * @param z the current redshift
 * @param volume the simulation volume
 */
INLINE static void feedback_logger_r_processes_log_data(
    const struct feedback_props *restrict feedback_properties, FILE *fp,
    struct feedback_history_r_processes *restrict r_process,
    struct feedback_history_accumulator *fha, const int step, const double time,
    const double a, const double z, const double volume) {

  /* Calculate Delta time */
  const double delta_time = time - fha->time_prev;

  const int N_r_processes = r_process->events;

  const double N_r_processes_p_time = (double)N_r_processes / delta_time;

  const double N_r_processes_p_time_p_volume = N_r_processes_p_time / volume;

  const double delta_mass = r_process->enrichment_mass;

  const double delta_mass_p_time = delta_mass / delta_time;

  const double delta_mass_p_time_p_volume = delta_mass_p_time / volume;

  /* Print the data to the file */
  fprintf(fp,
          "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f  %12.7e  %12.7e  "
          "%12.7e %7d      %12.7e %12.7e \n",
          step, fha->step_prev, time, fha->time_prev, a, fha->a_prev, z,
          fha->z_prev, delta_mass, delta_mass_p_time,
          delta_mass_p_time_p_volume, N_r_processes, N_r_processes_p_time,
          N_r_processes_p_time_p_volume);
}

/**
 * @brief function to reinitialize the r_processes logger.
 *
 * @param SNII the external variable that stores all the feedback information
 */
INLINE static void feedback_logger_r_processes_clear(
    struct feedback_history_r_processes *restrict r_process) {

  /* Set the enrichment mass to zero */
  r_process->enrichment_mass = 0.;

  /* Set the number of events to zero */
  r_process->events = 0;
}

/**
 * @brief log a r-process event
 *
 * @param r_processes the external variable that stores all the feedback
 * information
 * @param time the current simulation time
 * @param si the star particle
 * @param pj the gas particle
 * @param xpj the extra information of the gas particle
 * @param cosmo the cosmology struct
 * @param delta_mass the new r-processes mass in Europium
 */
INLINE static void feedback_logger_r_processes_log_event(
    struct feedback_history_r_processes *restrict r_processes,
    const double time, const struct spart *restrict si,
    const struct part *restrict pj, const struct xpart *restrict xpj,
    const struct cosmology *restrict cosmo, const double delta_mass) {

  if (lock_lock(&lock_r_processes) == 0) {
    r_processes->enrichment_mass += delta_mass;
    r_processes->events += 1;
  }
  if (lock_unlock(&lock_r_processes) != 0) error("Failed to unlock the lock");
}

/**
 * @brief log all the collected data to the logger file
 *
 * @param feedback_properties, the feedback structure
 * @param fp the file pointer to write to.
 * @param SNIa the external variable that stores all the feedback information
 * @param fha the feedback history accumulator struct
 * @param step the current simulation step
 * @param time the current simulation time
 * @param a the current scale factor
 * @param z the current redshift
 * @param volume the simulation volume
 */
INLINE static void feedback_logger_log_data(
    const struct feedback_props *restrict feedback_properties, FILE *fp_SNII,
    struct feedback_history_SNII *restrict SNII, FILE *fp_SNIa,
    struct feedback_history_SNIa *restrict SNIa, FILE *fp_r_processes,
    struct feedback_history_r_processes *restrict r_processes,
    struct feedback_history_accumulator *fha, const int step, const double time,
    const double a, const double z, const double volume) {

  /* Log the SNIa data */
  feedback_logger_SNIa_log_data(feedback_properties, fp_SNIa, SNIa, fha, step,
                                time, a, z, volume);
  feedback_logger_SNIa_clear(SNIa);
  fflush(fp_SNIa);

  /* Log the SNII data */
  feedback_logger_SNII_log_data(feedback_properties, fp_SNII, SNII, fha, step,
                                time, a, z, volume);
  feedback_logger_SNII_clear(SNII);
  fflush(fp_SNII);

  /* Log the r_processes data */
  feedback_logger_r_processes_log_data(feedback_properties, fp_r_processes,
                                       r_processes, fha, step, time, a, z,
                                       volume);
  feedback_logger_r_processes_clear(r_processes);
  fflush(fp_r_processes);

  /* Update the times in the gneral logger struct */
  feedback_logger_update_times(fha, step, time, a, z);
}

/**
 * @brief Initialize all the feedback log file 
 *
 * @param fp_SNII the SNII file pointer
 * @param fp_SNIa the SNIa file pointer
 * @param fp_r_processes the r_processes file pointer
 * @param us the used unit system
 * @param phys_const the physical constant struct
 */
INLINE static void feedback_logger_init_log_file(
    FILE *fp_SNII, FILE *fp_SNIa, FILE *fp_r_processes,
    const struct unit_system *restrict us,
    const struct phys_const *phys_const) {

  /* Initialize the SNII logger */
  feedback_logger_SNIa_init_log_file(fp_SNII, us, phys_const);
  fflush(fp_SNII);

  /* Initialize the SNIa logger */
  feedback_logger_SNIa_init_log_file(fp_SNIa, us, phys_const);
  fflush(fp_SNIa);

  /* Initialize the r_processes logger */
  feedback_logger_r_processes_init_log_file(fp_r_processes, us, phys_const);
  fflush(fp_SNIa);
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
INLINE static double feedback_logger_MPI(
    const struct engine e) {

  const int nodeID = e->nodeID;

  /* Send all the information, first make an array for all the int variables */
  int logger_ints_received[3];
  const int logger_ints_send[3] = {SNII->heating, SNIa->heating, r_processes->events};

  /* make an array of all the double variables */
  double logger_doubles_received[4];
  const double logger_doubles_send[4] = {SNII->SNII_energy, SNII->N_SNII,
                                   SNIa->SNIa_energy,
                                   r_processes->enrichment_mass};

  /* Do the MPI Reduce on the int array and the double array */
  MPI_Reduce(&logger_ints_send, &logger_ints_received, 3, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&logger_doubles_send, &logger_doubles_received, 4, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);

  /* Clear the variables or store them depending if we are node 0 or not */
  if (nodeID != 0) {
    feedback_logger_SNII_clear(SNII);
    feedback_logger_SNIa_clear(SNIa);
    feedback_logger_r_processes_clear(r_processes);
    return logger_time - delta_logger_time;
  }

  /* Get the values for the SNII */
  SNII->heating = logger_ints_received[0];
  SNII->SNII_energy = logger_doubles_received[0];
  SNII->N_SNII = logger_doubles_received[1];

  /* Get the values for the SNIa */
  SNIa->heating = logger_ints_received[1];
  SNIa->SNIa_energy = logger_doubles_received[2];

  /* Get the values for the r-processes */
  r_processes->events = logger_ints_received[2];
  r_processes->enrichment_mass = logger_doubles_received[3];

  /* Return the received time if we are node 0 */
  return logger_time;
}
#endif

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_H */
