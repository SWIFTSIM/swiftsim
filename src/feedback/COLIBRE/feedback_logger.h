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
    const struct feedback_history_accumulator *fha, const int step, const double time,
    const double a, const double z, const double volume) {

  if (step == 0) return;

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

INLINE static void feedback_logger_update_times(struct feedback_history_accumulator *fha, const int step, const double time,
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
    const double mass_init = pj->mass; // I first had this but this doesn't seem right: const double mass_init = si->mass_init;
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
    const struct feedback_history_accumulator *fha, const int step, const double time,
    const double a, const double z, const double volume) {

  if (step == 0) return;

  /* Calculate Delta time */
  const double delta_time = time - fha->time_prev;

  /* Get the total amount of SNIa energy */
  const double E_SNII = SNII->SNII_energy;

  /* Get the Energy of a single SNIa */
  const double E_single_SNII = feedback_properties->E_SNII;

  /* Calculate the number of SNIas in the simulation */
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
    const double mass_init = pj->mass; // I first had this but this doesn't seem right: const double mass_init = si->mass_init;
    const double delta_u = si->feedback_data.to_distribute.SNII_delta_u;
    const double deltaE = delta_u * mass_init;

    /* Update the total SNII energy */
    SNII->SNII_energy += deltaE;
    SNII->heating += 1;

    /* For the number of SNIIs first divide by the energy fraction, rest is done while writing the data */
    SNII->N_SNII += deltaE/f_E;
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

  /* Calculate the energy unit */
  const double E_unit = us->UnitMass_in_cgs * us->UnitLength_in_cgs *
                        us->UnitLength_in_cgs /
                        (us->UnitTime_in_cgs * us->UnitTime_in_cgs);

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
 * @brief log all the collected data to the r-processes logger file
 *
 * @param feedback_properties, the feedback structure
 * @param fp the file pointer to write to.
 * @param r_process the external variable that stores all the feedback information
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
    const double a, const double z, const double volume) {}

/**
 * @brief function to reinitialize the SNII logger.
 *
 * @param SNII the external variable that stores all the feedback information
 */
INLINE static void feedback_logger_r_processes_clear(
    struct feedback_history_r_processes *restrict r_process) {}

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
    feedback_logger_SNIa_log_data(feedback_properties, fp_SNIa, SNIa, fha, step, time, a, z, volume);
    feedback_logger_SNIa_clear(SNIa);
    fflush(fp_SNIa);

    /* Log the SNII data */
    feedback_logger_SNII_log_data(feedback_properties, fp_SNII, SNII, fha, step, time, a, z, volume);
    feedback_logger_SNII_clear(SNII);
    fflush(fp_SNII);

    /* Log the r_processes data */
    feedback_logger_r_processes_log_data(feedback_properties, fp_r_processes, r_processes, fha, step, time, a, z, volume);
    feedback_logger_r_processes_clear(r_processes);
    fflush(fp_r_processes);
  
    /* Update the times in the gneral logger struct */ 
    feedback_logger_update_times(fha, step, time, a, z);

    }

#ifdef WITH_MPI
INLINE static double feedback_logger_MPI(int nodeID, struct feedback_history_SNII *restrict SNII, struct feedback_history_SNIa *restrict SNIa, struct feedback_history_r_processes *restrict r_processes, double logger_time, double delta_logger_time) {
  
  /* Send around the SNII information */
  int total_SNII_events;
  double global_SNII_info[2];
  double send_SNII = [SNII->SNII_energy, SNII->N_SNII];
  
  /* Do the MPI reduce communication for the SNII logger */
  MPI_Reduce(&send_SNII, &global_SNII_info, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNII->heating, &total_SNII_events, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  /* Clear or store depending if we are node 0 or not */
  if (nodeID==0) {
    SNII->heating = total_SNII_events;
    SNII->SNII_energy = global_SNII_info[0];
    SNII->N_SNII = global_SNII_info[1];
  } else {
    feedback_logger_SNII_clear(SNII);
  }

  /* Send around the r-processes information */
  int total_r_processes_events;
  double total_r_processes_mass;
  
  /* Do the MPI reduce communication for the r-processes logger */
  MPI_Reduce(&r_processes->enrichement_mass, &total_r_processes_mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&r_processes->events, &total_r_processes_events, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  /* Clear or store depending if we are node 0 or not */
  if (nodeID==0) {
    r_processes->events = total_r_processes_events;
    r_processes->enrichement_mass = total_r_processes_mass;
  } else {
    feedback_logger_r_processes_clear(r_processes);
  }
  
  /* Send around the SNIa information */
  int total_SNIa_events;
  double global_SNIa_info;

  /* Send around the SNIa information */
  MPI_Reduce(&SNIa->heating, &total_SNIa_events, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SNIa->SNIa_energy, &global_SNIa_info, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* Clear or store depending if we are node 0 or not */
  if (nodeID==0) {
    SNIa->heating = total_SNIa_events;
    SNIa->SNIa_energy = global_SNIa_info;
    return logger_time;
  } else {
    feedback_logger_SNIa_clear(SNIa);
    return logger_time - delta_logger_time;
  }
}
#endif

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_H */
