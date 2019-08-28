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

#include "feedback_logger_struct.h"

/**
 * @brief Initialize the SFH logger file
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
  fprintf(fp, "# Stochastic SNIa Logger file\n");
  fprintf(fp, "######################################################\n");
  fprintf(fp, "# The quantities are all given in internal physical units!\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# (0) Simulation step\n");
  fprintf(fp, "# (1) Previous simulation step\n");
  fprintf(fp,
          "# (2) Time since Big Bang (cosmological run), Time since start of "
          "the simulation (non-cosmological run).\n");
  fprintf(fp,
          "# (3) Previous time since Big Bang (cosmological run), Time since start of "
          "the simulation (non-cosmological run).\n");
  fprintf(fp, "#     Unit = %e seconds\n", us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e yr or %e Myr\n", 1.f / phys_const->const_year,
          1.f / phys_const->const_year / 1e6);
  fprintf(fp, "# (4) Scale factor              (no unit)\n");
  fprintf(fp, "# (5) Previous scale factor     (no unit)\n");
  fprintf(fp, "# (6) Redshift                  (no unit)\n");
  fprintf(fp, "# (7) Previous Redshift         (no unit)\n");
  fprintf(fp, "# (8) Injected energy of SNIa events\n");
  fprintf(fp, "#     Unit = %e erg\n", E_unit);
  fprintf(fp, "#     Unit = %e x 10^51 erg\n", E_unit / 1e51);
  fprintf(fp, "# (9) Number of SNIa   (number, no unit)\n");
  fprintf(fp, "# (10) Number of SNIa per time (binned in time between current time and previous time).\n");
  fprintf(fp, "#     Unit = %e #/seconds\n", 1./us->UnitTime_in_cgs);
  fprintf(fp, "#     Unit = %e #/yr or %e #/Myr\n", phys_const->const_year,
          phys_const->const_year * 1e6);
  fprintf(fp, "# (11) Number of SNIa per time per volume (binned in time between current time and previous time).\n");
  fprintf(fp, "#     Unit = %e #/seconds/cm^3\n", 1./us->UnitTime_in_cgs/pow(us->UnitLength_in_cgs, 3));
  fprintf(fp, "#     Unit = %e #/yr/Mpc^3\n", phys_const->const_year * pow(phys_const->const_parsec*1e6,3));
  fprintf(fp, "#\n");
  fprintf(
      fp,
      "#  (0)      (1)         (2)              (3)           (4)           "
      " (5)         (6)          (7)          (8)            (9)         "
      "   (10)              (11) \n");
  fprintf(fp,
      "# step  prev. step      time          prev. time        a         "
      " prev a          z          prev z    injection E       Numb SNIa   "
      "   SNIa rate     SNIa rate / volume \n");
}

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

INLINE static void feedback_logger_SNIa_init(struct feedback_history_accumulator *fha, const double time, const double a, const double z) {
  
  fha->step_prev = 0;
  fha->time_prev = time;
  fha->z_prev = z;
  fha->a_prev = a;

}

INLINE static void feedback_logger_SNIa_log_data(FILE *fp, struct feedback_history_SNIa *restrict SNIa,
struct feedback_history_accumulator *fha, const int step,
const double time, const double a, const double z, const double volume) {

  if (step==0) return;

  const double delta_time = time - fha->time_prev;
  const double N_SNIa = SNIa->N_SNIa;
  const double N_SNIa_p_time = N_SNIa/delta_time; 
  const double E_SNIa = SNIa->SNIa_energy;
  const double N_SNIa_p_time_p_volume = N_SNIa_p_time / volume;

  fprintf(fp, "%7d %7d %16e %16e %12.7f %12.7f %12.7f %12.7f %14.7f %16.7f %18.7f %12.7f \n"
  , step, fha->step_prev, time, fha->time_prev, a,
  fha->a_prev, z, fha->z_prev, E_SNIa, N_SNIa, N_SNIa_p_time, N_SNIa_p_time_p_volume);

  /* Set the times to the new values */
  fha->step_prev = step;
  fha->time_prev = time;
  fha->a_prev = a;
  fha->z_prev = z;

  /* Reset the global struct */
  SNIa->N_SNIa = 0.;
  SNIa->SNIa_energy = 0.;

}

INLINE static void feedback_logger_SNIa_log_event(
    struct feedback_history_SNIa *restrict SNIa, const double time, const struct spart *restrict si,
    struct part *restrict pj, struct xpart *restrict xpj,
    const struct cosmology *restrict cosmo, const int step) {
  
  /* Get the injected energy */
  const double mass_init = si->mass_init;
  const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
  const double deltaE = delta_u * mass_init;
  
  SNIa->SNIa_energy += deltaE;
  SNIa->N_SNIa += deltaE * 1.9884e2;
  message("Event!!!! BAM!! E=%e N=%e", SNIa->SNIa_energy, SNIa->N_SNIa);
}

INLINE static void feedback_logger_SNIa_write_to_file(
    FILE *fp, 
    struct feedback_history_SNIa *restrict SNIa,
    const struct cosmology *restrict cosmo, const int step) {
  
}

INLINE static void feedback_logger_SNIa_log_event_debug(
    FILE *fp, const double time, const struct spart *restrict si,
    struct part *restrict pj, struct xpart *restrict xpj,
    const struct cosmology *restrict cosmo, const int step) {

  /* Get the times */
  const double a = cosmo->a;
  const double z = cosmo->z;

  /* Get the injected energy */
  const double mass_init = si->mass_init;
  const double delta_u = si->feedback_data.to_distribute.SNIa_delta_u;
  const double deltaE = delta_u * mass_init;

  fprintf(fp, "%6d %16e %12.7f %12.7f %14llu %14llu %16e %16e\n", step, time, a,
          z, si->id, pj->id, deltaE, deltaE * 1.9884e2);
}

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_H */
