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
}

INLINE static void feedback_logger_r_processes_init(const struct engine *restrict e) {
 
  /* Initialize the core variables */ 
  feedback_logger_core_init(e, &log_r_processes.core,3);

  /* Initialize the energy to zero */
  log_r_processes.enrichment_mass = 0.;

  /* Initialize the number of heating events to zero */
  log_r_processes.events = 0;
 
}

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_R_PROCESSES_H */
