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

/**
 * @brief Initialize the SFH logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void feedback_SNIa_logger_init_log_file(
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

INLINE static void feedback_SNIa_logger_write_to_log_file(
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
