/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_NEUTRINO_IO_H
#define SWIFT_DEFAULT_NEUTRINO_IO_H

#include "../../gravity.h"

/* Local includes */
#include "fermi_dirac.h"
#include "neutrino.h"
#include "neutrino_properties.h"

INLINE static void convert_gpart_pi(const struct engine* e,
                                    const struct gpart* gp, float* ret) {

  /* When we are running with the delta-f method, resample the momentum */
  if (e->neutrino_properties->use_delta_f) {
    /* Use a particle id dependent seed (sum of global seed and ID) */
    const long long neutrino_seed = e->neutrino_properties->neutrino_seed;
    const long long seed = gp->id_or_neg_offset + neutrino_seed;

    ret[0] = neutrino_seed_to_fermi_dirac(seed);  // eV
  } else {
    /* We don't know what the initial momentum was and we don't need it */
    ret[0] = 0.f;
  }
}

INLINE static void convert_gpart_mnu(const struct engine* e,
                                     const struct gpart* gp, float* ret) {

  /* When we are running with the delta-f method, resample the mass */
  if (e->neutrino_properties->use_delta_f) {

    /* Use a particle id dependent seed (sum of global seed and ID) */
    const long long neutrino_seed = e->neutrino_properties->neutrino_seed;
    const long long seed = gp->id_or_neg_offset + neutrino_seed;

    /* Fetch neutrino masses defined in the cosmology */
    const int N_nu = e->cosmology->N_nu;
    const double* m_eV_array = e->cosmology->M_nu_eV;

    ret[0] = neutrino_seed_to_mass(N_nu, m_eV_array, seed);  // eV
  } else {
    /* Otherwise, simply use the mass implied by the conversion factor */
    const double mass_factor = e->neutrino_mass_conversion_factor;

    ret[0] = gp->mass * mass_factor;  // eV
  }
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param gparts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int neutrino_write_particles(
    const struct gpart* gparts, struct io_props* list) {

  list[0] = io_make_output_field_convert_gpart(
      "InitialMomenta", FLOAT, 1, UNIT_CONV_NO_UNITS, 2.f, gparts,
      convert_gpart_pi,
      "Initial Fermi-Dirac momenta at infinity in electron-volts");

  list[1] = io_make_output_field_convert_gpart(
      "MicroscopicMasses", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, gparts,
      convert_gpart_mnu, "The microscopic neutrino masses in electron-volts");

  return 2;
}

#endif /* SWIFT_DEFAULT_NEUTRINO_IO_H */
