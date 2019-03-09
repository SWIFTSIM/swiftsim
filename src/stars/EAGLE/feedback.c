/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @file src/stars/EAGLE/feedback.c
 * @brief EAGLE functions related to stellar feedback.
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes */
#include "adiabatic_index.h"
#include "cosmology.h"
#include "error.h"
#include "stars.h"
#include "stars_part.h"

/**
 * @brief Return the change in temperature (in Kelvin) to apply to a
 * gas particle affected by SNII feedback.
 *
 * @param props The properties of the stellar model.
 */
double eagle_feedback_temperature_change(const struct stars_props* props) {

  /* In the EAGLE model, the change of temperature is constant */
  return props->feedback.delta_T;
}

/**
 * @brief Computes the number of super-novae exploding for a given
 * star particle.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNe(const struct spart* sp,
                                    const struct stars_props* props) {

  // MATTHIEU: Add IMF integration!
  return sp->mass_init * 1.;
}

/**
 * @brief Computes the fraction of the available super-novae energy to
 * inject for a given event.
 *
 * Note that the fraction can be > 1.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 */
double eagle_feedback_energy_fraction(const struct spart* sp,
                                      const struct stars_props* props) {

  /* const double f_E_max = props->feedback.f_E_max; */
  /* const double f_E_max = props->feedback.f_E_max; */
  /* const double sigma = props->feedback.sigma; */
  /* const double Z_0 = props->feedback.Z_0; */
  /* const double rho_0 = props->feedback.rho_0; */
  /* const double n_n = props->feedback.n_n; */

  // MATTHIEU: Add full EAGLE model.
  return 1.;
}

/**
 * @brief Compute the stellar properties required for feedback after
 * the density loop has finished.
 *
 * In EAGLE this function computes the energy to inject per event as
 * well as the probability to do so.
 *
 * @param sp The particle to act upon
 * @param stars_properties The #stars_props
 * @param hydro_properties The properties of the hydro scheme.
 * @param us The current system of units.
 * @param phys_const The physical constants in the internal system of units.
 * @param cosmo The current cosmological model.
 */
void stars_prepare_feedback(struct spart* sp,
                            const struct stars_props* star_props,
                            const struct hydro_props* hydro_props,
                            const struct unit_system* us,
                            const struct phys_const* phys_const,
                            const struct cosmology* cosmo) {

  /* Skip particles that were in the ICs. */
  if (sp->birth_time < 0.) return;

  /* Properties of the model (all in internal units) */
  const double deltaT = eagle_feedback_temperature_change(star_props);
  const double N_SNe = eagle_feedback_number_of_SNe(sp, star_props);
  const double f_E = eagle_feedback_energy_fraction(sp, star_props);
  const double E_SNe = star_props->feedback.E_SNe;

  /* Some constants (all in internal units) */
  const double mu_ionised = hydro_props->mu_ionised;
  const double m_p = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  /* Calculate the heating probability */
  double prob = E_SNe * mu_ionised * hydro_gamma_minus_one * m_p / k_B;
  prob *= N_SNe * f_E / (deltaT * sp->density.neighbour_mass);

  /* Calculate the change in internal energy of the gas particles that get
   * heated */
  double delta_u;
  if (prob > 1.) {

    /* Special case: we need to adjust the energy irrespective of the
       desired deltaT to ensure we inject all the available energy. */

    prob = 1.;

    delta_u = E_SNe * N_SNe / sp->density.neighbour_mass;

  } else {

    /* Normal case */
    delta_u = deltaT * k_B / (mu_ionised * hydro_gamma_minus_one * m_p);
  }

  message("Probability: %e delta_u: %e", prob, delta_u);

  /* Store all of this in the star particle for application in the feedback loop
   */
  sp->feedback.probability = prob;
  sp->feedback.delta_u = delta_u;
}
