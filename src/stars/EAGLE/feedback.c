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
  // Constant is for Chabrier IMF between 6 and 100 Msun.
  return sp->mass_init * 0.017362 * 1e10;
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

  /* Model parameters */
  const double f_E_max = props->feedback.f_E_max;
  const double f_E_min = props->feedback.f_E_min;
  const double Z_0 = props->feedback.Z_0;
  const double n_0 = props->feedback.n_0_cgs;
  const double n_Z = props->feedback.n_Z;
  const double n_n = props->feedback.n_n;

  /* Star properties */

  /* Smoothed metallicity (metal mass fraction) at birth time of the star */
  const double Z_smooth = sp->chemistry_data.smoothed_metal_mass_fraction_total;

  /* Physical density of the gas at the star's birth time */
  const double rho_birth = sp->birth_density;
  const double n_birth = rho_birth * props->feedback.conv_factor_rho_to_n_cgs;

  /* Calculate f_E */
  const double Z_term = pow(max(Z_smooth, 1e-6) / Z_0, n_Z);
  const double n_term = pow(n_birth / n_0, -n_n);
  const double denonimator = 1. + Z_term * n_term;

  return f_E_min + (f_E_max - f_E_min) / denonimator;
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
void stars_prepare_feedback(
    struct spart* sp, const struct stars_props* star_props,
    const struct hydro_props* hydro_props, const struct unit_system* us,
    const struct phys_const* phys_const, const struct cosmology* cosmo,
    const int with_cosmology, const integertime_t ti_current,
    const double time_base) {

  /* Skip particles that were in the ICs. */
  if (sp->birth_time < 0. || sp->feedback.probability == -1) return;

  /* Star and end of the time-step, star's birth time (all internal units) */
  double current_time, birth_time, old_time;
  if (with_cosmology) {
    current_time = cosmo->time;
    const integertime_t ti_begin =
        ti_current - get_integer_timestep(sp->time_bin);
    const double delta_time =
        cosmology_get_delta_time(cosmo, ti_begin, ti_current);
    old_time = current_time - delta_time;
    birth_time =
        cosmology_get_time_since_big_bang(cosmo, sp->birth_scale_factor);
  } else {
    current_time = cosmo->time;
    const integertime_t ti_begin = get_integer_timestep(sp->time_bin);
    const double delta_time = (ti_current - ti_begin) * time_base;
    old_time = current_time - delta_time;
    birth_time = sp->birth_time;
  }

  /* Is it time to go supernova? */
  if ((current_time - birth_time >= star_props->feedback.SNII_feedback_delay) &&
      (old_time - birth_time < star_props->feedback.SNII_feedback_delay)) {

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
    prob *= (N_SNe * f_E) / (deltaT * sp->density.neighbour_mass);

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

    message("ID=%lld Probability: %e delta_u: %e birth_time: %e", sp->id, prob,
            delta_u, sp->birth_time);

    /* Store all of this in the star particle for use in the feedback loop */
    sp->feedback.probability = prob;
    sp->feedback.delta_u = delta_u;
    sp->feedback.ti_current = ti_current;
  }
}
