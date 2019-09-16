/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COLIBRE_FEEDBACK_IACT_H
#define SWIFT_COLIBRE_FEEDBACK_IACT_H

/* Local includes */
#include "feedback_logger.h"
#include "random.h"

/* Define external variables */
extern FILE *SNIa_logger;

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float *dx,
                                    const float hi, const float hj,
                                    struct spart *restrict si,
                                    const struct part *restrict pj,
                                    const struct xpart *restrict xpj,
                                    const struct cosmology *restrict cosmo,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.to_collect.ngb_mass += mj;

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  const float rho = hydro_get_comoving_density(pj);
  if (rho != 0.f)
    si->feedback_data.to_collect.enrichment_weight_inv += wi / rho;
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 * @param time current physical time in the simulation
 * @param step current step counter
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float *dx, const float hi, const float hj,
    const struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, const struct cosmology *restrict cosmo,
    const integertime_t ti_current, const double time, const int step) {

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Gas particle density */
  const float rho_j = hydro_get_comoving_density(pj);

  /* Compute weighting for distributing feedback quantities */
  float Omega_frac;
  if (rho_j != 0.f) {
    Omega_frac = si->feedback_data.to_distribute.enrichment_weight * wi / rho_j;
  } else {
    Omega_frac = 0.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (Omega_frac < 0. || Omega_frac > 1.00001)
    error("Invalid fraction of material to distribute. Omega_frac=%e",
          Omega_frac);
#endif

  /* Update particle mass */
  const double current_mass = hydro_get_mass(pj);
  const double delta_mass = si->feedback_data.to_distribute.mass * Omega_frac;
  const double new_mass = current_mass + delta_mass;

  hydro_set_mass(pj, new_mass);

  /* Inverse of the new mass */
  const double new_mass_inv = 1. / new_mass;

  /* Update total metallicity */
  const double current_metal_mass_total =
      pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const double delta_metal_mass_total =
      si->feedback_data.to_distribute.total_metal_mass * Omega_frac;
  const double new_metal_mass_total =
      current_metal_mass_total + delta_metal_mass_total;

  pj->chemistry_data.metal_mass_fraction_total =
      new_metal_mass_total * new_mass_inv;

  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const double delta_metal_mass =
        si->feedback_data.to_distribute.metal_mass[elem] * Omega_frac;
    const double new_metal_mass = current_metal_mass + delta_metal_mass;

    pj->chemistry_data.metal_mass_fraction[elem] =
        new_metal_mass * new_mass_inv;
  }

  /* Update iron mass fraction from SNIa  */
  const double current_iron_from_SNIa_mass =
      pj->chemistry_data.iron_mass_fraction_from_SNIa * current_mass;
  const double delta_iron_from_SNIa_mass =
      si->feedback_data.to_distribute.Fe_mass_from_SNIa * Omega_frac;
  const double new_iron_from_SNIa_mass =
      current_iron_from_SNIa_mass + delta_iron_from_SNIa_mass;

  pj->chemistry_data.iron_mass_fraction_from_SNIa =
      new_iron_from_SNIa_mass * new_mass_inv;

  /* Update mass fraction from SNIa  */
  const double current_mass_from_SNIa =
      pj->chemistry_data.mass_from_SNIa * current_mass;
  const double delta_mass_from_SNIa =
      si->feedback_data.to_distribute.mass_from_SNIa * Omega_frac;
  const double new_mass_from_SNIa =
      current_mass_from_SNIa + delta_mass_from_SNIa;

  pj->chemistry_data.mass_from_SNIa = new_mass_from_SNIa * new_mass_inv;

  /* Update metal mass fraction from SNIa */
  const double current_metal_mass_from_SNIa =
      pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  const double delta_metal_mass_from_SNIa =
      si->feedback_data.to_distribute.metal_mass_from_SNIa * Omega_frac;
  const double new_metal_mass_from_SNIa =
      current_metal_mass_from_SNIa + delta_metal_mass_from_SNIa;

  pj->chemistry_data.metal_mass_fraction_from_SNIa =
      new_metal_mass_from_SNIa * new_mass_inv;

  /* Update mass fraction from SNII  */
  const double current_mass_from_SNII =
      pj->chemistry_data.mass_from_SNII * current_mass;
  const double delta_mass_from_SNII =
      si->feedback_data.to_distribute.mass_from_SNII * Omega_frac;
  const double new_mass_from_SNII =
      current_mass_from_SNII + delta_mass_from_SNII;

  pj->chemistry_data.mass_from_SNII = new_mass_from_SNII * new_mass_inv;

  /* Update metal mass fraction from SNII */
  const double current_metal_mass_from_SNII =
      pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  const double delta_metal_mass_from_SNII =
      si->feedback_data.to_distribute.metal_mass_from_SNII * Omega_frac;
  const double new_metal_mass_from_SNII =
      current_metal_mass_from_SNII + delta_metal_mass_from_SNII;

  pj->chemistry_data.metal_mass_fraction_from_SNII =
      new_metal_mass_from_SNII * new_mass_inv;

  /* Update mass fraction from AGB  */
  const double current_mass_from_AGB =
      pj->chemistry_data.mass_from_AGB * current_mass;
  const double delta_mass_from_AGB =
      si->feedback_data.to_distribute.mass_from_AGB * Omega_frac;
  const double new_mass_from_AGB = current_mass_from_AGB + delta_mass_from_AGB;

  pj->chemistry_data.mass_from_AGB = new_mass_from_AGB * new_mass_inv;

  /* Update metal mass fraction from AGB */
  const double current_metal_mass_from_AGB =
      pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  const double delta_metal_mass_from_AGB =
      si->feedback_data.to_distribute.metal_mass_from_AGB * Omega_frac;
  const double new_metal_mass_from_AGB =
      current_metal_mass_from_AGB + delta_metal_mass_from_AGB;

  pj->chemistry_data.metal_mass_fraction_from_AGB =
      new_metal_mass_from_AGB * new_mass_inv;

  /* Compute the current kinetic energy */
  const double current_v2 = xpj->v_full[0] * xpj->v_full[0] +
                            xpj->v_full[1] * xpj->v_full[1] +
                            xpj->v_full[2] * xpj->v_full[2];
  const double current_kinetic_energy_gas =
      0.5 * cosmo->a2_inv * current_mass * current_v2;

  /* Compute the current thermal energy */
  const double current_thermal_energy =
      current_mass * hydro_get_physical_internal_energy(pj, xpj, cosmo);

  /* Apply conservation of momentum */

  /* Update velocity following change in gas mass */
  xpj->v_full[0] *= current_mass * new_mass_inv;
  xpj->v_full[1] *= current_mass * new_mass_inv;
  xpj->v_full[2] *= current_mass * new_mass_inv;

  /* Update velocity following addition of mass with different momentum */
  xpj->v_full[0] += delta_mass * new_mass_inv * si->v[0];
  xpj->v_full[1] += delta_mass * new_mass_inv * si->v[1];
  xpj->v_full[2] += delta_mass * new_mass_inv * si->v[2];

  /* Compute the new kinetic energy */
  const double new_v2 = xpj->v_full[0] * xpj->v_full[0] +
                        xpj->v_full[1] * xpj->v_full[1] +
                        xpj->v_full[2] * xpj->v_full[2];
  const double new_kinetic_energy_gas = 0.5 * cosmo->a2_inv * new_mass * new_v2;

  /* Energy injected
   * (kinetic energy of ejecta + kinetic energy of star AGB) */
  const double injected_energy =
      si->feedback_data.to_distribute.energy * Omega_frac;

  /* Apply energy conservation to recover the new thermal energy of the gas */
  const double new_thermal_energy = current_kinetic_energy_gas +
                                    current_thermal_energy + injected_energy -
                                    new_kinetic_energy_gas;

  /* Convert this to a specific thermal energy */
  const double u_new_enrich = new_thermal_energy * new_mass_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (new_thermal_energy < 0.99 * current_thermal_energy)
    error("Enrichment is cooling the gas");
#endif

  /* Do the energy injection. */
  hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new_enrich);
  hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new_enrich);

  /* Finally, SNII stochastic feedback */

  /* Get the SNII feedback properties */
  const float prob = si->feedback_data.to_distribute.SNII_heating_probability;

  /* Are we doing some SNII (big boys) feedback? */
  if (prob > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval_two_IDs(
        si->id, pj->id, ti_current, random_number_stellar_feedback);
    /* Are we lucky? */
    if (rand < prob) {

      /* Compute new energy of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u = si->feedback_data.to_distribute.SNII_delta_u;
      const double u_new = u_init + delta_u;

//#ifdef SWIFT_DEBUG_CHECKS
      message("SNII event at star age [Myr]  = %.4f", si->feedback_data.to_distribute.SNII_star_age_Myr);
//#endif

      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* message( */
      /*     "We did some heating! id %llu star id %llu probability %.5e " */
      /*     "random_num %.5e du %.5e du/ini %.5e", */
      /*     pj->id, si->id, prob, rand, delta_u, delta_u / u_init); */
    }
  }

  /* Finally, SNIa stochastic feedback */

  /* Get the SNIa feedback properties */
  const float prob_SNIa =
      si->feedback_data.to_distribute.SNIa_heating_probability;

  /* Are we doing some SNIa (old boys) feedback? */
  if (prob_SNIa > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand_SNIa = random_unit_interval_two_IDs(
        si->id, pj->id, ti_current, random_number_SNIa_feedback);
    /* Are we lucky? */
    if (rand_SNIa < prob_SNIa) {

      /* Compute new energy of this particle */
      const double u_init_SNIa =
          hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u_SNIa = si->feedback_data.to_distribute.SNIa_delta_u;
      const double u_new = u_init_SNIa + delta_u_SNIa;

      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Write the event to the SNIa logger file */
      feedback_SNIa_logger_write_to_log_file(SNIa_logger, time, si, pj, xpj,
                                             cosmo, step);
      fflush(SNIa_logger);
    }
  }

  /* Kick gas particle away from the star using the momentum available in the
   * timestep.This is done stochastically.
   * However, if delta_v is small enough (or even negative), this translates
   * into kicking all neighboring particles away from the star. */

  const float delta_v = si->feedback_data.to_distribute.momentum_delta_v;
  const float momentum_prob =
      si->feedback_data.to_distribute.momentum_probability;

  const float HIIregion_prob = 
      si->feedback_data.to_distribute.HIIregion_probability;

  /* Draw a random number (Note mixing both IDs) */
  const float momentum_rand = random_unit_interval_two_IDs(
      si->id, pj->id, ti_current, random_number_stellar_winds);

  /* if lucky, perform the actual kick  */
  if (momentum_rand < momentum_prob) {

    /* Note that xpj->v_full = a^2 * dx/dt, with x the comoving coordinate.
     * Therefore, a physical kick, dv, gets translated into a
     * code velocity kick, a * dv */

    xpj->v_full[0] -= delta_v * dx[0] * r_inv * cosmo->a;
    xpj->v_full[1] -= delta_v * dx[1] * r_inv * cosmo->a;
    xpj->v_full[2] -= delta_v * dx[2] * r_inv * cosmo->a;

    /* Store how much physical momentum is received, so we don't care about
     * cosmology */
    xpj->tracers_data.momentum_received += delta_v * current_mass;
  }

  /* Draw a random number (Note mixing both IDs) */
  const float HIIregion_rand = random_unit_interval(
      si->id + pj->id, ti_current, random_number_HII_regions);

  /* if lucky, particle is now flagged as HII region  */
  if (HIIregion_rand < HIIregion_prob) {
   /* gas particle gets flagged as HII region */
    xpj->tracers_data.HIIregion_timer_gas = si->feedback_data.to_distribute.HIIregion_endtime;
    xpj->tracers_data.HIIregion_starid = si->feedback_data.to_distribute.HIIregion_starid;

   /* Impose maximal viscosity */
   /* hydro_diffusive_feedback_reset(pj); */
  }
}

#endif /* SWIFT_COLIBRE_FEEDBACK_IACT_H */
