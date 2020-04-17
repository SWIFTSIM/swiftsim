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
#include "event_logger.h"
#include "random.h"
#include "timestep_sync_part.h"
#include "tracers.h"

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
                                    struct spart *si, const struct part *pj,
                                    const struct xpart *xpj,
                                    const struct cosmology *cosmo,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r. */
  const float r = sqrtf(r2);

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
 * @param xpj Extra particle data.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 * @param time current physical time in the simulation
 * @param step current step counter
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  const struct spart *si, struct part *pj,
                                  struct xpart *xpj, const int with_cosmology,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current,
                                  const double time, const int step) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

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
  if (Omega_frac < 0. || Omega_frac > 1.01)
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

  /* Calculate mean metal weighted redshift */
  double delta_mass_times_time, delta_iron_mass_times_time;
  if (with_cosmology) {
    delta_mass_times_time = delta_metal_mass_total * cosmo->z;
  } else {
    delta_mass_times_time = delta_metal_mass_total * time;
  }
  /* Update metal mass tracker */
  pj->chemistry_data.metal_mass_tracker += delta_mass_times_time;

  /* Calculate mean metal weighted redshift */
  if (new_metal_mass_total > 0.f) {
    pj->chemistry_data.metal_weighted_redshift =
        pj->chemistry_data.metal_mass_tracker / new_metal_mass_total;
  }

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

  /* Calculate mean iron weighted redshift */
  if (with_cosmology) {
    delta_iron_mass_times_time = delta_iron_from_SNIa_mass * cosmo->z;
  } else {
    delta_iron_mass_times_time = delta_iron_from_SNIa_mass * time;
  }
  /* Update iron mass tracker */
  pj->chemistry_data.iron_mass_tracker += delta_iron_mass_times_time;

  /* Calculate mean iron weighted redshift */
  if (new_iron_from_SNIa_mass > 0.f) {
    pj->chemistry_data.iron_weighted_redshift =
        pj->chemistry_data.iron_mass_tracker / new_iron_from_SNIa_mass;
  }

  /* Update mass from SNIa  */
  const double delta_mass_from_SNIa =
      si->feedback_data.to_distribute.mass_from_SNIa * Omega_frac;

  pj->chemistry_data.mass_from_SNIa += delta_mass_from_SNIa;

  /* Update metal mass fraction from SNIa */
  const double current_metal_mass_from_SNIa =
      pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  const double delta_metal_mass_from_SNIa =
      si->feedback_data.to_distribute.metal_mass_from_SNIa * Omega_frac;
  const double new_metal_mass_from_SNIa =
      current_metal_mass_from_SNIa + delta_metal_mass_from_SNIa;

  pj->chemistry_data.metal_mass_fraction_from_SNIa =
      new_metal_mass_from_SNIa * new_mass_inv;

  /* Update mass from SNII  */
  const double delta_mass_from_SNII =
      si->feedback_data.to_distribute.mass_from_SNII * Omega_frac;

  pj->chemistry_data.mass_from_SNII += delta_mass_from_SNII;

  /* Update metal mass fraction from SNII */
  const double current_metal_mass_from_SNII =
      pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  const double delta_metal_mass_from_SNII =
      si->feedback_data.to_distribute.metal_mass_from_SNII * Omega_frac;
  const double new_metal_mass_from_SNII =
      current_metal_mass_from_SNII + delta_metal_mass_from_SNII;

  pj->chemistry_data.metal_mass_fraction_from_SNII =
      new_metal_mass_from_SNII * new_mass_inv;

  /* Update mass from AGB  */
  const double delta_mass_from_AGB =
      si->feedback_data.to_distribute.mass_from_AGB * Omega_frac;

  pj->chemistry_data.mass_from_AGB += delta_mass_from_AGB;

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

  /* Apply energy conservation to recover the new thermal energy of the gas
   * Note: in some specific cases the new_thermal_energy could be lower
   * than the current_thermal_energy, this is mainly the case if the change
   * in mass is relatively small and the velocity vectors between both the
   * gas particle and the star particle have a small angle. */
  const double new_thermal_energy = current_kinetic_energy_gas +
                                    current_thermal_energy + injected_energy -
                                    new_kinetic_energy_gas;

  /* Convert this to a specific thermal energy */
  const double u_new_enrich = new_thermal_energy * new_mass_inv;

  /* Do the energy injection. */
  hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new_enrich);
  hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new_enrich);

  /* Apply the different sources of energy and momentum feedback */

  int do_SNII_thermal = 0;
  int do_SNII_kinetic = 0;
  int do_SNIa = 0;

  /* SNII stochastic thermal feedback */

  /* Get the SNII thermal feedback properties */
  const float prob_SNII_thermal =
      si->feedback_data.to_distribute.SNII_heating_probability;

  /* Are we doing some SNII (big boys) thermal feedback? */
  if (prob_SNII_thermal > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand_SNII = random_unit_interval_two_IDs(
        si->id, pj->id, ti_current, random_number_stellar_feedback_1);

    /* Are we lucky? */
    do_SNII_thermal = (rand_SNII < prob_SNII_thermal);

    if (do_SNII_thermal) {

      /* Compute new energy of this particle */
      const double u_init_SNII =
          hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u_SNII = si->feedback_data.to_distribute.SNII_delta_u;
      const double u_new_SNII = u_init_SNII + delta_u_SNII;

#ifdef SWIFT_DEBUG_CHECKS
      message("SNII event at star age [Myr]  = %.4f",
              si->feedback_data.to_distribute.SNII_star_age_Myr);
#endif

      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new_SNII);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new_SNII);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Update cooling properties. */
      cooling_update_feedback_particle(xpj);

      /* Mark this particle has having been heated by supernova feedback */
      tracers_after_SNII_feedback(xpj, with_cosmology, cosmo->a, time);

      /* Write the event to the SNIII log file */
      event_logger_SNII_log_event(si, pj, xpj, cosmo, si->SNII_f_E);

#ifdef SWIFT_DEBUG_CHECKS
      event_logger_SNII_log_event_debug(time, si, pj, xpj, cosmo, step);
#endif
      /* message( */
      /*     "We did some heating! id %llu star id %llu probability %.5e " */
      /*     "random_num %.5e du %.5e du/ini %.5e", */
      /*     pj->id, si->id, prob, rand, delta_u, delta_u / u_init); */

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }

  /* SNII stochastic kinetic feedback */

  /* Get the SNII kinetic feedback properties */
  const float prob_SNII_kinetic =
      si->feedback_data.to_distribute.SNII_kick_probability;

  /* Are we doing some SNII (big boys) kinetic feedback? */
  if (prob_SNII_kinetic > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand_SNII = random_unit_interval_two_IDs(
        si->id, pj->id, ti_current, random_number_stellar_feedback_3);

    /* Are we lucky? */
    do_SNII_kinetic = (rand_SNII < prob_SNII_kinetic);

    if (do_SNII_kinetic) {

      /* Get the SNII feedback kick velocity */
      const float SNII_delta_v = si->feedback_data.to_distribute.SNII_delta_v;

      /* Note that xpj->v_full = a^2 * dx/dt, with x the comoving coordinate.
       * Therefore, a physical kick, dv, gets translated into a
       * code velocity kick, a * dv */

      xpj->v_full[0] -= SNII_delta_v * dx[0] * r_inv * cosmo->a;
      xpj->v_full[1] -= SNII_delta_v * dx[1] * r_inv * cosmo->a;
      xpj->v_full[2] -= SNII_delta_v * dx[2] * r_inv * cosmo->a;

      /* Update the signal velocity of the particle based on the velocity kick
       */
      hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, SNII_delta_v);

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }

  /* SNIa stochastic thermal feedback */

  /* Get the SNIa feedback properties */
  const float prob_SNIa =
      si->feedback_data.to_distribute.SNIa_heating_probability;

  /* Are we doing some SNIa feedback? */
  if (prob_SNIa > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand_SNIa = random_unit_interval_two_IDs(
        si->id, pj->id, ti_current, random_number_stellar_feedback_2);
    /* Are we lucky? */
    do_SNIa = (rand_SNIa < prob_SNIa);

    if (do_SNIa) {

      /* Compute new energy of this particle */
      const double u_init_SNIa =
          hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u_SNIa = si->feedback_data.to_distribute.SNIa_delta_u;
      const double u_new_SNIa = u_init_SNIa + delta_u_SNIa;

      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new_SNIa);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new_SNIa);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Update cooling properties. */
      cooling_update_feedback_particle(xpj);

      /* Mark this particle has having been heated by supernova feedback */
      tracers_after_SNIa_feedback(xpj, with_cosmology, cosmo->a, time);

      /* Write the event to the SNIa log file */
      event_logger_SNIa_log_event(si, pj, xpj, cosmo);

#ifdef SWIFT_DEBUG_CHECKS
      event_logger_SNIa_log_event_debug(time, si, pj, xpj, cosmo, step);
#endif

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }

  /* Kick gas particle away from the star using the momentum available in the
   * timestep.This is done stochastically.
   * However, if delta_v is small enough (or even negative), this translates
   * into kicking all neighboring particles away from the star. */

  const float delta_v = si->feedback_data.to_distribute.momentum_delta_v;
  const float momentum_prob =
      si->feedback_data.to_distribute.momentum_probability;

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

    /* Mark this particle has having been heated by supernova feedback */
    tracers_after_momentum_feedback(xpj, with_cosmology, cosmo->a, time);

    /* Update the signal velocity of the particle based on the velocity kick */
    hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, delta_v);

    /* Synchronize the particle on the timeline */
    timestep_sync_part(pj);
  }

  /* Put particles into HII regions */

  const float HIIregion_prob =
      si->feedback_data.to_distribute.HIIregion_probability;

  /* Draw a random number (Note mixing both IDs) */
  const float HIIregion_rand = random_unit_interval_two_IDs(
      si->id, pj->id, ti_current, random_number_HII_regions);

  /* if lucky, particle is now flagged as HII region  */
  if (HIIregion_rand < HIIregion_prob) {

    /* gas particle gets flagged as HII region */
    xpj->tracers_data.HIIregion_timer_gas =
        si->feedback_data.to_distribute.HIIregion_endtime;
    xpj->tracers_data.HIIregion_starid =
        si->feedback_data.to_distribute.HIIregion_starid;

    /* If the particle hasn't received feedback from other
       sources, heat it up to the HII region energy level */
    if (!do_SNII_thermal && !do_SNIa) {

      /* Compute new energy of this particle */
      const double u_init_HII =
          hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float u_HII = si->feedback_data.to_distribute.HII_u;

      /* Put the particle on the HII temperature floor or leave it
         if it was already above */
      if (u_init_HII < u_HII) {

        /* Inject energy into the particle */
        hydro_set_physical_internal_energy(pj, xpj, cosmo, u_HII);
        hydro_set_drifted_physical_internal_energy(pj, cosmo, u_HII);

        /* Make sure the particle does not cool any more */
        hydro_set_physical_internal_energy_dt(pj, cosmo, 0.f);
      }

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }
}

#endif /* SWIFT_COLIBRE_FEEDBACK_IACT_H */
