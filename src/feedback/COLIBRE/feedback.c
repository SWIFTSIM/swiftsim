/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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

/* This file's header */

/* This file's header */
#include "feedback.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "event_logger.h"
#include "feedback_tables.h"
#include "hydro_properties.h"
#include "imf.h"
#include "inline.h"
#include "interpolate.h"
#include "physical_constants.h"
#include "timers.h"
#include "yield_tables.h"
#include "dust_properties.h"
#include "dust.h"

/**
 * @brief Return the change in temperature (in internal units) to apply to a
 * gas particle affected by SNII feedback.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
INLINE static double eagle_SNII_feedback_temperature_change(
    const struct spart* sp, const struct feedback_props* props) {

  /* In the EAGLE REF model, the change of temperature is constant */
  return props->SNII_deltaT_desired;
}

/**
 * @brief Computes the number of supernovae of type II exploding for a given
 * star particle.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 * @param min_dying_mass_Msun stellar mass that dies at the end of this timestep
 * @param max_dying_mass_Msun stellar mass that dies at the beginning of this
 * timestep
 */

INLINE static double colibre_feedback_number_of_SNII(
    const struct spart* sp, const struct feedback_props* fp,
    const float min_dying_mass_Msun, const float max_dying_mass_Msun) {

  /* Calculate how many supernovae have exploded in this timestep */

  /* minimum star mass that dies is between SNII min and max mass */
  if (min_dying_mass_Msun <= exp10(fp->log10_SNII_max_mass_msun) &&
      min_dying_mass_Msun > exp10(fp->log10_SNII_min_mass_msun)) {
    /* both minimum and maximum star masses that die are within the SNII
     * bounds*/
    /* standard case: integrate imf from min_dying_mass_Msun to
     * max_dying_mass_Msun */
    if (max_dying_mass_Msun <= exp10(fp->log10_SNII_max_mass_msun))
      return integrate_imf(log10(min_dying_mass_Msun),
                           log10(max_dying_mass_Msun),
                           eagle_imf_integration_no_weight, NULL, fp) *
             sp->mass_init * fp->mass_to_solar_mass;
    /* minimum star mass that dies is within SNII bounds but maximum star mass
     *      that dies this timestep is higher then the maximum SNII mass*/
    /* integrate imf from min_dying_mass_Msun to SNII_max_mass_msun */
    else
      return integrate_imf(log10(min_dying_mass_Msun),
                           fp->log10_SNII_max_mass_msun,
                           eagle_imf_integration_no_weight, NULL, fp) *
             sp->mass_init * fp->mass_to_solar_mass;
    /* minimum star mass that dies is below the SNII mass range */
  } else if (min_dying_mass_Msun <= exp10(fp->log10_SNII_min_mass_msun)) {
    /* maximum star mass that dies is within the SNII range */
    /* integrate imf from SNII_min_mass_msun to max_dying_mass_Msun*/
    if (max_dying_mass_Msun > exp10(fp->log10_SNII_min_mass_msun) &&
        max_dying_mass_Msun <= exp10(fp->log10_SNII_max_mass_msun))
      return integrate_imf(fp->log10_SNII_min_mass_msun,
                           log10(max_dying_mass_Msun),
                           eagle_imf_integration_no_weight, NULL, fp) *
             sp->mass_init * fp->mass_to_solar_mass;
    /* minimum star mass that dies this timestep is below minimum SNII mass,
     * maximum star mass that dies this timestep is above the maximum SNII
     * mass*/
    /* integrate imf over full SNII range, from SNII_min_mass_msun to
     * SNII_max_mass_msun */
    else if (max_dying_mass_Msun > exp10(fp->log10_SNII_min_mass_msun))
      return integrate_imf(fp->log10_SNII_min_mass_msun,
                           fp->log10_SNII_max_mass_msun,
                           eagle_imf_integration_no_weight, NULL, fp) *
             sp->mass_init * fp->mass_to_solar_mass;
    /* both minimum and maximum star mass that dies this timestep are below
     * the SNII mass range; no SNII anymore */
    else
      return 0.0;
    /* minimum star mass that dies this timestep is higher than the maximum SNII
     * mass; no SNII yet */
  } else
    return 0.;
}

/**
 * @brief Computes the number of supernovae of type II exploding for a given
 * star particle.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 */
INLINE static double eagle_feedback_number_of_SNII(
    const struct spart* sp, const struct feedback_props* props) {

  /* Note: For a Chabrier 2003 IMF and SNII going off between 6 and 100
   * M_sun, the first term is 0.017362 M_sun^-1 */
  return props->num_SNII_per_msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the fraction of the available super-novae energy to
 * inject for a given event.
 *
 * Note that the fraction can be > 1.
 *
 * We use equation 7 of Schaye et al. 2015.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
INLINE static double eagle_SNII_feedback_energy_fraction(
    const struct spart* sp, const struct feedback_props* props) {

  /* Model parameters */
  const double f_E_max = props->f_E_max;
  const double f_E_min = props->f_E_min;
  const double Z_0 = props->Z_0;
  const double n_0 = props->n_0_cgs;
  const double n_Z = props->n_Z;
  const double n_n = props->n_n;

  /* Star properties */

  /* Metallicity (metal mass fraction) at birth time of the star */
  const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* Physical density of the gas at the star's birth time */
  const double rho_birth = sp->sf_data.birth_density;
  const double n_birth = rho_birth * props->rho_to_n_cgs;

  /* Calculate f_E */
  const double Z_term = pow(max(Z, 1e-6) / Z_0, n_Z);
  const double n_term = pow(n_birth / n_0, -n_n);
  const double denonimator = 1. + Z_term * n_term;

  return f_E_min + (f_E_max - f_E_min) / denonimator;
}

/**
 * @brief Compute the properties of the SNII stochastic feedback energy
 * injection.
 *
 * Only does something if the particle reached the SNII age during this time
 * step.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in internal units.
 * @param dt Length of time-step in internal units.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel.
 * @param ngb_gas_N Total un-weighted number of gas particles in the star's
 * kernel
 * @param feedback_props The properties of the feedback model.
 * @param age of star particle at the beginning of the timestep
 * @param timestep in Gyr
 * @param min_dying_mass_Msun stellar mass that dies at the end of this timestep
 * @param max_dying_mass_Msun stellar mass that dies at the beginning of this
 * timestep
 * @param ti_begin Integer time value at the beginning of timestep
 */
INLINE static void compute_SNII_feedback(
    struct spart* sp, const double star_age, const double dt,
    const float ngb_gas_mass, const int ngb_gas_N,
    const struct feedback_props* feedback_props,
    const float min_dying_mass_Msun, const float max_dying_mass_Msun,
    const integertime_t ti_begin) {

  /* Time after birth considered for SNII feedback (internal units) */
  const double SNII_wind_delay = feedback_props->SNII_wind_delay;

  /* Are we doing feedback this step?
   * Note that since the ages are calculated using an interpolation table we
   * must allow some tolerance here*/
  /* If SNII_wind_delay < 0, then use timed feedback */
  if ((star_age <= SNII_wind_delay &&
       (star_age + 1.001 * dt) > SNII_wind_delay) ||
      SNII_wind_delay < 0.) {

    /* Make sure a star does not do feedback twice in the delay time case */
    if (sp->SNII_f_E != -1.f && SNII_wind_delay >= 0.) {
#ifdef SWIFT_DEBUG_CHECKS
      message("Star has already done feedback! sp->id=%lld age=%e d=%e", sp->id,
              star_age, dt);
#endif
      return;
    }

    /* Properties of the model (all in internal units) */
    const double delta_T =
        eagle_SNII_feedback_temperature_change(sp, feedback_props);
    const double N_SNe_eagle =
        eagle_feedback_number_of_SNII(sp, feedback_props);
    const double N_SNe_colibre = colibre_feedback_number_of_SNII(
        sp, feedback_props, min_dying_mass_Msun, max_dying_mass_Msun);
    const double f_kin = feedback_props->SNII_f_kinetic;
    double delta_v = feedback_props->SNII_delta_v;
    double N_SNe;

    /* Total mass ejected at this time-step by the stellar particle sp */
    const float M_ej = sp->feedback_data.to_distribute.mass;

    /* Total un-weighted mass in the star's kernel including the ejecta mass */
    const float ngb_gas_mass_new = ngb_gas_mass + M_ej;

    /* Number of SNe detonated during this time-step */
    if (SNII_wind_delay > 0.) {
      N_SNe = N_SNe_eagle;
    } else {
      N_SNe = N_SNe_colibre;
    }

    /* Energy per SN */
    const double E_SNe = feedback_props->E_SNII;
    const double f_E = eagle_SNII_feedback_energy_fraction(sp, feedback_props);

    /* Total SN energy released by all detonated SNe in the stellar particle sp
     * during this time-step */
    const double E_SN_total = f_E * E_SNe * N_SNe;

    /* Conversion factor from T to internal energy */
    const double conv_factor = feedback_props->temp_to_u_factor;

    /* Calculate the default heating and kick probabilities */
    double prob_thermal =
        (1.0 - f_kin) * E_SN_total / (conv_factor * delta_T * ngb_gas_mass_new);

    /* Note that in the denominator we have ngb_gas_mass * delta_v * delta_v
     * and not 0.5 ngb_gas_mass * delta_v * delta_v. That is because in our
     * method, if we have a kick event then we kick two particles instead of
     * one. This implies that we need twice as much energy and the probability
     * must be lowered by a factor of 2. Futhermore, unlike in the thermal
     * feedback here we do not take into account the ejecta mass becasue the
     * kicks are done before the mass transfer. That is, we take ngb_gas_mass
     * and not ngb_gas_mass_new. */
    double prob_kinetic =
        f_kin * E_SN_total / (ngb_gas_mass * delta_v * delta_v);

    /* Calculate the change in internal energy of the gas particles that get
     * heated */
    double delta_u;

    /* Total kinetic energy used to kick gas particles by this stellar particle
    at this time-step */
    double E_kinetic;

    /* Numbers of heating and kick events at this time step
    for this stellar particle. In each kick event we kick two particles in the
    oposite directions */
    int number_of_th_SN_events = 0;
    int number_of_kin_SN_events = 0;

    /* If the heating probability is less than 1, compute the number of heating
     * events by drawing a random number ngb_gas_N times where ngb_gas_N
     * is the number of gas particles within the star's kernel */
    if (prob_thermal <= 1.) {

      /* Normal case */
      delta_u = delta_T * conv_factor;

      for (int i = 0; i < ngb_gas_N; i++) {
        const float rand_thermal = random_unit_interval_star_ID_and_ray_idx(
            sp->id, i, ti_begin, random_number_stellar_feedback_1);

        if (rand_thermal < prob_thermal) number_of_th_SN_events++;
      }

      /* If the probability is larger than or equal to 1 then adjust the energy
       */
    } else {

      /* Every gas neighbour is heated */
      number_of_th_SN_events = ngb_gas_N;

      /* Special case: we need to adjust the thermal energy per unit mass
       * irrespective of the desired delta T to ensure we inject all the
       * available SN energy. */

      prob_thermal = 1.;
      delta_u = (1.0 - f_kin) * E_SN_total / ngb_gas_mass_new;
    }

    /* IMPORTANT. If we have more heating events than the maximum number of
     * rays (colibre_feedback_number_of_rays), then obviously we cannot
     * distribute all of the heating events (since 1 event = 1 ray), so we need
     * to increase the thermal energy per ray and make the number of events
     * equal to the number of rays */
    if (number_of_th_SN_events > colibre_feedback_number_of_rays) {
      const double alpha_thermal = (double)number_of_th_SN_events /
                                   (double)colibre_feedback_number_of_rays;
      delta_u *= alpha_thermal;
      number_of_th_SN_events = colibre_feedback_number_of_rays;
    }

    /* Repeat the above steps for kick events in SNII feedback */
    if (prob_kinetic <= 1.) {

      for (int i = 0; i < ngb_gas_N; i++) {
        const float rand_kinetic = random_unit_interval_star_ID_and_ray_idx(
            sp->id, i, ti_begin, random_number_stellar_feedback_2);

        if (rand_kinetic < prob_kinetic) number_of_kin_SN_events++;
      }

      /* Total kinetic energy needed = Kinetic energy to kick one pair of two
       * particles, each of mean mass ngb_gas_mass_new/ngb_gas_N, with the kick
       * velocity delta_v \times the number of kick events
       * ( = the number of pairs) */
      E_kinetic = ngb_gas_mass / ngb_gas_N * delta_v * delta_v *
                  number_of_kin_SN_events;
    } else {

      number_of_kin_SN_events = ngb_gas_N;

      /* Special case: we need to adjust the kick velocity irrespective of the
       * desired delta v to ensure we inject all the available SN energy. */
      prob_kinetic = 1.;
      E_kinetic = f_kin * E_SN_total;
    }

    /* The number of kick events cannot be greater than the number of rays */
    number_of_kin_SN_events =
        min(number_of_kin_SN_events, colibre_feedback_number_of_rays);

#ifdef SWIFT_DEBUG_CHECKS
    if (f_E < feedback_props->f_E_min || f_E > feedback_props->f_E_max)
      error("f_E is not in the valid range! f_E=%f sp->id=%lld", f_E, sp->id);
#endif

    /* Store all of this in the star for delivery onto the gas */
    sp->SNII_f_E = f_E;
    sp->feedback_data.to_distribute.SNII_heating_probability = prob_thermal;
    sp->feedback_data.to_distribute.SNII_kick_probability = prob_kinetic;
    sp->feedback_data.to_distribute.SNII_delta_u = delta_u;
    sp->feedback_data.to_distribute.SNII_E_kinetic = E_kinetic;
    sp->feedback_data.to_distribute.SNII_number_of_heating_events =
        number_of_th_SN_events;
    sp->feedback_data.to_distribute.SNII_number_of_kick_events =
        number_of_kin_SN_events;
  }
}

/**
 * @brief Return the change in temperature (in internal units) to apply to a
 * gas particle affected by SNIa feedback.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
INLINE static double eagle_SNIa_feedback_temperature_change(
    const struct spart* sp, const struct feedback_props* props) {

  /* In the EAGLE REF model, the change of temperature is constant */
  return props->SNIa_deltaT_desired;
}

/**
 * @brief Computes the fraction of the available super-novae Ia energy to
 * inject for a given event.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
INLINE static double eagle_SNIa_feedback_energy_fraction(
    const struct spart* sp, const struct feedback_props* props) {

  /* Model parameters */
  const double SNIa_f_E = props->SNIa_f_E;

  return SNIa_f_E;
}

/**
 * @brief Compute the properties of the SNIa stochastic feedback energy
 * injection.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in internal units.
 * @param dt Length of time-step in internal units.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel.
 * @param feedback_props The properties of the feedback model.
 */
INLINE static void compute_SNIa_feedback(
    struct spart* sp, const double star_age, const double dt,
    const float ngb_gas_mass, const struct feedback_props* feedback_props,
    const double dt_Gyr, const double star_age_Gyr) {

  /* Properties of the model (all in internal units) */
  const double delta_T =
      eagle_SNIa_feedback_temperature_change(sp, feedback_props);
  const double N_SNe = dtd_number_of_SNIa(
      sp, star_age_Gyr, star_age_Gyr + dt_Gyr, feedback_props);
  const double E_SNe = feedback_props->E_SNIa;
  const double f_E = eagle_SNIa_feedback_energy_fraction(sp, feedback_props);

  /* Total mass ejected at this time-step by the stellar particle sp */
  const float M_ej = sp->feedback_data.to_distribute.mass;

  /* Total un-weighted mass in the star's kernel including the ejecta mass */
  const float ngb_gas_mass_new = ngb_gas_mass + M_ej;

  /* Conversion factor from T to internal energy */
  const double conv_factor = feedback_props->temp_to_u_factor;

  /* Calculate the default heating probability */
  double prob =
      f_E * E_SNe * N_SNe / (conv_factor * delta_T * ngb_gas_mass_new);

  /* Calculate the change in internal energy of the gas particles that get
   * heated */
  double delta_u;
  if (prob <= 1.) {

    /* Normal case */
    delta_u = delta_T * conv_factor;

  } else {

    /* Special case: we need to adjust the energy irrespective of the
       desired deltaT to ensure we inject all the available energy. */

    prob = 1.;
    delta_u = f_E * E_SNe * N_SNe / ngb_gas_mass_new;
  }

  /* Store all of this in the star for delivery onto the gas */
  sp->SNIa_f_E = f_E;
  sp->feedback_data.to_distribute.SNIa_heating_probability = prob;
  sp->feedback_data.to_distribute.SNIa_delta_u = delta_u;
}

/**
 * @brief Find the bins and offset along the metallicity dimension of the
 * AGB yields table.
 *
 * @param iz_low (return) Lower index along the metallicity dimension.
 * @param iz_high (return) High index along the metallicity dimension.
 * @param dz (return) Offset between the metallicity bin and Z.
 * @param log10_Z log10 of the star metallicity (metal mass fraction).
 * @param props The properties of the feedback model.
 */
INLINE static void determine_bin_yield_AGB(int* iz_low, int* iz_high, float* dz,
                                           const float log10_Z,
                                           const struct feedback_props* props) {

  const double* AGB_Z = props->yield_AGB.metallicity;
  const int N_bins = eagle_feedback_AGB_N_metals;

  if (log10_Z > log10_min_metallicity) {

    /* Find metallicity bin which contains the star's metallicity */
    int j;
    for (j = 0; j < (N_bins - 1) && log10_Z > AGB_Z[j + 1]; j++)
      ;

    /* Store the indices */
    *iz_low = j;
    *iz_high = *iz_low + 1;

    *iz_high = min(*iz_high, N_bins - 1);

    /* Compute offset */
    if ((log10_Z >= AGB_Z[0]) && (log10_Z <= AGB_Z[N_bins - 1])) {

      *dz = log10_Z - AGB_Z[*iz_low];
    } else {
      *dz = 0.f;
    }

    /* Normalize offset */
    const float delta_Z = AGB_Z[*iz_high] - AGB_Z[*iz_low];

    if (delta_Z > 0.f)
      *dz /= delta_Z;
    else
      *dz = 0.f;

  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0.f;
  }
}

/**
 * @brief Find the bins and offset along the metallicity dimension of the
 * SNII yields table.
 *
 * @param iz_low (return) Lower index along the metallicity dimension.
 * @param iz_high (return) High index along the metallicity dimension.
 * @param dz (return) Offset between the metallicity bin and Z.
 * @param log10_Z log10 of the star metallicity (metal mass fraction).
 * @param props The properties of the feedback model.
 */
INLINE static void determine_bin_yield_SNII(
    int* iz_low, int* iz_high, float* dz, const float log10_Z,
    const struct feedback_props* props) {

  const double* SNII_Z = props->yield_SNII.metallicity;
  const int N_bins = eagle_feedback_SNII_N_metals;

  if (log10_Z > log10_min_metallicity) {

    /* Find metallicity bin which contains the star's metallicity */
    int j;
    for (j = 0; j < (N_bins - 1) && log10_Z > SNII_Z[j + 1]; j++)
      ;

    /* Store the indices */
    *iz_low = j;
    *iz_high = *iz_low + 1;

    *iz_high = min(*iz_high, N_bins - 1);

    /* Compute offset */
    if ((log10_Z >= SNII_Z[0]) && (log10_Z <= SNII_Z[N_bins - 1])) {

      *dz = log10_Z - SNII_Z[*iz_low];
    } else {
      *dz = 0.f;
    }

    /* Normalize offset */
    const float delta_Z = SNII_Z[*iz_high] - SNII_Z[*iz_low];

    if (delta_Z > 0.f)
      *dz = *dz / delta_Z;
    else
      *dz = 0.f;

  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0.f;
  }
}

/**
 * @brief Computes the number of neutron star-neutron star mergers
 * per yr for a given star particle since it was formed.
 *
 *
 * @param sp The #spart.
 * @param t Elapsed time (in Gyr).
 * @param props The properties of the stellar model.
 */
double integrate_rate_of_NSM(const struct spart* sp, const double t0,
                             const double t1,
                             const struct feedback_props* props) {

  /* The calculation is written as the integral between t0 and t1 */
  double num_NSM_per_Msun = props->NSM_per_Msun * log(t1 / t0);
  return num_NSM_per_Msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Stochastic implementation of enrichment of r-process elements
 * due to neutron star mergers (NSM).
 * To do this compute the number of NSM that occur during the timestep
 * and multiply by constants.
 *
 * The number of NSM events is drawn from a 1/t DTD.
 *
 * The Eu mass ejected per event is independent of Z.
 *
 * @param props properties of the feedback model.
 * @param sp #spart we are computing feedback from.
 * @param star_age_Gyr age of star in Gyr.
 * @param dt_Gyr timestep dt in Gyr.
 * @param ti_current Current integer time (for random numbers).
 * @param cosmo The cosmological model (for logging).
 */
INLINE static void evolve_NSM_stochastic(const struct feedback_props* props,
                                         struct spart* sp,
                                         const double star_age_Gyr,
                                         const double dt_Gyr,
                                         const integertime_t ti_current,
                                         const struct cosmology* cosmo) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dt_Gyr < 0.) error("Negative time-step length!");
  if (star_age_Gyr < 0.) error("Negative age!");
#endif

  /* First we check that the amount of time since star was formed
   is larger than NSM_t_delay_Gyr = 30Myr */
  if (star_age_Gyr < 0.03) return;

  /* Compute the number of NS merger events in timestep by
   integrating the rate of NS merger events (number per yr) */
  float num_NSM =
      integrate_rate_of_NSM(sp, star_age_Gyr, star_age_Gyr + dt_Gyr, props);

  int num_events = 0;
  if (num_NSM > 1.f) {
    num_events = floor(num_NSM);
    num_NSM -= num_events;
  }

  /* I define my probability as the left-over after we removed the integer
   * number of events */
  float prob_num = num_NSM;

  /* Draw a random number */
  const float rand =
      random_unit_interval(sp->id, ti_current, random_number_enrichment_1);

  /* Are we lucky? If so we have 1 more event */
  if (prob_num > rand) num_events++;

  if (num_events > 0) {

    /* Compute the mass produced by NSM events */
    const float delta_mass =
        num_events * props->yield_Eu_from_NSM * props->solar_mass_to_mass;

    /* compute mass of Europium */
    sp->feedback_data.to_distribute.metal_mass[chemistry_element_Eu] +=
        delta_mass;

    /* Mass to distribute */
    sp->feedback_data.to_distribute.mass_from_NSM += delta_mass;

    /* Write the event to the r-process log file */
    event_logger_r_processes_log_event(sp, cosmo, delta_mass, num_events,
                                       /*flag=*/0);
  }
}

/**
 * @brief Stochastic implementation of enrichment of r-process elements
 * due to common envelope jets supernovae (CEJSN, type of rare core-collapse
 * SN).
 *
 * To do this compute the number of CEJSN that occur during the timestep
 * and multiply by constants.
 * This assumes a uniform probability over the Z-dependent lifetime range of
 * SNII stars.
 *
 * The Eu mass ejected per event is independent of Z.
 *
 * @param log10_min_mass log10 of the minimal mass of stars dying in this step
 * in solar masses.
 * @param log10_max_mass log10 of the maximal mass of stars dying in this step
 * in solar masses.
 * @param props properties of the feedback model.
 * @param sp #spart we are computing feedback for.
 * @param Z The star's metallicity.
 * @param star_age_Gyr age of star in Gyr.
 * @param dt_Gyr timestep dt in Gyr.
 * @param ti_current Current integer time (for random numbers).
 * @param cosmo The cosmological model (for logging).
 */
INLINE static void evolve_CEJSN_stochastic(
    const float log10_min_mass, const float log10_max_mass,
    const struct feedback_props* props, struct spart* sp, const float Z,
    const double star_age_Gyr, const double dt_Gyr,
    const integertime_t ti_current, const struct cosmology* cosmo) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dt_Gyr < 0.) error("Negative time-step length!");
  if (star_age_Gyr < 0.) error("Negative age!");
#endif

  /* Abort early if the star is clearly too old or clearly too young */
  if (log10_max_mass < props->log10_SNII_min_mass_msun) return;
  if (log10_min_mass > props->log10_SNII_max_mass_msun) return;

  /* Compute the lifetime of stars with the min and max SNII masses
   * and metallicities */
  const float CEJSN_max_mass = exp10f(props->log10_SNII_max_mass_msun);
  const float CEJSN_min_mass = exp10f(props->log10_SNII_min_mass_msun);
  const float lifetime_Gyr_max_mass = lifetime_in_Gyr(CEJSN_max_mass, Z, props);
  const float lifetime_Gyr_min_mass = lifetime_in_Gyr(CEJSN_min_mass, Z, props);

  /* Compute the age of the star at the beginning and end of step */
  double t_start_Gyr = star_age_Gyr;
  double t_end_Gyr = star_age_Gyr + dt_Gyr;

  /* Do we need to correct because the age window overalps with the
     minimal or maximal age? */
  if (t_end_Gyr > lifetime_Gyr_min_mass) t_end_Gyr = lifetime_Gyr_min_mass;
  if (t_start_Gyr < lifetime_Gyr_max_mass) t_start_Gyr = lifetime_Gyr_max_mass;

  /* The length of the step used to determine probabilities */
  const double step_length_Gyr = t_end_Gyr - t_start_Gyr;

  /* Range over which the CEJSN are sampled */
  const float delta_lifetime_Gyr =
      lifetime_Gyr_min_mass - lifetime_Gyr_max_mass;

  /* Compute probability based on step_length / lifetime */

  /* Number of CEJSN events in timestep */
  float num_CEJSN = props->CEJSN_per_Msun * sp->mass_init *
                    props->mass_to_solar_mass *
                    (step_length_Gyr / delta_lifetime_Gyr);

  int num_events = 0;
  if (num_CEJSN > 1) {
    num_events = floor(num_CEJSN);
    num_CEJSN -= num_events;
  }

  /* I define my probability as the left-over after we removed the integer
   * number of events */
  const float prob_num = num_CEJSN;

  /* Draw a random number */
  const float rand =
      random_unit_interval(sp->id, ti_current, random_number_enrichment_2);

  /* Are we lucky? If so we have 1 more event */
  if (prob_num > rand) num_events++;

  if (num_events > 0) {

    /* Compute the mass produced by CEJSN */
    const float delta_mass =
        num_events * props->yield_Eu_from_CEJSN * props->solar_mass_to_mass;

    /* compute mass of Europium */
    sp->feedback_data.to_distribute.metal_mass[chemistry_element_Eu] +=
        delta_mass;

    /* Mass to distribute */
    sp->feedback_data.to_distribute.mass_from_CEJSN += delta_mass;

    /* Write the event to the r-process log file */
    event_logger_r_processes_log_event(sp, cosmo, delta_mass, num_events,
                                       /*flag=*/1);
  }
}

/**
 * @brief Stochastic implementation of enrichment of r-process elements
 * due to collapsars (rare type core-collapse SN).
 *
 * To do this compute the number of collapsars that occur during the timestep
 * and multiply by constants.
 * This assumes a uniform probability over the Z-dependent lifetime range of
 * collapsars.
 *
 * The Eu mass ejected per event is independent of Z.
 *
 * @param log10_min_mass log10 of the minimal mass of stars dying in this step
 * in solar masses.
 * @param log10_max_mass log10 of the maximal mass of stars dying in this step
 * in solar masses.
 * @param props properties of the feedback model.
 * @param sp #spart we are computing feedback for.
 * @param Z The star's metallicity.
 * @param star_age_Gyr age of star in Gyr.
 * @param dt_Gyr timestep dt in Gyr.
 * @param ti_current Current integer time (for random numbers).
 * @param cosmo The cosmological model (for logging).
 */
INLINE static void evolve_collapsar_stochastic(
    const float log10_min_mass, const float log10_max_mass,
    const struct feedback_props* props, struct spart* sp, const float Z,
    const double star_age_Gyr, const double dt_Gyr,
    const integertime_t ti_current, const struct cosmology* cosmo) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dt_Gyr < 0.) error("Negative time-step length!");
  if (star_age_Gyr < 0.) error("Negative age!");
#endif

  /* Abort early if the star is clearly too old or clearly too young */
  if (log10_max_mass < props->log10_collapsar_min_mass_msun) return;
  if (log10_min_mass > props->log10_collapsar_max_mass_msun) return;

  /* Compute the lifetime of stars with the min and max collapsar masses
   * and metallicities */
  const float collapsar_max_mass = exp10f(props->log10_collapsar_max_mass_msun);
  const float collapsar_min_mass = exp10f(props->log10_collapsar_min_mass_msun);
  const float lifetime_Gyr_max_mass =
      lifetime_in_Gyr(collapsar_max_mass, Z, props);
  const float lifetime_Gyr_min_mass =
      lifetime_in_Gyr(collapsar_min_mass, Z, props);

  /* Compute the age of the star at the beginning and end of step */
  double t_start_Gyr = star_age_Gyr;
  double t_end_Gyr = star_age_Gyr + dt_Gyr;

  /* Do we need to correct because the age window overalps with the
     minimal or maximal age? */
  if (t_end_Gyr > lifetime_Gyr_min_mass) t_end_Gyr = lifetime_Gyr_min_mass;
  if (t_start_Gyr < lifetime_Gyr_max_mass) t_start_Gyr = lifetime_Gyr_max_mass;

  /* The length of the step used to determine probabilities */
  const double step_length_Gyr = t_end_Gyr - t_start_Gyr;

  /* Range over which the collapsars are sampled */
  const float delta_lifetime_Gyr =
      lifetime_Gyr_min_mass - lifetime_Gyr_max_mass;

  /* Compute probability based on step_length / lifetime */

  /* Number of collapsar events in timestep */
  float num_collapsar = props->collapsar_per_Msun * sp->mass_init *
                        props->mass_to_solar_mass *
                        (step_length_Gyr / delta_lifetime_Gyr);

  int num_events = 0;
  if (num_collapsar > 1.f) {
    num_events = floor(num_collapsar);
    num_collapsar -= num_events;
  }

  /* I define my probability as the left-over after we removed the integer
   * number of events */
  const float prob_num = num_collapsar;

  /* Draw a random number */
  const float rand =
      random_unit_interval(sp->id, ti_current, random_number_enrichment_3);

  /* Are we lucky? If so we have 1 more event */
  if (prob_num > rand) num_events++;

  if (num_events > 0) {

    /* Compute the mass produced by collapsar */
    const float delta_mass =
        num_events * props->yield_Eu_from_collapsar * props->solar_mass_to_mass;

    /* compute mass of Europium */
    sp->feedback_data.to_distribute.metal_mass[chemistry_element_Eu] +=
        delta_mass;

    /* Mass to distribute */
    sp->feedback_data.to_distribute.mass_from_collapsar += delta_mass;

    /* Write the event to the r-process log file */
    event_logger_r_processes_log_event(sp, cosmo, delta_mass, num_events,
                                       /*flag=*/2);
  }
}

/**
 * @brief compute enrichment and feedback due to SNIa. To do this compute the
 * number of SNIa that occur during the timestep, multiply by constants read
 * from tables.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param props properties of the feedback model
 * @param sp #spart we are computing feedback from
 * @param star_age_Gyr age of star in Gyr
 * @param dt_Gyr timestep dt in Gyr
 */
INLINE static void evolve_SNIa(const float log10_min_mass,
                               const float log10_max_mass,
                               const struct feedback_props* props,
                               struct spart* sp, double star_age_Gyr,
                               double dt_Gyr) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dt_Gyr < 0.) error("Negative time-step length!");
  if (star_age_Gyr < 0.) error("Negative age!");
#endif

  /* Compute the number of SNIa */
  const float num_SNIa =
      dtd_number_of_SNIa(sp, star_age_Gyr, star_age_Gyr + dt_Gyr, props);

  /* compute mass of each metal */
  for (int i = 0; i < enrichment_of_N_elements_from_yield_tables; i++) {
    sp->feedback_data.to_distribute.metal_mass[i] +=
        num_SNIa * props->yield_SNIa_IMF_resampled[i] *
        props->solar_mass_to_mass;
  }

  /* Update the metallicity of the material released */
  sp->feedback_data.to_distribute.metal_mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Update the metal mass produced */
  sp->feedback_data.to_distribute.total_metal_mass +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Compute the mass produced by SNIa
   * Note: SNIa do not inject H or He so the mass injected is the same
   * as the metal mass injected. */
  sp->feedback_data.to_distribute.mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Compute the iron mass produced */
  sp->feedback_data.to_distribute.Fe_mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_IMF_resampled[chemistry_element_Fe] *
      props->solar_mass_to_mass;
}

/**
 * @brief compute enrichment and feedback due to SNII. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * Note for Matthieu: This function is poorly written and needs improving.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stellar_yields array to store calculated yields for passing to
 * integrate_imf
 * @param props properties of the feedback model.
 * @param sp spart we are computing feedback from
 */
INLINE static void evolve_SNII(float log10_min_mass, float log10_max_mass,
                               float* stellar_yields,
                               const struct feedback_props* props,
			       const struct dustevo_props* dp,
                               struct spart* sp) {

  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* Metallicity (metal mass fraction) at birth time of the star */
  const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* If mass at beginning of step is less than tabulated lower bound for IMF,
   * limit it.*/
  if (log10_min_mass < props->log10_SNII_min_mass_msun)
    log10_min_mass = props->log10_SNII_min_mass_msun;

  /* If mass at end of step is greater than tabulated upper bound for IMF, limit
   * it.*/
  if (log10_max_mass > props->log10_SNII_max_mass_msun)
    log10_max_mass = props->log10_SNII_max_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, props);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low = 0, iz_high = 0, low_index_3d, high_index_3d, low_index_2d,
      high_index_2d;
  float dz = 0.;
  determine_bin_yield_SNII(&iz_low, &iz_high, &dz, log10(Z), props);

  /* compute metals produced */
  float metal_mass_released[enrichment_of_N_elements_from_yield_tables],
      metal_mass_released_total;
  for (int elem = 0; elem < enrichment_of_N_elements_from_yield_tables;
       elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          enrichment_of_N_elements_from_yield_tables,
          eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          enrichment_of_N_elements_from_yield_tables,
          eagle_feedback_N_imf_bins);
      low_index_2d = row_major_index_2d(iz_low, mass_bin_index,
                                        eagle_feedback_SNII_N_metals,
                                        eagle_feedback_N_imf_bins);
      high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                         eagle_feedback_SNII_N_metals,
                                         eagle_feedback_N_imf_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) *
              (props->yield_SNII.yield_IMF_resampled[low_index_3d] +
               sp->chemistry_data.metal_mass_fraction[elem] *
                   props->yield_SNII.ejecta_IMF_resampled[low_index_2d]) +
          dz * (props->yield_SNII.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    props->yield_SNII.ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_SNII_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                       eagle_feedback_SNII_N_metals,
                                       eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * (props->yield_SNII.total_metals_IMF_resampled[low_index_2d] +
                    sp->chemistry_data.metal_mass_fraction_total *
                        props->yield_SNII.ejecta_IMF_resampled[low_index_2d]) +
        dz * (props->yield_SNII.total_metals_IMF_resampled[high_index_2d] +
              sp->chemistry_data.metal_mass_fraction_total *
                  props->yield_SNII.ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* compute dust produced */
  float dust_mass_released[grain_species_count];

  for (int grain = 0; grain < grain_species_count;
       grain++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, grain, mass_bin_index, eagle_feedback_SNII_N_metals,
	  grain_species_count,
          eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, grain, mass_bin_index, eagle_feedback_SNII_N_metals,
          grain_species_count,
          eagle_feedback_N_imf_bins);

      stellar_yields[mass_bin_index] =
          (1 - dz) * (dp->dyield_SNII.yield_IMF_resampled[low_index_3d]) +
          dz * (dp->dyield_SNII.yield_IMF_resampled[high_index_3d]);
    }

    dust_mass_released[grain] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }


  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < enrichment_of_N_elements_from_yield_tables; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  for (int grain = 0; grain < grain_species_count; grain++)
    dust_mass_released[grain] = max(dust_mass_released[grain], 0.f);  

  /* compute the total mass ejected from the star*/
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_SNII_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                       eagle_feedback_SNII_N_metals,
                                       eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * props->yield_SNII.ejecta_IMF_resampled[low_index_2d] +
        dz * props->yield_SNII.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* compute the total mass released */
  mass_released = metal_mass_released_total +
                  metal_mass_released[chemistry_element_H] +
                  metal_mass_released[chemistry_element_He];

  /* normalize the yields */
  if (mass_released > 0) {
    /* Set normalisation factor. Note additional multiplication by the star
     * initial mass as tables are per initial mass */
    const float norm_factor = sp->mass_init * mass_ejected / mass_released;

    for (int i = 0; i < enrichment_of_N_elements_from_yield_tables; i++) {
      sp->feedback_data.to_distribute.metal_mass[i] +=
          metal_mass_released[i] * norm_factor;
    }
    for (int i = 0; i < enrichment_of_N_elements_from_yield_tables; i++) {
      sp->feedback_data.to_distribute.mass_from_SNII +=
          sp->feedback_data.to_distribute.metal_mass[i];
    }
    sp->feedback_data.to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->feedback_data.to_distribute.metal_mass_from_SNII +=
        metal_mass_released_total * norm_factor;

    for (int grain = 0; grain < grain_species_count; grain++) {
      sp->feedback_data.to_distribute.dust_mass[grain] +=
          dust_mass_released[grain] * norm_factor;
    }

  } else {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
}

/**
 * @brief compute enrichment and feedback due to AGB. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * Note for Matthieu: This function is poorly written and needs improving.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stellar_yields array to store calculated yields for passing to
 * integrate_imf
 * @param props Properties of the feedback model.
 * @param sp spart we are computing feedback for.
 */
INLINE static void evolve_AGB(const float log10_min_mass, float log10_max_mass,
                              float* stellar_yields,
                              const struct feedback_props* props,
			      const struct dustevo_props* dp,
                              struct spart* sp) {

  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* Metallicity (metal mass fraction) at birth time of the star */
  const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* If mass at end of step is greater than tabulated lower bound for IMF, limit
   * it.*/
  if (log10_max_mass > props->log10_SNII_min_mass_msun)
    log10_max_mass = props->log10_SNII_min_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, props);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low = 0, iz_high = 0, low_index_3d, high_index_3d, low_index_2d,
      high_index_2d;
  float dz = 0.f;
  determine_bin_yield_AGB(&iz_low, &iz_high, &dz, log10(Z), props);

  /* compute metals produced */
  float metal_mass_released[enrichment_of_N_elements_from_yield_tables],
      metal_mass_released_total;
  for (int elem = 0; elem < enrichment_of_N_elements_from_yield_tables;
       elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          enrichment_of_N_elements_from_yield_tables,
          eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          enrichment_of_N_elements_from_yield_tables,
          eagle_feedback_N_imf_bins);
      low_index_2d = row_major_index_2d(iz_low, mass_bin_index,
                                        eagle_feedback_AGB_N_metals,
                                        eagle_feedback_N_imf_bins);
      high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                         eagle_feedback_AGB_N_metals,
                                         eagle_feedback_N_imf_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) * (props->yield_AGB.yield_IMF_resampled[low_index_3d] +
                      sp->chemistry_data.metal_mass_fraction[elem] *
                          props->yield_AGB.ejecta_IMF_resampled[low_index_2d]) +
          dz * (props->yield_AGB.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    props->yield_AGB.ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d =
        row_major_index_2d(iz_high, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * (props->yield_AGB.total_metals_IMF_resampled[low_index_2d] +
                    sp->chemistry_data.metal_mass_fraction_total *
                        props->yield_AGB.ejecta_IMF_resampled[low_index_2d]) +
        dz * (props->yield_AGB.total_metals_IMF_resampled[high_index_2d] +
              sp->chemistry_data.metal_mass_fraction_total *
                  props->yield_AGB.ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);


  /* compute dust produced */
  float dust_mass_released[grain_species_count];

  for (int grain = 0; grain < grain_species_count;
       grain++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, grain, mass_bin_index, eagle_feedback_AGB_N_metals,
	  grain_species_count,
          eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, grain, mass_bin_index, eagle_feedback_AGB_N_metals,
          grain_species_count,
          eagle_feedback_N_imf_bins);

      stellar_yields[mass_bin_index] =
          (1 - dz) * (dp->dyield_AGB.yield_IMF_resampled[low_index_3d]) +
          dz * (dp->dyield_AGB.yield_IMF_resampled[high_index_3d]);
    }

    dust_mass_released[grain] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < enrichment_of_N_elements_from_yield_tables; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  for (int grain = 0; grain < grain_species_count; grain++)
    dust_mass_released[grain] = max(dust_mass_released[grain], 0.f);

  /* compute the total mass ejected from the star */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d =
        row_major_index_2d(iz_high, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * props->yield_AGB.ejecta_IMF_resampled[low_index_2d] +
        dz * props->yield_AGB.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* compute the total mass released */
  mass_released = metal_mass_released_total +
                  metal_mass_released[chemistry_element_H] +
                  metal_mass_released[chemistry_element_He];

  /* normalize the yields */
  if (mass_released > 0) {

    /* Set normalisation factor. Note additional multiplication by the stellar
     * initial mass as tables are per initial mass */
    const float norm_factor = sp->mass_init * mass_ejected / mass_released;

    for (int i = 0; i < enrichment_of_N_elements_from_yield_tables; i++) {
      sp->feedback_data.to_distribute.metal_mass[i] +=
          metal_mass_released[i] * norm_factor;
      sp->feedback_data.to_distribute.mass_from_AGB +=
          metal_mass_released[i] * norm_factor;
    }
    sp->feedback_data.to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->feedback_data.to_distribute.metal_mass_from_AGB +=
        metal_mass_released_total * norm_factor;

    for (int grain = 0; grain < grain_species_count; grain++) {
      sp->feedback_data.to_distribute.dust_mass[grain] +=
          dust_mass_released[grain] * norm_factor;
    }

  } else {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
}

/**
 * @brief Gets interpolated cumulative stellar wind momentum input from table
 *
 * @param props feedback_props structure for getting model parameters for
 * coefficients
 * @param t_Myr stellar age in Myr
 * @param logZ log10 of (stellar) metal mass fraction Z
 */
double get_cumulative_stellarwind_momentum(const struct feedback_props* fp,
                                           float t_Myr, float log10_Z) {
  float d_age, d_met;
  int met_index, age_index;

  if (t_Myr < fp->HII_agebins[0]) return 0.;

  /* Get index alongside the metallicity dimension */
  get_index_1d(fp->HII_log10_Zbins, fp->HII_nr_metbins, log10_Z, &met_index,
               &d_met);

  /* Get index alongside the age dimension */
  get_index_1d(fp->HII_agebins, fp->HII_nr_agebins, t_Myr, &age_index, &d_age);

  /* Interpolate! */
  const float log10_Pcum_loc =
      interpolation_2d_flat(fp->SW_log10_Pcum, met_index, age_index, d_met,
                            d_age, fp->HII_nr_metbins, fp->HII_nr_agebins);

  return exp10(log10_Pcum_loc);
}

/**
 * @brief Gets interpolated cumulative ionizing photons from table
 * @param props feedback_props structure for getting model parameters for
 * coefficients
 * @param t_Myr stellar age in Myr
 * @param logZ log10 of (stellar) metal mass fraction Z
 */
double get_cumulative_ionizing_photons(const struct feedback_props* fp,
                                       float t_Myr, float log10_Z) {
  float d_age, d_met;
  int met_index, age_index;

  if (t_Myr < fp->HII_agebins[0]) return 0.;

  /* Get index alongside the metallicity dimension */
  get_index_1d(fp->HII_log10_Zbins, fp->HII_nr_metbins, log10_Z, &met_index,
               &d_met);

  /* Get index alongside the age dimension */
  get_index_1d(fp->HII_agebins, fp->HII_nr_agebins, t_Myr, &age_index, &d_age);

  /* Interpolate! */
  const float log10_Qcum_loc =
      interpolation_2d_flat(fp->HII_log10_Qcum, met_index, age_index, d_met,
                            d_age, fp->HII_nr_metbins, fp->HII_nr_agebins);

  return exp10(log10_Qcum_loc);
}

/**
 * @brief Calculates the average stellar wind momentum input between t1 and t2
 * for an initial metallicity of Z
 *
 * @param props feedback_props structure for getting model parameters
 * @param t1 initial time in Myr
 * @param t2 final time in Myr
 * @param Z metal mass fraction
 * @param Pbar average momentum input of a star with metallicity Z over this
 * period of time (t1 - t2)
 */
double compute_average_stellarwind_momentum(const struct feedback_props* fp,
                                            float t1, float t2, float Z) {

  const float log10_Z = log10(Z);
  const double P_t1 = get_cumulative_stellarwind_momentum(fp, t1, log10_Z);
  const double P_t2 = get_cumulative_stellarwind_momentum(fp, t2, log10_Z);

  return ((P_t2 - P_t1) / (t2 - t1)) * fp->sec_to_Myr;
}

/**
 * @brief Calculates the average ionizing luminosity between t1 and t2 for an
 * initial metallicity of Z
 *
 * @param props feedback_props structure for getting model parameters
 * @param t1 initial time in Myr
 * @param t2 final time in Myr
 * @param Z metal mass fraction
 * @param Qbar photoionizing luminosity of a star with metallicity Z over this
 * period of time (t1 - t2)
 */
double compute_average_photoionizing_luminosity(const struct feedback_props* fp,
                                                float t1, float t2, float Z) {

  /* No luminosity at t=0 */
  if (t2 <= 0.f) return 0.f;

  const float log10_Z = log10(Z + FLT_MIN);
  const double Q_t1 = get_cumulative_ionizing_photons(fp, t1, log10_Z);
  const double Q_t2 = get_cumulative_ionizing_photons(fp, t2, log10_Z);

  return ((Q_t2 - Q_t1) / (t2 - t1)) * fp->sec_to_Myr;
}

/**
 * @brief Calculates the amount of momentum available for this star
 * from Starburst 99. Fitting function taken from Agertz et al. (2013)
 *
 * @param sp spart that we're evolving
 * @param us unit_system structure for unit conversion
 * @param props feedback_props structure for getting model parameters
 * @param star_age_Gyr Age of star in Gyr
 * @param dt current timestep in internal units
 * @param ngb_gas_mass mass within the sph kernel
 */
INLINE static void compute_stellar_momentum(struct spart* sp,
                                            const struct unit_system* us,
                                            const struct feedback_props* props,
                                            const double star_age_Gyr,
                                            const double dt,
                                            const float ngb_gas_mass) {

  const double tw = props->SW_max_age_Myr;      /* Myr */
  const double delta_v_km_p_s = props->delta_v; /* km s^-1 */

  /* delta_v in code units */
  double delta_v =
      delta_v_km_p_s /
      (units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) * 1.0e-5);

  /* Unit conversion constant */
  const double Myr_in_s = 1.0e6 * 365 * 24 * 3600.;

  /* Convert the times to the units used by the model */
  const double star_age_Myr = star_age_Gyr * 1e3;
  double dt_cgs = dt * us->UnitTime_in_cgs;
  const double dt_Myr = dt_cgs / Myr_in_s;

  /* Prevent star particle from injecting momentum for longer than tw */
  float dt_new = dt;
  if (star_age_Myr + dt_Myr > tw) {
    dt_new = (tw - star_age_Myr) * Myr_in_s / us->UnitTime_in_cgs;
    dt_cgs = dt_new * us->UnitTime_in_cgs;
  }

  /* Star too old to do any momentum injection or time-step is 0 */
  if (star_age_Myr > tw || dt == 0.) return;

  /* Mass within the SPH kernel in grams */
  const double ngb_gas_mass_in_g = ngb_gas_mass * us->UnitMass_in_cgs;

  /* Star metallicity (metal mass fraction) at birth */
  double Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* Bring the metallicity in the range covered by the model */
  Z = max(Z, props->Zmin_early_fb);
  Z = min(Z, props->Zmax_early_fb);

  /* Get the average momentum input from stellar winds during this timestep
   * from the BPASS tables */
  const float t1_Myr = (float)star_age_Gyr * 1.e3;
  const float t2_Myr = t1_Myr + (float)dt_Myr;
  const double Pbar =
      compute_average_stellarwind_momentum(props, t1_Myr, t2_Myr, Z);
  const double P_cgs = Pbar * sp->mass_init * us->UnitMass_in_cgs * dt_cgs;

  /* Velocity kick */
  const double delta_v_cgs = delta_v_km_p_s * 1e5;

  /* Get the available momentum in code units at the given dt and store it */
  const double momentum = P_cgs / props->Momentum_to_cgs;

  /* Now compute the probability of kicking particle with given delta_v
   * in the current timestep.
   * Note that this could be prop > 1 if there are no enough particles in the
   * kernel to distribute the amount of momentum available in
   * the timestep, but this is OK because we then adjust delta_v */

  double prob = 0.;

  /* We want the code to decide the velocity kick for us */
  if (delta_v < 0.) {

    /* We kick all particles in the kernel */
    prob = 1.;

    /* Kick velocity (in code units) needed so that we can inject
     * all the momentum inside the kernel */
    delta_v = momentum / ngb_gas_mass;

    /* User decided velocity kick */
  } else {

    /* The probability of kicking a particle is given
     * by the kick velocity chosen in the parameter file
     * and is normalised by the mass available in the kernel */
    prob = (P_cgs / delta_v_cgs) / ngb_gas_mass_in_g;

    /* Mass inside the kernel too small makes prob > 1 */
    if (prob > 1.) {
      message("Not enough parts in the kernel to distribute momentum...");

      /* Correct the kick (in code units) to be consistent with the mass within
       * the kernel and amount of momentum available */
      delta_v = momentum / ngb_gas_mass;

      message("Wind speed set to delta_v = %e km/s",
              delta_v * units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) *
                  1.0e-5);
    }
  }

  /* Store values for perfoming the actual kick later on */
  sp->feedback_data.to_distribute.momentum_probability = prob;
  sp->feedback_data.to_distribute.momentum_delta_v = delta_v;
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 * @param time_beg_of_step time at the beginning of timestep
 * @param ti_begin The integer time at the beginning of the step (for random
 * numbers).
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo, 
			       const struct dustevo_props* dustevo_props,
			       struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const double time_beg_of_step,
                               const integertime_t ti_begin) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");
#endif

  /* Allocate temporary array for calculating imf weights */
  float stellar_yields[eagle_feedback_N_imf_bins];

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_in_cgs = 1e9 * 365.25 * 24. * 3600.;
  const double Myr_in_cgs = 1e6 * 365.25 * 24. * 3600.;
  const double time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const double conversion_factor = time_to_cgs / Gyr_in_cgs;
  const double dt_Gyr = dt * conversion_factor;
  const double dt_Myr = dt * conversion_factor * 1e3;
  const double star_age_Gyr = age * conversion_factor;
  const double star_age_Myr = age * conversion_factor * 1e3;

  /* Get the total metallicity (metal mass fraction) at birth time and impose a
   * minimum */
  const float Z = max(chemistry_get_total_metal_mass_fraction_for_feedback(sp),
                      exp10(log10_min_metallicity));

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;
  const int ngb_Number = sp->feedback_data.to_collect.ngb_N;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.enrichment_weight_inv < 0.)
    error("Negative inverse weight ! %e id=%lld",
          sp->feedback_data.to_collect.enrichment_weight_inv, sp->id);
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight = 1.f / enrichment_weight_inv;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_weight < 0.)
    error("Negative weight after reset!");
#endif

  /* Compute amount of momentum available for this stars, given its mass and
     age, only if needed */
  if (feedback_props->with_StellarWinds) {
    compute_stellar_momentum(sp, us, feedback_props, star_age_Gyr, dt,
                             ngb_gas_mass);
  }

  sp->star_timestep = dt;

  /* Compute ionizing photons for HII regions only if needed*/
  if (feedback_props->with_HIIRegions &&
      star_age_Myr <= feedback_props->HIIregion_max_age_Myr) {

    /* only rebuild every HIIregion_dtMyr and at the first timestep the star was
     * formed*/
    if ((((sp->HIIregion_last_rebuild + feedback_props->HIIregion_dt_Myr) >=
          star_age_Myr) &&
         ((sp->HIIregion_last_rebuild + feedback_props->HIIregion_dt_Myr) <
          star_age_Myr + dt_Myr)) ||
        (sp->HIIregion_last_rebuild < 0.)) {

      /* log when this HII region was (re)built */
      double old_star_age_Myr;
      if (sp->HIIregion_last_rebuild >= 0.) {
        old_star_age_Myr = sp->HIIregion_last_rebuild;
      } else {
        old_star_age_Myr = 0.;
      }

      sp->HIIregion_last_rebuild = star_age_Myr;

      const double rho_birth = (double)sp->sf_data.birth_density;
      const double n_birth = rho_birth * feedback_props->rho_to_n_cgs;
      const double alpha_B = (double)feedback_props->alpha_caseb_recomb;
      const double t_half = age * time_to_cgs + 0.5 * dt_Myr * Myr_in_cgs;

      double Qbar;
      float t1_Myr = (float)old_star_age_Myr;
      float t2_Myr = (float)star_age_Myr;

      Qbar = compute_average_photoionizing_luminosity(feedback_props, t1_Myr,
                                                      t2_Myr, Z);

      /* Time-dependent solution of the Stromgren sphere */
      /* R(t) = R_S (1 - e^(-t/t_rec) )^(1/3) */
      /*    with the recombination timescale t_rec = 1/ (n*alpha_B) */
      /*    and R_S the Stromgren radius R_S = (3/(4 pi alpha_B) * Q(t) /
       * n^2)^(1/3) */
      /* masses in system units */
      /* [n_birth] = cm-3 */
      /* [alpha_B] = cgs */
      /* [t_half]  = s */
      /* [Qbar]    = average number of ionizing photons per second per g stellar
       * mass */
      sp->HIIregion_mass_to_ionize =
          (float)(0.84 * (double)sp->mass_init *
                  (1. - exp(-alpha_B * n_birth * t_half)) * (10. / n_birth) *
                  (Qbar / 1.e12));

#ifdef SWIFT_DEBUG_CHECKS
      if (sp->HIIregion_mass_to_ionize > 1.e10 ||
          sp->HIIregion_mass_to_ionize < 0.) {
        message("sp->mass_init = %.4e", sp->mass_init);
        message("alpha_B = %.4e", alpha_B);
        message("n_birth = %.4e", n_birth);
        message("age = %.4e", age);
        message("time_to_cgs = %.4e", time_to_cgs);
        message("Qbar = %.4e", Qbar);
        message("time term = %.4e", (1. - exp(-alpha_B * n_birth * t_half)));
        message("sp->HIIregion_mass_to_ionize = %.4e",
                sp->HIIregion_mass_to_ionize);

        error("Weird values for HII mass. Stopping.");
      }
#endif

      sp->HIIregion_mass_in_kernel = ngb_gas_mass;
      sp->feedback_data.to_distribute.HIIregion_probability =
          sp->HIIregion_mass_to_ionize / ngb_gas_mass;

      /* convert dtMyr to dt (SU) */
      const float HIIregion_dt =
          feedback_props->HIIregion_dt_Myr * Myr_in_cgs / time_to_cgs;
      sp->feedback_data.to_distribute.HIIregion_endtime =
          time_beg_of_step + HIIregion_dt;
      sp->feedback_data.to_distribute.HIIregion_starid = sp->id;

      /* Energy to distribute */
      sp->feedback_data.to_distribute.HII_u = feedback_props->HII_u;
    } else {
      sp->feedback_data.to_distribute.HIIregion_probability = -1.;
      sp->feedback_data.to_distribute.HIIregion_endtime = -1.;
      sp->feedback_data.to_distribute.HIIregion_starid = -1;
      sp->feedback_data.to_distribute.HII_u = 0.f;
    }

  } else if (feedback_props->HIIregion_max_age_Myr > 0.) {
    sp->HIIregion_last_rebuild = -1.;
    sp->feedback_data.to_distribute.HIIregion_probability = -1.;
    sp->feedback_data.to_distribute.HIIregion_endtime = -1.;
    sp->feedback_data.to_distribute.HIIregion_starid = -1;
    sp->feedback_data.to_distribute.HII_u = 0.f;
    sp->HIIregion_mass_to_ionize = 0.f;
    sp->HIIregion_mass_in_kernel = -1.f;
  }

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const float max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const float min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment and SN feedback. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

#ifdef SWIFT_DEBUG_CHECK
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Life is better in log */
  const float log10_max_dying_mass_Msun = log10f(max_dying_mass_Msun);
  const float log10_min_dying_mass_Msun = log10f(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  four channels SNIa, SNII, AGB and r-processes */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
                feedback_props, sp, star_age_Gyr, dt_Gyr);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
                stellar_yields, feedback_props, dustevo_props, sp);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
               stellar_yields, feedback_props, dustevo_props, sp);
  }
  if (feedback_props->with_r_process_enrichment) {
    evolve_NSM_stochastic(feedback_props, sp, star_age_Gyr, dt_Gyr, ti_begin,
                          cosmo);

    evolve_CEJSN_stochastic(log10_min_dying_mass_Msun,
                            log10_max_dying_mass_Msun, feedback_props, sp, Z,
                            star_age_Gyr, dt_Gyr, ti_begin, cosmo);

    evolve_collapsar_stochastic(log10_min_dying_mass_Msun,
                                log10_max_dying_mass_Msun, feedback_props, sp,
                                Z, star_age_Gyr, dt_Gyr, ti_begin, cosmo);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He + metals) */
  sp->feedback_data.to_distribute.mass =
      sp->feedback_data.to_distribute.total_metal_mass +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass *
      feedback_props->AGB_ejecta_specific_kinetic_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  /* Compute properties of the stochastic SNII feedback model. */
  if (feedback_props->with_SNII_feedback) {
    compute_SNII_feedback(sp, age, dt, ngb_gas_mass, ngb_Number, feedback_props,
                          min_dying_mass_Msun, max_dying_mass_Msun, ti_begin);
  }

  /* Compute properties of the stochastic SNIa feedback model. */
  if (feedback_props->with_SNIa_feedback) {
    compute_SNIa_feedback(sp, age, dt, ngb_gas_mass, feedback_props, dt_Gyr,
                          star_age_Gyr);
  }

  /* Star age in Myr to store in case an SNII event occurs */
  sp->feedback_data.to_distribute.SNII_star_age_Myr = (float)star_age_Myr;

  /* Modify the HII-region probability to include the ejecta mass */
  if (sp->feedback_data.to_distribute.HIIregion_probability != -1.) {
    sp->feedback_data.to_distribute.HIIregion_probability /=
        (1.f + sp->feedback_data.to_distribute.mass / ngb_gas_mass);
  }

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props* fp,
                         const struct phys_const* phys_const,
                         const struct unit_system* us,
                         struct swift_params* params,
                         const struct hydro_props* hydro_props,
                         const struct cosmology* cosmo/*, struct dustevo_props* dp*/) {

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNII_feedback =
      parser_get_param_int(params, "COLIBREFeedback:use_SNII_feedback");

  fp->with_SNIa_feedback =
      parser_get_param_int(params, "COLIBREFeedback:use_SNIa_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "COLIBREFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "COLIBREFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "COLIBREFeedback:use_SNIa_enrichment");

  fp->with_r_process_enrichment =
      parser_get_param_int(params, "COLIBREFeedback:with_r_process_enrichment");

  fp->with_HIIRegions =
      parser_get_param_int(params, "COLIBREFeedback:with_HIIRegions");

  fp->with_StellarWinds =
      parser_get_param_int(params, "COLIBREFeedback:with_StellarWinds");

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "COLIBREFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "COLIBREFeedback:IMF_min_mass_Msun");

  /* Check that it makes sense. */
  if (fp->imf_max_mass_msun < fp->imf_min_mass_msun) {
    error("Can't have the max IMF mass smaller than the min IMF mass!");
  }

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);

  /* Properties of the SNII energy feedback model ------------------------- */

  /* Set the delay time before SNII occur */
  const double Gyr_in_cgs = 1.0e9 * 365.25 * 24. * 3600.;
  fp->SNII_wind_delay =
      parser_get_param_double(params, "COLIBREFeedback:SNII_wind_delay_Gyr") *
      Gyr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Read the temperature change to use in stochastic heating */
  fp->SNII_deltaT_desired =
      parser_get_param_float(params, "COLIBREFeedback:SNII_delta_T_K");
  fp->SNII_deltaT_desired /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Energy released by supernova type II */
  fp->E_SNII_cgs =
      parser_get_param_double(params, "COLIBREFeedback:SNII_energy_erg");
  fp->E_SNII =
      fp->E_SNII_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Fraction of SNII energy injected in kinetic form */
  fp->SNII_f_kinetic =
      parser_get_param_float(params, "COLIBREFeedback:SNII_f_kinetic");

  /* Kick velocity used by supernova type II */
  fp->SNII_delta_v =
      parser_get_param_float(params, "COLIBREFeedback:SNII_delta_v_km_p_s");
  fp->SNII_delta_v *= 1e5;
  fp->SNII_delta_v /= units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Stellar mass limits for SNII feedback */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "COLIBREFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "COLIBREFeedback:SNII_max_mass_Msun");

  /* Check that it makes sense. */
  if (SNII_max_mass_msun < SNII_min_mass_msun) {
    error("Can't have the max SNII mass smaller than the min SNII mass!");
  }

  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);

  /* Properties of the energy fraction model */
  fp->f_E_min = parser_get_param_double(
      params, "COLIBREFeedback:SNII_energy_fraction_min");
  fp->f_E_max = parser_get_param_double(
      params, "COLIBREFeedback:SNII_energy_fraction_max");
  fp->Z_0 = parser_get_param_double(params,
                                    "COLIBREFeedback:SNII_energy_fraction_Z_0");
  fp->n_0_cgs = parser_get_param_double(
      params, "COLIBREFeedback:SNII_energy_fraction_n_0_H_p_cm3");
  fp->n_n = parser_get_param_double(params,
                                    "COLIBREFeedback:SNII_energy_fraction_n_n");
  fp->n_Z = parser_get_param_double(params,
                                    "COLIBREFeedback:SNII_energy_fraction_n_Z");

  /* Properties of the stochastic SNIa model */
  fp->SNIa_deltaT_desired =
      parser_get_param_double(params, "COLIBREFeedback:SNIa_delta_T_K");

  fp->SNIa_f_E =
      parser_get_param_double(params, "COLIBREFeedback:SNIa_energy_fraction");

  /* Check that it makes sense. */
  if (fp->f_E_max < fp->f_E_min) {
    error("Can't have the maximal energy fraction smaller than the minimal!");
  }

  /* Properties of the SNII enrichment model -------------------------------- */

  /* Set factors for each element adjusting SNII yield */
  for (int elem = 0; elem < enrichment_of_N_elements_from_yield_tables;
       ++elem) {
    char buffer[50];
    sprintf(buffer, "COLIBREFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Load the SNIa model */
  dtd_init(fp, phys_const, us, params);

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "COLIBREFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Properties of the r-process enrichment model
   * -------------------------------- */

  /* Properties of neutron star mergers model */
  if (fp->with_r_process_enrichment)
    message(
        "Running COLIBRE with r-process enrichment produced by: Neutron star "
        "mergers, Common envelope jets SN and Collapsars");

  fp->NSM_per_Msun =
      parser_get_param_double(params, "COLIBREFeedback:num_of_NSM_per_Msun");
  fp->yield_Eu_from_NSM = parser_get_param_double(
      params, "COLIBREFeedback:yield_Eu_from_NSM_event_Msun");
  fp->CEJSN_per_Msun =
      parser_get_param_double(params, "COLIBREFeedback:num_of_CEJSN_per_Msun");
  fp->yield_Eu_from_CEJSN = parser_get_param_double(
      params, "COLIBREFeedback:yield_Eu_from_CEJSN_event_Msun");
  fp->collapsar_per_Msun = parser_get_param_double(
      params, "COLIBREFeedback:num_of_collapsar_per_Msun");
  fp->yield_Eu_from_collapsar = parser_get_param_double(
      params, "COLIBREFeedback:yield_Eu_from_collapsar_event_Msun");

  /* Stellar mass limits for collapsar events */
  const double collapsar_min_mass_msun = parser_get_param_double(
      params, "COLIBREFeedback:collapsar_min_mass_Msun");
  const double collapsar_max_mass_msun = parser_get_param_double(
      params, "COLIBREFeedback:collapsar_max_mass_Msun");

  /* Check that it makes sense. */
  if (collapsar_max_mass_msun < collapsar_min_mass_msun) {
    error(
        "Can't have the max collapsar mass smaller than the min collapsar "
        "mass!");
  }

  fp->log10_collapsar_min_mass_msun = log10(collapsar_min_mass_msun);
  fp->log10_collapsar_max_mass_msun = log10(collapsar_max_mass_msun);

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Read AGB ejecta velocity */
  const float ejecta_velocity_km_p_s = parser_get_param_float(
      params, "COLIBREFeedback:AGB_ejecta_velocity_km_p_s");

  /* Convert to internal units */
  const float ejecta_velocity_cgs = ejecta_velocity_km_p_s * 1e5;
  const float ejecta_velocity =
      ejecta_velocity_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Convert to specific thermal energy */
  fp->AGB_ejecta_specific_kinetic_energy =
      0.5f * ejecta_velocity * ejecta_velocity;

  /* Properties of the HII regions and stellar winds model ------------------ */

  if (fp->with_HIIRegions)
    message("Running COLIBRE feedback with early feedback: HII regions");

  if (fp->with_StellarWinds)
    message("Running COLIBRE feedback with early feedback: Stellar Winds");

  if (fp->with_HIIRegions || fp->with_StellarWinds) {
    parser_get_param_string(params, "COLIBREFeedback:earlyfb_filename",
                            fp->early_feedback_table_path);

    /* Read the HII tables */
    read_feedback_tables(fp);

    /* get the optional timescales, or set them from the tables by default */

    fp->HIIregion_max_age_Myr = parser_get_opt_param_float(
        params, "COLIBREFeedback:HIIregion_maxage_Myr", default_maxage_Myr_HII);

    if (fp->HIIregion_max_age_Myr == 0)
      error(
          "HIIregion_maxage_Myr can't be 0. Consider turning off HII regions "
          "by setting with_HIIRegions: 0"
          "in the parameter file");

    fp->SW_max_age_Myr = parser_get_opt_param_float(
        params, "COLIBREFeedback:stellarwind_maxage_Myr",
        default_maxage_Myr_SW);

    if (fp->SW_max_age_Myr == 0)
      error(
          "stellarwind_maxage_Myr can't be 0. Consider turning off stellar "
          "winds by setting with_StellarWinds: 0"
          "in the parameter file");

    fp->delta_v = parser_get_param_double(
        params, "COLIBREFeedback:Momentum_desired_delta_v");

    fp->HIIregion_dt_Myr = parser_get_param_float(
        params, "COLIBREFeedback:HIIregion_rebuild_dt_Myr");

    /* set the minimum and maximum metallicities */
    fp->Zmin_early_fb = exp10(fp->HII_log10_Zbins[0]);
    fp->Zmax_early_fb = exp10(fp->HII_log10_Zbins[fp->HII_nr_metbins - 1]);

    /* check that we won't exceed the maximum stellar age in the table */
    if (fp->HIIregion_max_age_Myr > fp->HII_agebins[fp->HII_nr_agebins - 1])
      error(
          "HIIregion_maxage_Myr (%.2f Myr) exceeds maximum age in table (%.2f "
          "Myr)",
          fp->HIIregion_max_age_Myr, fp->HII_agebins[fp->HII_nr_agebins - 1]);

    if (fp->SW_max_age_Myr > fp->HII_agebins[fp->HII_nr_agebins - 1])
      error(
          "stellarwind_maxage_Myr (%.2f Myr) exceeds maximum age in table "
          "(%.2f Myr)",
          fp->SW_max_age_Myr, fp->HII_agebins[fp->HII_nr_agebins - 1]);

  } else {

    /* Initialize to zero if run without early feedback */
    fp->HIIregion_max_age_Myr = 0.;
    fp->HIIregion_dt_Myr = 0.;
    fp->delta_v = 0.;
    fp->SW_max_age_Myr = 0.;
  }

  /* Energy injected by the HII regions */
  const double HII_region_temp =
      parser_get_param_double(params, "COLIBREFeedback:HIIregion_temperature");
  const double HIIregion_fion = parser_get_param_double(
      params, "COLIBREFeedback:HIIregion_ionization_fraction");
  const double XH = 1. - phys_const->const_primordial_He_fraction;
  const double mu_HII = 4.0 / ((1.0 + HIIregion_fion) * (1.0 + (3.0 * XH)));
  fp->HII_u = HII_region_temp * phys_const->const_boltzmann_k /
              (phys_const->const_proton_mass * hydro_gamma_minus_one * mu_HII);

  /* Properties of the enrichment down-sampling ----------------------------- */

  const double stellar_evolution_age_cut_Gyr = parser_get_param_double(
      params, "COLIBREFeedback:stellar_evolution_age_cut_Gyr");

  fp->stellar_evolution_age_cut =
      stellar_evolution_age_cut_Gyr * Gyr_in_cgs /
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  fp->stellar_evolution_sampling_rate = parser_get_param_double(
      params, "COLIBREFeedback:stellar_evolution_sampling_rate");

  if (fp->stellar_evolution_sampling_rate < 1 ||
      fp->stellar_evolution_sampling_rate >= (1 << (8 * sizeof(char) - 1)))
    error("Stellar evolution sampling rate too large. Must be >0 and <%d",
          (1 << (8 * sizeof(char) - 1)));

  /* Make sure the cut does not happen before the end of the early feedback */
  if (stellar_evolution_age_cut_Gyr * 1000. < fp->HIIregion_max_age_Myr)
    error(
        "Downsampling the stellar evolution calculation at ages lower than the "
        "max HII region age is forbidden. Increase the "
        "stellar_evolution_age_cut_Gyr parameter.");

  if (stellar_evolution_age_cut_Gyr * 1000. < fp->SW_max_age_Myr)
    error(
        "Downsampling the stellar evolution calculation at ages lower than the "
        "max stellar wind age is forbidden. Increase the "
        "stellar_evolution_age_cut_Gyr parameter.");

  /* Gather common conversion factors --------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Conversion factor to CGS for momentum */
  const double momentum_factor =
      units_cgs_conversion_factor(us, UNIT_CONV_MOMENTUM);
  fp->Momentum_to_cgs = momentum_factor;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Get recombination coefficient in cgs units */
  const float dimension_alphaB[5] = {0, 3, -1, 0, 0}; /* [cm^3 s^-1] */
  fp->alpha_caseb_recomb =
      phys_const->const_caseb_recomb *
      units_general_cgs_conversion_factor(us, dimension_alphaB);

  fp->Myr_to_sec = 1.e6 * phys_const->const_year *
                   units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  fp->sec_to_Myr = 1. / fp->Myr_to_sec;

  /* Initialise the IMF ------------------------------------------------- */

  init_imf(fp);

  /* Calculate number of type II SN per unit solar mass based on our choice
   * of IMF and integration limits for type II SNe.
   * Note: No weighting by yields here. */
  fp->num_SNII_per_msun =
      integrate_imf(fp->log10_SNII_min_mass_msun, fp->log10_SNII_max_mass_msun,
                    eagle_imf_integration_no_weight,
                    /*(stellar_yields=)*/ NULL, fp);

  /* Initialise the yields ---------------------------------------------- */

  /* Read yield table filepath  */
  parser_get_param_string(params, "COLIBREFeedback:filename",
                          fp->yield_table_path);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the tables  */
  read_yield_tables(fp);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp); // <HERE>

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);

  message("initialized stellar feedback");
}

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param feedback_props the feedback data structure.
 */
void feedback_clean(struct feedback_props* fp) {

  free(fp->HII_log10_Zbins);
  free(fp->HII_agebins);
  free(fp->HII_log10_Qcum);
  free(fp->SW_log10_Pcum);
  swift_free("imf-tables", fp->imf);
  swift_free("imf-tables", fp->imf_mass_bin);
  swift_free("imf-tables", fp->imf_mass_bin_log10);
  swift_free("feedback-tables", fp->yields_SNIa);
  swift_free("feedback-tables", fp->yield_SNIa_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.mass);
  swift_free("feedback-tables", fp->yield_AGB.metallicity);
  swift_free("feedback-tables", fp->yield_AGB.yield);
  swift_free("feedback-tables", fp->yield_AGB.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.ejecta);
  swift_free("feedback-tables", fp->yield_AGB.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.total_metals);
  swift_free("feedback-tables", fp->yield_AGB.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.mass);
  swift_free("feedback-tables", fp->yield_SNII.metallicity);
  swift_free("feedback-tables", fp->yield_SNII.yield);
  swift_free("feedback-tables", fp->yield_SNII.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.ejecta);
  swift_free("feedback-tables", fp->yield_SNII.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.total_metals);
  swift_free("feedback-tables", fp->yield_SNII.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->lifetimes.mass);
  swift_free("feedback-tables", fp->lifetimes.metallicity);
  swift_free("feedback-tables", fp->yield_mass_bins);
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    free(fp->lifetimes.dyingtime[i]);
  }
  free(fp->lifetimes.dyingtime);
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    free(fp->SNIa_element_names[i]);
  }
  free(fp->SNIa_element_names);
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    free(fp->SNII_element_names[i]);
  }
  free(fp->SNII_element_names);
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    free(fp->AGB_element_names[i]);
  }
  free(fp->AGB_element_names);
}

/**
 * @brief Zero pointers in yield_table structs
 *
 * @param table yield_table struct in which pointers to tables
 * set to NULL
 */
void zero_yield_table_pointers(struct yield_table* table) {

  table->mass = NULL;
  table->metallicity = NULL;
  table->yield_IMF_resampled = NULL;
  table->yield = NULL;
  table->ejecta_IMF_resampled = NULL;
  table->ejecta = NULL;
  table->total_metals_IMF_resampled = NULL;
  table->total_metals = NULL;
}

/**
 * @brief Restore feedback tables (if applicable) after
 * restart
 *
 * @param fp the #feedback_props structure
 */
void feedback_restore_tables(struct feedback_props* fp) {

  init_imf(fp);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the yield tables  */
  read_yield_tables(fp);

  /* Read the HII tables */
  if (fp->with_HIIRegions || fp->with_StellarWinds) {
    message("Reading early feedback tables");
    read_feedback_tables(fp);
  }

  if (fp->with_HIIRegions)
    message("Running COLIBRE feedback with early feedback: HII regions");

  if (fp->with_StellarWinds)
    message("Running COLIBRE feedback with early feedback: Stellar Winds");

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* zero AGB and SNII table pointers */
  zero_yield_table_pointers(&feedback_copy.yield_AGB);
  zero_yield_table_pointers(&feedback_copy.yield_SNII);

  /* zero SNIa table pointers */
  feedback_copy.yield_SNIa_IMF_resampled = NULL;
  feedback_copy.yields_SNIa = NULL;
  feedback_copy.yield_SNIa_total_metals_IMF_resampled = 0;

  /* zero element name tables */
  feedback_copy.SNIa_element_names = NULL;
  feedback_copy.SNII_element_names = NULL;
  feedback_copy.AGB_element_names = NULL;

  /* zero mass bins table */
  feedback_copy.yield_mass_bins = NULL;

  /* zero lifetime tracks */
  feedback_copy.lifetimes.mass = NULL;
  feedback_copy.lifetimes.metallicity = NULL;
  feedback_copy.lifetimes.dyingtime = NULL;

  /* zero IMF tables */
  feedback_copy.imf = NULL;
  feedback_copy.imf_mass_bin = NULL;
  feedback_copy.imf_mass_bin_log10 = NULL;

  /* zero the HII tables */
  feedback_copy.HII_log10_Zbins = NULL;
  feedback_copy.HII_log10_Qcum = NULL;
  feedback_copy.SW_log10_Pcum = NULL;
  feedback_copy.HII_agebins = NULL;

  restart_write_blocks((void*)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {
  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  feedback_restore_tables(feedback);
}
