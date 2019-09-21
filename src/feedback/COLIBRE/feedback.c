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
#include "feedback.h"

/* Local includes. */
#include "hydro_properties.h"
#include "imf.h"
#include "inline.h"
#include "interpolate.h"
#include "physical_constants.h"
#include "timers.h"
#include "yield_tables.h"

/* Minimal/maximal value of the metallicity (metal mass fraction)
 * available in the Starburst99 model */
static const double colibre_feedback_momentum_SB99_Z_min = 0.001;
static const double colibre_feedback_momentum_SB99_Z_max = 0.04;

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
  const double rho_birth = sp->birth_density;
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
 * @param feedback_props The properties of the feedback model.
 * @param age of star particle at the beginning of the timestep
 * @param timestep in Gyr
 * @param min_dying_mass_Msun stellar mass that dies at the end of this timestep
 * @param max_dying_mass_Msun stellar mass that dies at the beginning of this
 * timestep
 */
INLINE static void compute_SNII_feedback(
    struct spart* sp, const double star_age, const double dt,
    const float ngb_gas_mass, const struct feedback_props* feedback_props,
    const float min_dying_mass_Msun, const float max_dying_mass_Msun) {

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
    double N_SNe;

    if (SNII_wind_delay > 0.) {
      N_SNe = N_SNe_eagle;
    } else {
      N_SNe = N_SNe_colibre;
    }

    const double E_SNe = feedback_props->E_SNII;
    const double f_E = eagle_SNII_feedback_energy_fraction(sp, feedback_props);

    /* Conversion factor from T to internal energy */
    const double conv_factor = feedback_props->temp_to_u_factor;

    /* Calculate the default heating probability */
    double prob = f_E * E_SNe * N_SNe / (conv_factor * delta_T * ngb_gas_mass);

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
      delta_u = f_E * E_SNe * N_SNe / ngb_gas_mass;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (f_E < feedback_props->f_E_min || f_E > feedback_props->f_E_max)
      error("f_E is not in the valid range! f_E=%f sp->id=%lld", f_E, sp->id);
#endif

    /* Store all of this in the star for delivery onto the gas */
    sp->SNII_f_E = f_E;
    sp->feedback_data.to_distribute.SNII_heating_probability = prob;
    sp->feedback_data.to_distribute.SNII_delta_u = delta_u;
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

  /* Conversion factor from T to internal energy */
  const double conv_factor = feedback_props->temp_to_u_factor;

  /* Calculate the default heating probability */
  double prob = f_E * E_SNe * N_SNe / (conv_factor * delta_T * ngb_gas_mass);

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
    delta_u = f_E * E_SNe * N_SNe / ngb_gas_mass;
    if (ngb_gas_mass == 0.) {
      prob = 0.;
      delta_u = 0.;
    } else {
      message("WOW the probability is so high! %e %e %e %e %e %e %llu", f_E,
              E_SNe, N_SNe, conv_factor, delta_T, ngb_gas_mass, sp->id);
    }
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

  /* Check if we're outside the mass range for SNIa */
  if (log10_min_mass >= props->log10_SNIa_max_mass_msun) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (dt_Gyr < 0.) error("Negative time-step length!");
  if (star_age_Gyr < 0.) error("Negative age!");
#endif

  /* If the max mass is outside the mass range update it to be the maximum
   * and use updated values for the star's age and timestep in this function */
  if (log10_max_mass > props->log10_SNIa_max_mass_msun) {

    const float Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);
    const float max_mass = exp10f(props->log10_SNIa_max_mass_msun);
    const float lifetime_Gyr = lifetime_in_Gyr(max_mass, Z, props);

    dt_Gyr = max(star_age_Gyr + dt_Gyr - lifetime_Gyr, 0.);
    star_age_Gyr = lifetime_Gyr;
  }

  /* Compute the number of SNIa */
  const float num_SNIa =
      dtd_number_of_SNIa(sp, star_age_Gyr, star_age_Gyr + dt_Gyr, props);

  /* compute mass of each metal */
  for (int i = 0; i < chemistry_element_count; i++) {
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
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
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

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total metal mass ejected from the star*/
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

    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.metal_mass[i] +=
          metal_mass_released[i] * norm_factor;
    }
    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.mass_from_SNII +=
          sp->feedback_data.to_distribute.metal_mass[i];
    }
    sp->feedback_data.to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->feedback_data.to_distribute.metal_mass_from_SNII +=
        metal_mass_released_total * norm_factor;
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
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
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

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total metal mass ejected from the star */
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

    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.metal_mass[i] +=
          metal_mass_released[i] * norm_factor;
      sp->feedback_data.to_distribute.mass_from_AGB +=
          metal_mass_released[i] * norm_factor;
    }
    sp->feedback_data.to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->feedback_data.to_distribute.metal_mass_from_AGB +=
        metal_mass_released_total * norm_factor;
  } else {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
}

/**
 * @brief Gets interpolated cumulative ionizing photons from table
 * @param props feedback_props structure for getting model parameters for
 * coefficients
 * @param t_Myr stellar age in Myr
 * @param logZ log10 of (stellar) metal mass fraction Z
 */
double get_cumulative_ionizing_photons(const struct feedback_props* fp,
                                       float t_Myr, float logZ) {
  float logQcum_loc;
  float d_age, d_met;
  int met_index, age_index;

  if (t_Myr < fp->HII_agebins[0]) return 0.;

  get_index_1d(fp->HII_logZbins, fp->HII_nr_metbins, logZ, &met_index, &d_met);

  get_index_1d(fp->HII_agebins, fp->HII_nr_agebins, t_Myr, &age_index, &d_age);

  logQcum_loc =
      interpolation_2d_flat(fp->HII_logQcum, met_index, age_index, d_met, d_age,
                            fp->HII_nr_metbins, fp->HII_nr_agebins);

  return exp10(logQcum_loc);
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

  double Q_t1, Q_t2, Qbar;
  /* find a way to get that from constants, if possible without passing the
     whole structure through everything...*/
  const double Myr_inv = 3.1688e-14;

  Q_t1 = get_cumulative_ionizing_photons(fp, t1, log10(Z));
  Q_t2 = get_cumulative_ionizing_photons(fp, t2, log10(Z));

  Qbar = (Q_t2 - Q_t1) / (t2 - t1) * Myr_inv;
  return Qbar;
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

  /* From starburst 99 */
  const double p1 = props->p1;            /* g cm s^-1 Mo^-1 */
  const double p2 = props->p2;            /* Metallicity mass faction*/
  const double p3 = props->p3;            /* exponent of metallicity*/
  const double tw = props->tw;            /* Myr */
  double delta_v_km_p_s = props->delta_v; /* km s^-1 */

  /* delta_v in code units */
  double delta_v =
      delta_v_km_p_s /
      (units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) * 1.0e-5);

  /* Unit conversion constant */
  const double Myr_in_s = 1.0e6 * 365 * 24 * 60 * 60.;

  /* Convert the times to the units used my the model */
  const double star_age_Myr = star_age_Gyr * 1e3;
  double dt_cgs = dt * us->UnitTime_in_cgs;
  double dt_Myr = dt_cgs / Myr_in_s;

  /* Prevent star particle from injecting momentum for longer than tw */
  float dt_new = dt;
  if (star_age_Myr + dt_Myr > tw) {
    dt_new = (tw - star_age_Myr) * Myr_in_s / us->UnitTime_in_cgs;
    dt_cgs = dt_new * us->UnitTime_in_cgs;
  }

  /* Star too old to do any momentum injection? */
  if (star_age_Myr > tw) {
    sp->feedback_data.to_distribute.momentum = 0.0;
    sp->feedback_data.to_distribute.momentum_probability = -1;
    sp->feedback_data.to_distribute.momentum_weight = 0.0;
    sp->feedback_data.to_distribute.momentum_delta_v = 0.0;
    return;
  }

  /* Initial stellar mass in Msun */
  const double mstr_in_Msun = sp->mass_init * props->mass_to_solar_mass;

  /* Mass within the SPH kernel in grams */
  const double ngb_gas_mass_in_g = ngb_gas_mass * us->UnitMass_in_cgs;

  /* Star metallicity (metal mass fraction) at birth */
  double Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* Bring the metallicity in the range covered by the model */
  Z = min(Z, colibre_feedback_momentum_SB99_Z_max);
  Z = max(Z, colibre_feedback_momentum_SB99_Z_min);

  /* Apply the Agert+13 model (all in CGS) */

  /* Maximum time over which star can inject momentum */
  const double tw_cgs = tw * Myr_in_s;
  const double p1_cgs = p1 * mstr_in_Msun * pow(Z / p2, p3);

  /* Momentum rate */
  const double dp_dt_cgs = p1_cgs / tw_cgs;

  /* Velocity kick */
  const double delta_v_cgs = delta_v_km_p_s * 1e5;

  /* Get the momentum rate in code units and store it */
  sp->feedback_data.to_distribute.momentum =
      dp_dt_cgs * dt_cgs / props->Momentum_to_cgs;
  sp->feedback_data.to_distribute.momentum_weight = ngb_gas_mass;

  /* Now compute the robability of kicking particle with given delta_v
   * in the current timestep.
   * Note that this could be prop > 1 if there are no enough particles in the
   * kernel to distribute the amount of momentum available in
   * the timestep, but this is OK. */

  double prob = 0.;

  /* We want the code to decide the velocity kick for us */
  if (delta_v < 0.) {

    /* We kick all particles in the kernel */
    prob = 1.;

    /* Kick velocity (in code units) needed so that we can inject
     * all the momentum inside the kernel */
    delta_v = sp->feedback_data.to_distribute.momentum /
              sp->feedback_data.to_distribute.momentum_weight;

    /* User decided velocity kick */
  } else {

    /* The probability of kicking a particle is given
     * by the kicking velocity chosen in the parameter file
     * and is normalised by the mass available in the kernel */
    prob = (dp_dt_cgs * dt_cgs / delta_v_cgs) / ngb_gas_mass_in_g;

    /* Mass inside the kernel too small makes prob > 1 */
    if (prob > 1.) {

      message("Not enough parts in the kernel to distribute momentum...");

      /* Correct the kick (in code units) to be consistent with the mass within
       * the kernel and amount of momentum available */
      delta_v = (sp->feedback_data.to_distribute.momentum /
                 sp->feedback_data.to_distribute.momentum_weight);

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
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const double time_beg_of_step) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");
#endif

  /* Allocate temporary array for calculating imf weights */
  float stellar_yields[eagle_feedback_N_imf_bins];

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_in_cgs = 1e9 * 365. * 24. * 3600.;
  const double Myr_in_cgs = 1e6 * 365. * 24. * 3600.;
  const double time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const double conversion_factor = time_to_cgs / Gyr_in_cgs;
  const double dt_Gyr = dt * conversion_factor;
  const double dt_Myr = dt * conversion_factor * 1e3;
  const double star_age_Gyr = age * conversion_factor;
  const double star_age_Myr = age * conversion_factor * 1e3;

  /* Get the metallicity */
  const float Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

  /* Compute amount of momentum available for this stars, given its mass and age
   */
  compute_stellar_momentum(sp, us, feedback_props, star_age_Gyr, dt,
                           ngb_gas_mass);

  sp->star_timestep = dt;

  /* Compute ionizing photons for HII regions */
  if (feedback_props->with_HIIregions &&
      star_age_Myr <= feedback_props->HIIregion_maxageMyr) {
    /* only rebuild every HIIregion_dtMyr and at the first timestep the star was
     * formed*/
    if ((((sp->HIIregion_last_rebuild + feedback_props->HIIregion_dtMyr) >=
          star_age_Myr) &&
         ((sp->HIIregion_last_rebuild + feedback_props->HIIregion_dtMyr) <
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

      const double rho_birth = (double)sp->birth_density;
      const double n_birth = rho_birth * feedback_props->rho_to_n_cgs;
      const double alpha_B = (double)feedback_props->alpha_caseb_recomb;
      const double t_half = age * time_to_cgs + 0.5 * dt_Myr * Myr_in_cgs;

      double Qbar;
      float t1_Myr = (float)old_star_age_Myr;
      float t2_Myr = (float)star_age_Myr;

      Qbar = compute_average_photoionizing_luminosity(feedback_props, t1_Myr,
                                                      t2_Myr, Z);

      /* masses in system units */
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
          feedback_props->HIIregion_dtMyr * Myr_in_cgs / time_to_cgs;
      sp->feedback_data.to_distribute.HIIregion_endtime =
          time_beg_of_step + HIIregion_dt;
      sp->feedback_data.to_distribute.HIIregion_starid = sp->id;

    } else {
      sp->feedback_data.to_distribute.HIIregion_probability = -1.;
      sp->feedback_data.to_distribute.HIIregion_endtime = -1.;
      sp->feedback_data.to_distribute.HIIregion_starid = -1;
    }

  } else if (feedback_props->with_HIIregions) {
    sp->HIIregion_last_rebuild = -1.;
    sp->feedback_data.to_distribute.HIIregion_probability = -1.;
    sp->feedback_data.to_distribute.HIIregion_endtime = -1.;
    sp->feedback_data.to_distribute.HIIregion_starid = -1;
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
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Compute properties of the stochastic SNII feedback model. */
  if (feedback_props->with_SNII_feedback) {
    compute_SNII_feedback(sp, age, dt, ngb_gas_mass, feedback_props,
                          min_dying_mass_Msun, max_dying_mass_Msun);
  }
  if (feedback_props->with_SNIa_feedback) {
    compute_SNIa_feedback(sp, age, dt, ngb_gas_mass, feedback_props, dt_Gyr,
                          star_age_Gyr);
  }

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
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
                feedback_props, sp, star_age_Gyr, dt_Gyr);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
                stellar_yields, feedback_props, sp);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
               stellar_yields, feedback_props, sp);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He  metals) */
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

  /* Star age in Myr to store in case an SNII event occurs */
  sp->feedback_data.to_distribute.SNII_star_age_Myr = (float)star_age_Myr;

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
                         const struct cosmology* cosmo) {

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

  fp->with_HIIregions =
      parser_get_param_int(params, "COLIBREFeedback:use_HIIregions");

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
  const double Gyr_in_cgs = 1.0e9 * 365. * 24. * 3600.;
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
  fp->p1 = parser_get_param_double(params,
                                   "COLIBREFeedback:Momentum_per_StellarMass");
  fp->p2 = parser_get_param_double(params,
                                   "COLIBREFeedback:Momentum_Metallicity_norm");
  fp->p3 = parser_get_param_double(
      params, "COLIBREFeedback:Momentum_Metallicity_exponent");
  fp->tw =
      parser_get_param_double(params, "COLIBREFeedback:Momentum_time_scale");
  fp->delta_v = parser_get_param_double(
      params, "COLIBREFeedback:Momentum_desired_delta_v");

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
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "COLIBREFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  const double SNIa_max_mass_msun =
      parser_get_param_double(params, "COLIBREFeedback:SNIa_max_mass_Msun");
  fp->log10_SNIa_max_mass_msun = log10(SNIa_max_mass_msun);

  /* Load the SNIa model */
  dtd_init(fp, phys_const, us, params);

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "COLIBREFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

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

  /* Properties of the HII region model ------------------------------------- */
  fp->HIIregion_maxageMyr =
      parser_get_param_float(params, "COLIBREFeedback:HIIregion_maxage_Myr");

  fp->HIIregion_dtMyr = parser_get_param_float(
      params, "COLIBREFeedback:HIIregion_rebuild_dt_Myr");

  /* Read the HII table */
  if (fp->with_HIIregions) {

    /* Read yield table filepath  */
    parser_get_param_string(params, "COLIBREFeedback:earlyfb_filename",
                            fp->early_feedback_table_path);
#ifdef HAVE_HDF5
    hid_t dataset;
    herr_t status;

    hid_t tempfile_id =
        H5Fopen(fp->early_feedback_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (tempfile_id < 0)
      error("unable to open file %s\n", fp->early_feedback_table_path);

    /* read sizes of array dimensions */
    dataset = H5Dopen(tempfile_id, "Header/NMETALLICITIES", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &fp->HII_nr_metbins);
    if (status < 0) error("error reading number of metallicities \n");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset: number of metallicities \n");

    dataset = H5Dopen(tempfile_id, "Header/NAGES", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &fp->HII_nr_agebins);
    if (status < 0) error("error reading number of ages \n");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset: number of ages\n");

    /* allocate arrays */
    if (posix_memalign((void**)&fp->HII_logZbins, SWIFT_STRUCT_ALIGNMENT,
                       fp->HII_nr_metbins * sizeof(float)) != 0)
      error("Failed to allocate metallicity bins\n");
    if (posix_memalign((void**)&fp->HII_agebins, SWIFT_STRUCT_ALIGNMENT,
                       fp->HII_nr_agebins * sizeof(float)) != 0)
      error("Failed to allocate age bins\n");
    if (posix_memalign(
            (void**)&fp->HII_logQcum, SWIFT_STRUCT_ALIGNMENT,
            fp->HII_nr_metbins * fp->HII_nr_agebins * sizeof(float)) != 0)
      error("Failed to allocate Q array\n");

    /* read in the metallicity bins */
    dataset = H5Dopen(tempfile_id, "MetallicityBins", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     fp->HII_logZbins);
    if (status < 0) error("error reading metallicity bins\n");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset: metallicity bins");

    /* read in the stellar ages bins */
    dataset = H5Dopen(tempfile_id, "AgeBins", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     fp->HII_agebins);
    if (status < 0) error("error reading age bins\n");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset: age bins");

    /* read in cumulative ionizing photons */
    dataset = H5Dopen(tempfile_id, "logQcumulative", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     fp->HII_logQcum);
    if (status < 0) error("error reading Q\n");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset: logQcumulative");

    if (fp->HIIregion_maxageMyr > fp->HII_agebins[fp->HII_nr_agebins - 1])
      error("Stopping for now (%.4f is larger than %.4f)",
            fp->HIIregion_maxageMyr, fp->HII_agebins[fp->HII_nr_agebins - 1]);
#else
    error("Need HDF5 to read early feedback tables");
#endif
  }

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
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);

  message("initialized stellar feedback");
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
