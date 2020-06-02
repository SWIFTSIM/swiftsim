/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Joel Pfeffer (j.l.pfeffer@ljmu.ac.uk)
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
#ifndef SWIFT_MOSAICS_CLEVO_H
#define SWIFT_MOSAICS_CLEVO_H

#include <math.h>
#include <float.h>

/* Local includes */
#include "dsyevj3.h"

/**
 * @brief Return the cluster evaporation timescale, t0
 *
 * @param sp The #spart to act upon.
 * @param props Properties of the stars model.
 * @param tensors_in_Gyr Convert tidal tensor to Gyr^-2
 */
__attribute__((always_inline)) INLINE static double get_local_evaporation_t0(
    struct spart* restrict sp, const struct stars_props* props,
    const double tensors_in_Gyr) {

  /* The lower elements are never referenced by dsyevj3 */
  double tensor[3][3];
  tensor[0][0] = (double)sp->tidal_tensor[2][0] * tensors_in_Gyr;
  tensor[0][1] = (double)sp->tidal_tensor[2][1] * tensors_in_Gyr;
  tensor[0][2] = (double)sp->tidal_tensor[2][2] * tensors_in_Gyr;
  tensor[1][1] = (double)sp->tidal_tensor[2][3] * tensors_in_Gyr;
  tensor[1][2] = (double)sp->tidal_tensor[2][4] * tensors_in_Gyr;
  tensor[2][2] = (double)sp->tidal_tensor[2][5] * tensors_in_Gyr;

  /* get the eigenvalues, then sort */
  double tide_val[3];
  dsyevj3(tensor, tide_val);
  sort3(tide_val);

  //TODO option to exclude Omega ?
  /* Strength of the tidal field, T + Omega^2 */
  double tidal_strength = tide_val[2] +
          fabs(-tide_val[2] - tide_val[1] - tide_val[0]) / 3.f;

  /* Timescale of evaporation in the tidal field */
  double t0evap;

  if (tidal_strength > 0.) {

    /* Dissolution timescales for circular orbit at Solar radius */
    /* King profile W0=5, 21.3e6 yr in Gyr */
    const double t0evapsunW5 = 0.0213;

    /* King profile W0=7, 10.7e6 yr in Gyr */
    const double t0evapsunW7 = 0.0107;

    /* Approximate tidal field strength at solar galactocentric radius *
     * 7.0342e-31 s^-2 */
    const double tidesun = 700.5234;

    /* Interpolate between King profiles */
    double t0evapsun = t0evapsunW5 + 
        0.5 * (props->W0 - 5.) * (t0evapsunW7 - t0evapsunW5);

    t0evap = t0evapsun / sqrt(tidal_strength / tidesun);

  } else {
    /* Compressive field */
    t0evap = FLT_MAX;
  }

  return t0evap;
}

/**
 * @brief Calculate the tidal heating parameters
 *
 * @param sp The particle to act upon.
 * @param age Age of star.
 * @param dt Particle timestep.
 * @param time Current time.
 * @param tensors_in_Gyr Convert tidal tensor to Gyr^-2
 * @param Mlheatsum (return) Tidal heating for shocks at this timestep.
 */
__attribute__((always_inline)) INLINE static void calc_heating(
    struct spart* restrict sp, const double age, const double dt, 
    const double time, const double tensors_in_Gyr, double* Mlheatsum) {

  /* fraction of maximum at which duration is measured */
  const double sigma1_width = 0.88;
  const double sigma2_width = 0.61;

  /* Time at the previous step */
  const float tlast = (float)(time - dt);

  /* loop over tidal tensor components, to determine minima, maxima, tshock 
   * and tidal heating for each component.
   * Note we only store upper six components of tidal tensor.
   */
  for (int i = 0; i < 6; i++) {
    const double tide1 = (double)sp->tidal_tensor[0][i] * tensors_in_Gyr;
    const double tide2 = (double)sp->tidal_tensor[1][i] * tensors_in_Gyr;
    const double tide3 = (double)sp->tidal_tensor[2][i] * tensors_in_Gyr;

    int apply_shock = 0;

    /* Is the first point for this particle is a maximum? */
    if (age == dt && fabs(tide3) < fabs(tide2)) {

      /* Time of last maximum is previous time step */
      sp->tmaxsh[i] = tlast;
      sp->tidmax[i] = tide2;
      sp->extreme_max[i] = 1;
      sp->shock_duration[i] = 0.f;
      sp->shock_indicator[i] = 0;
    }

    /* If previous time step marks a valid minimum */
    if (fabs(tide2) < fabs(tide1) && fabs(tide2) < fabs(tide3) &&
        fabs(tide2) < 0.88 * fabs(sp->tidmax[i])) {

      /* if previous extreme was a maximum, or this minimum is lower than
       * previous */
      if (sp->extreme_max[i] || (fabs(tide2) < fabs(sp->tidmin[i]))) {

        /* no shock duration is present (shallow max) */
        if (sp->shock_indicator[i] == 0)
          sp->shock_duration[i] = 0.5 * (tlast - sp->tminsh[i]);

        /* Time of last minimum was previous time step */
        apply_shock = 1;
        sp->tminsh[i] = tlast;
        sp->tidmin[i] = tide2;
        sp->extreme_max[i] = 0;
      }
    }

    /* If previous time step marks a valid maximum */
    if (fabs(tide2) > fabs(tide1) && fabs(tide2) > fabs(tide3) &&
        0.88*fabs(tide2) > fabs(sp->tidmin[i])) {

      //if previous extreme was a minimum, or this maximum is higher than previous
      if (!sp->extreme_max[i] ||
          (fabs(tide2) > fabs(sp->tidmax[i])) ) {

        /* Time of last maximum was previous time step */
        sp->tmaxsh[i] = tlast;
        sp->tidmax[i] = tide2;
        sp->extreme_max[i] = 1;
        sp->shock_duration[i] = 0.f;
        sp->shock_indicator[i] = 0;
      }
    }

    /* First order shock duration */
    if (sp->extreme_max[i] && sp->shock_indicator[i] == 0 &&
          fabs(tide3) < sigma1_width * fabs(sp->tidmax[i])) {

      sp->shock_indicator[i] = 1;

      /* time interpolation fraction */
      const double dt_fraction =
          (sigma1_width * fabs(sp->tidmax[i]) - fabs(tide2)) /
          (fabs(tide3) - fabs(tide2));

      /* interpolate to time of shock duration measurement */
      const double shock_time = time + (dt_fraction - 1.f) * dt;

      /* set shock duration indicator to total width of 1 sigma */
      /* 2 * dt / sigma(=1) */
      sp->shock_duration[i] = 2.f * (shock_time - sp->tmaxsh[i]);
    }

    /* Second order shock duration */
    if (sp->extreme_max[i] && sp->shock_indicator[i] < 2 &&
          fabs(tide3) < sigma2_width * fabs(sp->tidmax[i])) {

      sp->shock_indicator[i] = 2;

      /* time interpolation fraction */
      const double dt_fraction =
          (sigma2_width * fabs(sp->tidmax[i]) - fabs(tide2)) /
          (fabs(tide3) - fabs(tide2));

      /* interpolate to time of shock duration measurement */
      const double shock_time = time + (dt_fraction - 1.f) * dt;

      /* set shock duration indicator to total width of 2 sigma */
      /* 2 * dt / sigma(=2) */
      sp->shock_duration[i] = shock_time - sp->tmaxsh[i];
    }

    /* If a shock was completed during the previous time step and mass loss
     * should be applied */
    if (apply_shock) {
      /* Tidal heating for the soon-to-occur application of mass loss
       * Note this could be negative, but is later applied as heatsum^2 */
      Mlheatsum[i] = sp->heatsum[i];

      /* reset and start integrating heating for next shock */
      sp->shock_indicator[i] = 0;
      sp->heatsum[i] = 0.5 * (tide2 + tide3) * dt;
    } else {
      /* Keep summing tidal heating */
      sp->heatsum[i] += 0.5 * (tide2 + tide3) * dt;
    }
  } /* For each tensor component */
}

/**
 * @brief Compute the Weinberg adiabatic correction
 *
 * @param mass Cluster mass in Solar masses.
 * @param rh Cluster half-mass radius in pc.
 * @param tshock Tidal shock duration.
 * @param W0 King profile parameter.
 * @param etapot Potential energy constant.
 */
__attribute__((always_inline)) INLINE static double adiabatic_correction(
    const double mass, const double rh, const double tshock, const double W0,
    const double etapot) {

  /* Gravity pc^3 Msun^-1 Gyr^-2 */
  const double G = 4498.568;

  /* constant in omega for King profile W0=[5,7] */
  const double cst5 = 0.68;
  const double cst7 = 0.82;

  /* Interpolate */
  const double cst = cst5 + 0.5 * (W0 - 5.f) * (cst7 - cst5);

  /* angular frequency of stars */
  const double omega = cst * 
      sqrt(8.f * etapot * etapot * etapot * G * mass / (rh * rh * rh));

  /* number of stellar revolutions per passage */
  const double x = omega * tshock;

  double Aw = 1.f / ((1. + x * x) * sqrt(1. + x * x));

  if (!isfinite(Aw)) {
    error("Adiabatic correction is not finite!");
  }

  return Aw;
}

/**
 * @brief Compute star cluster mass loss by tidal shocks.
 *
 * @param sp The #spart.
 * @param props Properties of the stars model.
 * @param Mlheatsum Tidal heating to apply at this timestep.
 * @param mass Star cluster mass in Solar masses.
 */
__attribute__((always_inline)) INLINE static double mass_loss_shocks(
    struct spart* restrict sp, const struct stars_props* props,
    const double* Mlheatsum, const double mass) {

  /* Anything to do? */

  double max_heatsum = 0.f;
  for (int i = 0; i < 6; i++) {
    if (fabs(Mlheatsum[i]) > max_heatsum)
      max_heatsum = fabs(Mlheatsum[i]);
  }

  if (max_heatsum == 0.)
    return 0.f;

  double max_shock_duration = 0.f;
  for (int i = 0; i < 6; i++) {
    if (fabs(sp->shock_duration[i]) > max_shock_duration)
      max_shock_duration = fabs(sp->shock_duration[i]);
  }

  if (max_shock_duration == 0.)
    return 0.f;

  /* Ok, now complete the shock */

  /* Gravity pc^3 Msun^-1 Gyr^-2 */
  const double G = 4498.568;

  /* potential energy constant */
  const double etapot = 0.4;

  /* fraction of shock energy input used in mass loss */
  const double efrac = 0.25;

  /* mean(M(r)*r)/(Mtot*rh) for King profile W0=[5,7] */
  const double zeta5 = 0.81034;
  const double zeta7 = 1.0273;

  /* ratio of square radii for King profile W0=[5,7] */
  const double r2avrh25 = 2.;
  const double r2avrh27 = 3.5;

  /* Now interpolate */
  const double zeta = zeta5 + 0.5 * (props->W0 - 5.f) * (zeta7 - zeta5);
  const double r2avrh2 =
      r2avrh25 + 0.5 * (props->W0 - 5.f) * (r2avrh27 - r2avrh25);

  /* Half-mass radius in pc */
  const double rh = props->rh;

  /* Sum up the tidal heating */
  double heatparam = 0.f;
  for (int i = 0; i < 6; i++) {
    /* Weinberg adiabatic correction */
    double Aw = adiabatic_correction(mass, rh, (double)sp->shock_duration[i], props->W0, 
        etapot);

    /* Add the mirror for off-diagonal tensor terms */
    double factor = 1.f;
    if (i==1 || i==2 || i==5)
      factor = 2.f;

    heatparam += factor * Mlheatsum[i] * Mlheatsum[i] * Aw;
  }

  /* first order tsh/P */
  const double t1 = 3. * etapot * G * mass / 
      (efrac * rh * rh * rh * r2avrh2 * heatparam);

  /* factor to include second order energy gain */
  const double InclE2 = 1.f + (12.f * zeta) / (5.f * etapot * r2avrh2);

  /* Disruption time / shock period: tsh_int = t1 / InclE2 */
  return -mass * InclE2 / t1;
}

/**
 * @brief Compute star cluster mass loss by evaporation.
 *
 * @param props Properties of the stars model.
 * @param mass Star cluster mass in Solar masses.
 * @param t0 Evaporation timescale in Gyr.
 * @param dt Particle timestep in Gyr.
 */
__attribute__((always_inline)) INLINE static double mass_loss_evaporation(
    const struct stars_props* props, const double mass, const double t0, 
    const double dt) {

  /* Gravity in pc^3 Msun^-1 Gyr^-2 */
  const double G = 4498.568;

  /* Dissolution law exponent for W0 = [5,7] King profile */
  const double gamma5 = 0.62;
  const double gamma7 = 0.70;

  /* adapt dissolution law exponent to King profile parameter */
  const double gamma = gamma5 + 0.5 * (props->W0 - 5.) * (gamma7 - gamma5);

  double dm = 0.f;

  /* Add the 'isolated' evaporation term, following Gieles & Baumgardt 08 */
  if (props->spitzer_evap_term) {

    /* Mass loss fraction per relaxation time in 'isolated' regime
     * (Spitzer 87, Gieles & Baumgardt 08) */
    const double xi0 = 0.0074;

    /* Mean star mass for Kroupa IMF */
    const double N = mass / 0.64;

    /* Cluster half-mass radius in pc */
    const double rh = props->rh;

    /* Relaxation time (Spitzer & Hart 71; Gierz & Heggie 94) */
    const double trh = 0.138 * sqrt(N) * rh * sqrt(rh) / 
        ( sqrt(0.64 * G) * log(0.11 * N) );

    dm += -xi0 * mass / trh * dt;
  }

  const double dt_t0 = dt / t0;

  /* Was it mass loss? Set to first order if mass gain.
   * dm>0 only occurs just before disruption because we don't consider
   * 3rd order effects
   */
  double dm_order1 = -pow(mass, 1. - gamma) * dt_t0;
  double dm_order2 = (1. - gamma) * pow(mass, 1. - 2. * gamma) * dt_t0 * dt_t0;

  if (dm_order1 + dm_order2 > 0.)
    dm += dm_order1;
  else
    dm += dm_order1 + dm_order2;

  return dm;
}

/**
 * @brief Do star cluster evolution for the mosaics subgrid model
 *
 * @param sp The #spart to act upon.
 * @param e The #engine.
 * @param with_cosmology Are we running with cosmological time integration.
 */
__attribute__((always_inline)) INLINE static void mosaics_clevo(
    struct spart* restrict sp, const struct engine* e, 
    const int with_cosmology) {

  const struct stars_props* props = e->stars_properties;
  const struct cosmology *cosmo = e->cosmology;

  /* Initial number of clusters, up to max array length */
  const int ninit = min(sp->initial_num_clusters_evo, MOSAICS_MAX_CLUSTERS);

  /* Proportional change in mass from stellar evolution */
  const float stellar_evo_factor = sp->mass / sp->mass_prev_timestep;

  if (ninit == 0) {
    /* No clusters formed in this particle, so just handle the field */
    sp->field_mass *= stellar_evo_factor;

    /* For stellar evo. mass loss at the next timestep */
    sp->mass_prev_timestep = sp->mass;
    return;
  }

  /* Any surviving clusters? */
  if (sp->num_clusters > 0) {

    const integertime_t ti_step = get_integer_timestep(sp->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current - 1, sp->time_bin);

    /* Calculate age of the star at current time */
    double star_age;
    if (with_cosmology) {
      star_age = cosmology_get_delta_time_from_scale_factors(
              cosmo, (double)sp->birth_scale_factor, cosmo->a);
    } else {
      star_age = e->time - (double)sp->birth_time;
    }

    /* Get particle time-step */
    double dt_star;
    if (with_cosmology) {
      dt_star = cosmology_get_delta_time(e->cosmology, ti_begin,
                                         ti_begin + ti_step);
    } else {
      dt_star = get_timestep(sp->time_bin, e->time_base);
    }

    /* Convert dt and stellar age from internal units to Gyr. */
    const double Gyr_in_cgs = 1e9 * 365.25 * 24. * 3600.;
    const double conversion_factor = props->time_to_cgs / Gyr_in_cgs;
    const double dt_star_Gyr = dt_star * conversion_factor;
    const double star_age_Gyr = star_age * conversion_factor;
    const double time = e->time * conversion_factor;

    /* Unit conversion for tensors and factor out cosmology, in Gyr */
    const double tensors_in_Gyr = cosmo->a3_inv *
        props->tidal_tensor_to_cgs * Gyr_in_cgs * Gyr_in_cgs;

    /* Evaporation timescale constant */
    double t0_evap;
    if (props->fixed_t0evap > 0) {
      /* Use a fixed value? */
      t0_evap = props->fixed_t0evap;
    } else {
      /* Get t0 from local tidal field */
      t0_evap = get_local_evaporation_t0(sp, props, tensors_in_Gyr);
    }

    /* Calculate the tidal heating parameters for this timestep */
    double Mlheatsum[6] = {0.f};
    calc_heating(sp, star_age_Gyr, dt_star_Gyr, time, tensors_in_Gyr, 
        Mlheatsum);

    /* Evolve the clusters */
    for (int i = 0; i < ninit; i++) {

      if (sp->clusters.mass[i] == 0.f)
        continue;

      const double mass = 
          (double)sp->clusters.mass[i] * props->mass_to_solar_mass;

      double dmevap = 0.f, dmsh = 0.f;

      /* Mass loss by tidal shocks */
      if (props->cluster_shocks)
        dmsh = mass_loss_shocks(sp, props, Mlheatsum, mass);

      /* Mass loss by evaporation */
      if (props->cluster_evap)
        dmevap = mass_loss_evaporation(props, mass, t0_evap, dt_star_Gyr);

      /* Total mass loss for this cluster */
      double dmdis = dmevap + dmsh;

      /* Now back to internal units */
      dmdis *= props->solar_mass_to_mass;
      dmevap *= props->solar_mass_to_mass;
      dmsh *= props->solar_mass_to_mass;

      /* Update cluster mass */
      if (sp->clusters.mass[i] + (float)dmdis < props->clMF_min) {
        /* Cluster went below minimum mass limit, consider fully disrupted */

        /* fraction of mass loss */
        dmevap *= (double)(-sp->clusters.mass[i]) / dmdis;
        dmsh *= (double)(-sp->clusters.mass[i]) / dmdis;
  
        /* new mass loss (can't lose more than it has) */
        dmdis = (double)(-sp->clusters.mass[i]);

        /* We're done with this one */
        sp->clusters.mass[i] = 0.f;
        sp->num_clusters--;

        /* time of cluster disruption */
        if (with_cosmology) {
          sp->clusters.disruption_time[i] = cosmo->a;
        } else {
          sp->clusters.disruption_time[i] = e->time;
        }

      } else {
  
        /* Cluster still survives, for now */
        sp->clusters.mass[i] += (float)dmdis;
      }

      //TODO
      if (dmevap > 0 || dmsh > 0) {
        printf("Positive mass loss?? dmevp = %g; dmsh = %g\n", dmevap, dmsh);
      }

      /* Cluster mass lost to field */
      sp->clusters.dmevap[i] += (float)dmevap;
      sp->clusters.dmshock[i] += (float)dmsh;
      sp->field_mass += (float)fabs(dmdis);

      /* Stellar evolutionary mass loss for clusters */
      sp->clusters.mass[i] *= stellar_evo_factor;

    } /* Cluster loop */

    /* Re-sum the cluster masses */
    sp->cluster_mass_total = 0.f;
    for (int i = 0; i < ninit; i++) {
      sp->cluster_mass_total += sp->clusters.mass[i];
    }

  } /* surviving clusters */

  /* Apply stellar evolution to the 'field' populations */
  sp->field_mass *= stellar_evo_factor;
  for (int i = 0; i < ninit; i++) {
    sp->clusters.dmevap[i] *= stellar_evo_factor;
    sp->clusters.dmshock[i] *= stellar_evo_factor;
  }

  /* For stellar evo. mass loss at the next timestep */
  sp->mass_prev_timestep = sp->mass;
}

#endif /* SWIFT_MOSAICS_CLEVO_H */
