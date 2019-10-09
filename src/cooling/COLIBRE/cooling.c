/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 * @file src/cooling/COLIBRE/cooling.c
 * @brief COLIBRE cooling functions
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling.h"
#include "cooling_rates.h"
#include "cooling_struct.h"
#include "cooling_tables.h"
#include "entropy_floor.h"
#include "error.h"
#include "exp10.h"
#include "hydro.h"
#include "interpolate.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/* Maximum number of iterations for
 * bisection integration schemes */
static const int bisection_max_iterations = 150;

/* Tolerances for termination criteria. */
static const float explicit_tolerance = 0.05;
static const float bisection_tolerance = 1.0e-6;
static const double bracket_factor = 1.5; /* sqrt(1.1) */

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * Also calls the additional H reionisation energy injection if need be.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The space data, including a pointer to array of particles
 */
void cooling_update(const struct cosmology *cosmo,
                    struct cooling_function_data *cooling, struct space *s) {

  /* Extra energy for reionization? */
  if (!cooling->H_reion_done) {

    /* Does this timestep straddle Hydrogen reionization? If so, we need to
     * input extra heat */
    if (cosmo->z <= cooling->H_reion_z && cosmo->z_old > cooling->H_reion_z) {

      if (s == NULL) error("Trying to do H reionization on an empty space!");

      /* Inject energy to all particles */
      cooling_Hydrogen_reionization(cooling, cosmo, s);

      /* Flag that reionization happened */
      cooling->H_reion_done = 1;
    }
  }
}

/**
 * @brief Bisection integration scheme
 *
 * @param u_ini_cgs Internal energy at beginning of hydro step in CGS.
 * @param n_H_cgs Hydrogen number density in CGS.
 * @param redshift Current redshift.
 * @param n_H_index Particle hydrogen number density index.
 * @param d_n_H Particle hydrogen number density offset.
 * @param met_index Particle metallicity index.
 * @param d_met Particle metallicity offset.
 * @param red_index Redshift index.
 * @param d_red Redshift offset.
 * @param Lambda_He_reion_cgs Cooling rate coming from He reionization.
 * @param ratefact_cgs Multiplication factor to get a cooling rate.
 * @param cooling #cooling_function_data structure.
 * @param abundance_ratio Array of ratios of metal abundance to solar.
 * @param dt_cgs timestep in CGS.
 * @param ID ID of the particle (for debugging).
 */
static INLINE double bisection_iter(
    const double u_ini_cgs, const double n_H_cgs, const double redshift,
    int n_H_index, float d_n_H, int met_index, float d_met, int red_index,
    float d_red, double Lambda_He_reion_cgs, double ratefact_cgs,
    const struct cooling_function_data *restrict cooling,
    const float abundance_ratio[colibre_cooling_N_elementtypes], double dt_cgs,
    long long ID) {

  /* Bracketing */
  double u_lower_cgs = u_ini_cgs;
  double u_upper_cgs = u_ini_cgs;

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      colibre_cooling_rate(log10(u_ini_cgs), redshift, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, met_index, d_met, red_index, d_red,
                           cooling, 0, 0, 0, 0);

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  if (LambdaNet_cgs < 0) {

    /* we're cooling! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, 0, 0, 0, 0);

    int i = 0;
    while (u_lower_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs >
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs /= bracket_factor;
      u_upper_cgs /= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          colibre_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, 0, 0, 0, 0);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling \n more info: n_H_cgs = %.4e, u_ini_cgs = %.4e, redshift = "
          "%.4f\n"
          "n_H_index = %i, d_n_H = %.4f\n"
          "met_index = %i, d_met = %.4f, red_index = %i, d_red = %.4f, initial "
          "Lambda = %.4e",
          ID, n_H_cgs, u_ini_cgs, redshift, n_H_index, d_n_H, met_index, d_met,
          red_index, d_red,
          colibre_cooling_rate(log10(u_ini_cgs), redshift, n_H_cgs,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, 0, 0, 0, 0));
    }
  } else {

    /* we are heating! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, 0, 0, 0, 0);

    int i = 0;
    while (u_upper_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs <
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs *= bracket_factor;
      u_upper_cgs *= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          colibre_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, 0, 0, 0, 0);
      i++;
    }

    if (i >= bisection_max_iterations) {
      message("Aborting...");
      message("particle %llu", ID);
      message("n_H_cgs = %.4e", n_H_cgs);
      message("u_ini_cgs = %.4e", u_ini_cgs);
      message("redshift = %.4f", redshift);
      message("indices nH, met, red = %i, %i, %i", n_H_index, met_index,
              red_index);
      message("index weights nH, met, red = %.4e, %.4e, %.4e", d_n_H, d_met,
              d_red);
      fflush(stdout);
      message(
          "cooling rate = %.4e",
          colibre_cooling_rate(log10(u_ini_cgs), redshift, n_H_cgs,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, 0, 0, 0, 0));
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling",
          ID);
    }
  }

  /********************************************/
  /* We now have an upper and lower bound.    */
  /* Let's iterate by reducing the bracketing */
  /********************************************/

  /* bisection iteration */
  int i = 0;
  double u_next_cgs;

  do {

    /* New guess */
    u_next_cgs = 0.5 * (u_lower_cgs + u_upper_cgs);

    /* New rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_next_cgs), redshift, n_H_cgs,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, 0, 0, 0, 0);

    /* Where do we go next? */
    if (u_next_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs > 0.0) {
      u_upper_cgs = u_next_cgs;
    } else {
      u_lower_cgs = u_next_cgs;
    }

    i++;
  } while (fabs(u_upper_cgs - u_lower_cgs) / u_next_cgs > bisection_tolerance &&
           i < bisection_max_iterations);

  if (i >= bisection_max_iterations)
    error("Particle id %llu failed to converge", ID);

  return u_upper_cgs;
}

/**
 * @brief Set the subgrid properties of the gas particle
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param starform the star formation law properties to initialize
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void set_subgrid_part(const struct phys_const *phys_const,
                      const struct unit_system *us,
                      const struct cosmology *cosmo,
                      const struct hydro_props *hydro_props,
                      const struct entropy_floor_properties *floor_props,
                      const struct cooling_function_data *cooling,
                      struct part *restrict p, struct xpart *restrict xp) {

  const double rho = hydro_get_physical_density(p, cosmo);

  /* Get the EOS temperature from the entropy floor */
  const double temperature_eos =
      entropy_floor_temperature(p, cosmo, floor_props);
  const float logT_EOS_max = (float)log10(temperature_eos) + cooling->dlogT_EOS;

  const float temp = cooling_get_temperature(phys_const, hydro_props, us, cosmo,
                                             cooling, p, xp);

  const float logT = log10(temp);

  /* Get internal energy at the last kick step */
  const float u_start = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Get this particle's abundance ratios compared to solar
   * Note that we need to add S and Ca that are in the tables but not tracked
   * by the particles themselves.
   * The order is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, OA] */
  float abundance_ratio[colibre_cooling_N_elementtypes];
  float logZZsol = abundance_ratio_to_solar(p, cooling, abundance_ratio);

  /* Get the Hydrogen and Helium mass fractions */
  float const *metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H =
      hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* Get index along the different table axis */
  int ired, imet, iden;
  float dred, dmet, dden;
  if (cosmo->z < cooling->H_reion_z) {
    get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
                 &ired, &dred);
  } else {
    ired = colibre_cooling_N_redshifts - 2;
    dred = 1.f;
  }
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &imet, &dmet);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &iden,
               &dden);

  if (logT < logT_EOS_max) {

    /* below entropy floor: use subgrid properties */
    const float pres = gas_pressure_from_internal_energy(rho, u_start);
    const double pres_to_cgs =
        units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);

    const float logP = (float)log10((double)pres * pres_to_cgs);

    /* check what would be the maximum Peq from the table for given redshift and
     * metalllicity */
    const float logPeq_max = interpolation_3d_no_z(
        cooling->table.logPeq, ired, imet, colibre_cooling_N_density - 1, dred,
        dmet, 0., colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
        colibre_cooling_N_density);

    float logHI, logHII, logH2;
    float logn_at_Peq = -100.f, logT_at_Peq = -100.f;

    if (logP >= logPeq_max) {
      /* EOS pressure (logP) is larger than maximum Peq (can happen for very
       * steep EOS) use Teq, mu, and HI fractions from the highest density bin,
       * but calculate n from P, T, and mu*/

      logT_at_Peq = interpolation_3d_no_z(
          cooling->table.logTeq, ired, imet, colibre_cooling_N_density - 1,
          dred, dmet, 0., colibre_cooling_N_redshifts,
          colibre_cooling_N_metallicity, colibre_cooling_N_density);

      const float mu_at_Peq = interpolation_3d_no_z(
          cooling->table.meanpartmass_Teq, ired, imet,
          colibre_cooling_N_density - 1, dred, dmet, 0.,
          colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
          colibre_cooling_N_density);

      /*const float logkB = (float) log10(const_boltzmann_k_cgs); doesn't work*/
      const double log10_kB = cooling->log10_kB_cgs;

      logn_at_Peq =
          logP - logT_at_Peq + log10(XH) + log10(mu_at_Peq) - log10_kB;

      logHI = interpolation_4d_no_z_no_w(
          cooling->table.logHfracs_Teq, ired, imet,
          colibre_cooling_N_density - 1, neutral, dred, dmet, 0., 0.,
          colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
          colibre_cooling_N_density, 3);

      logHII = interpolation_4d_no_z_no_w(
          cooling->table.logHfracs_Teq, ired, imet,
          colibre_cooling_N_density - 1, ionized, dred, dmet, 0., 0.,
          colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
          colibre_cooling_N_density, 3);

      logH2 = interpolation_4d_no_z_no_w(
          cooling->table.logHfracs_Teq, ired, imet,
          colibre_cooling_N_density - 1, molecular, dred, dmet, 0., 0.,
          colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
          colibre_cooling_N_density, 3);

    } else {

      /* need to find thermal equilibrium state with the same pressure *
       * logPeq is neither equally spaced, nor necessarily increasing
       * monotonically * simple solution: loop over densities and pick the first
       * one where logP = logPeq * start with the resolved density index (iden),
       * as the subgrid density will always be higher than the resolved density
       *
       * In cases where the solution cannot be found (e.g. the equilibrium
       * temperature solution intersects the entropy floor), we revert to the
       * normal (non-sugrid) density.
       */

      int iden_eq = iden;
      float dden_eq = dden;

      for (int i = iden; i < colibre_cooling_N_density; i++) {
        float logPeq_interp = interpolation_3d_no_z(
            cooling->table.logPeq, ired, imet, i, dred, dmet, 0.,
            colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
            colibre_cooling_N_density);
        if (logPeq_interp > logP) {
          const float logPeq_prev = interpolation_3d_no_z(
              cooling->table.logPeq, ired, imet, i - 1, dred, dmet, 0.,
              colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
              colibre_cooling_N_density);

          logn_at_Peq = (logP - logPeq_prev) / (logPeq_interp - logPeq_prev) *
                            (cooling->nH[i] - cooling->nH[i - 1]) +
                        cooling->nH[i - 1];

          iden_eq = i - 1;
          dden_eq = (logn_at_Peq - cooling->nH[i - 1]) /
                    (cooling->nH[i] - cooling->nH[i - 1]);
          break;
        }
      }

      logHI = interpolation_4d_no_w(
          cooling->table.logHfracs_Teq, ired, imet, iden_eq, neutral, dred,
          dmet, dden_eq, 0., colibre_cooling_N_redshifts,
          colibre_cooling_N_metallicity, colibre_cooling_N_density, 3);

      logHII = interpolation_4d_no_w(
          cooling->table.logHfracs_Teq, ired, imet, iden_eq, ionized, dred,
          dmet, dden_eq, 0., colibre_cooling_N_redshifts,
          colibre_cooling_N_metallicity, colibre_cooling_N_density, 3);

      logH2 = interpolation_4d_no_w(
          cooling->table.logHfracs_Teq, ired, imet, iden_eq, molecular, dred,
          dmet, dden_eq, 0., colibre_cooling_N_redshifts,
          colibre_cooling_N_metallicity, colibre_cooling_N_density, 3);

      logT_at_Peq = interpolation_3d(
          cooling->table.logTeq, ired, imet, iden_eq, dred, dmet, dden_eq,
          colibre_cooling_N_redshifts, colibre_cooling_N_metallicity,
          colibre_cooling_N_density);
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (logn_at_Peq == -100.f) error("Did not find a solution!");
    if (logT_at_Peq == -100.f) error("Did not find a solution!");
#endif

    xp->tracers_data.nHI_over_nH = exp10(logHI);
    xp->tracers_data.nHII_over_nH = exp10(logHII);
    xp->tracers_data.nH2_over_nH = 0.5 * exp10(logH2);
    xp->tracers_data.subgrid_temp = exp10(logT_at_Peq);
    /* convert log nH to comoving density in SU */
    xp->tracers_data.subgrid_dens =
        exp10(logn_at_Peq) / cooling->number_density_to_cgs / XH *
        phys_const->const_proton_mass / cosmo->a3_inv;

  } else {

    /* above entropy floor: use table properties */
    const double crho = hydro_get_comoving_density(p);
    /* subgrid_dens should be the same as p->rho */
    xp->tracers_data.subgrid_dens = crho;

    /* interpolate the tables for H fractions */
    /* check if in an HII region */
    if (xp->tracers_data.HIIregion_timer_gas > 0.) {
      xp->tracers_data.subgrid_temp = cooling->HIIregion_temp;
      xp->tracers_data.nHI_over_nH = 1. - cooling->HIIregion_fion;
      xp->tracers_data.nHII_over_nH = cooling->HIIregion_fion;
      xp->tracers_data.nH2_over_nH = 0.;
    } else {
      xp->tracers_data.subgrid_temp = temp;

      int item;
      float dtem;

      get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log10(temp),
                   &item, &dtem);

      /* necessary to define to use the interpolation routine */
      const float weights[3] = {1.0, 1.0, 1.0};

      /* "sum" from element 0 to 0 (neutral hydrogen) */
      xp->tracers_data.nHI_over_nH = interpolation4d_plus_summation(
          cooling->table.logHfracs_all, weights, neutral, neutral, ired, item,
          imet, iden, dred, dtem, dmet, dden, colibre_cooling_N_redshifts,
          colibre_cooling_N_temperature, colibre_cooling_N_metallicity,
          colibre_cooling_N_density, 3);

      /* "sum" from element 1 to 1 (ionized hydrogen) */
      xp->tracers_data.nHII_over_nH = interpolation4d_plus_summation(
          cooling->table.logHfracs_all, weights, ionized, ionized, ired, item,
          imet, iden, dred, dtem, dmet, dden, colibre_cooling_N_redshifts,
          colibre_cooling_N_temperature, colibre_cooling_N_metallicity,
          colibre_cooling_N_density, 3);

      /* "sum" from element 2 to 2 (molecular hydrogen) */
      xp->tracers_data.nH2_over_nH =
          0.5 * interpolation4d_plus_summation(
                    cooling->table.logHfracs_all, weights, molecular, molecular,
                    ired, item, imet, iden, dred, dtem, dmet, dden,
                    colibre_cooling_N_redshifts, colibre_cooling_N_temperature,
                    colibre_cooling_N_metallicity, colibre_cooling_N_density,
                    3);
    }
  }
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * We want to compute u_new such that u_new = u_old + dt * du/dt(u_new, X),
 * where X stands for the metallicity, density and redshift. These are
 * kept constant.
 *
 * We first compute du/dt(u_old). If dt * du/dt(u_old) is small enough, we
 * use an explicit integration and use this as our solution.
 *
 * Otherwise, we try to find a solution to the implicit time-integration
 * problem. This leads to the root-finding problem:
 *
 * f(u_new) = u_new - u_old - dt * du/dt(u_new, X) = 0
 *
 * A bisection scheme is used.
 * This is done by first bracketing the solution and then iterating
 * towards the solution by reducing the window down to a certain tolerance.
 * Note there is always at least one solution since
 * f(+inf) is < 0 and f(-inf) is > 0.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_properties the hydro_props struct
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The cooling time-step of this particle.
 * @param dt_therm The hydro time-step of this particle.
 * @param time Time since Big Bang
 */
void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm,
                       const double time) {

  /* No cooling happens over zero time */
  if (dt == 0.) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--cooling runtime flag?");
#endif

  /* Get internal energy at the last kick step */
  const float u_start = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Get the change in internal energy due to hydro forces */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Get internal energy at the end of the next kick step (assuming dt does not
   * increase) */
  double u_0 = (u_start + hydro_du_dt * dt_therm);

  /* Check for minimal energy */
  u_0 = max(u_0, hydro_properties->minimal_internal_energy);

  /* Convert to CGS units */
  const double u_0_cgs = u_0 * cooling->internal_energy_to_cgs;
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Change in redshift over the course of this time-step
     (See cosmology theory document for the derivation) */
  const double delta_redshift = -dt * cosmo->H * cosmo->a_inv;

  /* Get this particle's abundance ratios compared to solar
   * Note that we need to add S and Ca that are in the tables but not tracked
   * by the particles themselves.
   * The order is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, OA] */
  float abundance_ratio[colibre_cooling_N_elementtypes];
  float logZZsol = abundance_ratio_to_solar(p, cooling, abundance_ratio);

  /* Get the Hydrogen and Helium mass fractions */
  float const *metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H =
      hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* ratefact = n_H * n_H / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  const double ratefact_cgs = n_H_cgs * (XH * cooling->inv_proton_mass_cgs);

  /* compute hydrogen number density, metallicity and redshift indices and
   * offsets (These are fixed for any value of u, so no need to recompute them)
   */

  float d_red, d_met, d_n_H;
  int red_index, met_index, n_H_index;

  if (cosmo->z < cooling->H_reion_z) {
    get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
                 &red_index, &d_red);
  } else {
    red_index = colibre_cooling_N_redshifts - 2;
    d_red = 1.0;
  }
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Start by computing the cooling (heating actually) rate from Helium
     re-ionization as this needs to be added on no matter what */

  /* Get helium and hydrogen reheating term */
  const double Helium_reion_heat_cgs =
      eagle_helium_reionization_extraheat(cosmo->z, delta_redshift, cooling);

  /* Convert this into a rate */
  const double Lambda_He_reion_cgs =
      Helium_reion_heat_cgs / (dt_cgs * ratefact_cgs);

  /* Let's compute the internal energy at the end of the step */
  double u_final_cgs;

  /* First try an explicit integration (note we ignore the derivative) */
  const double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      colibre_cooling_rate(log10(u_0_cgs), cosmo->z, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, met_index, d_met, red_index, d_red,
                           cooling, 0, 0, 0, 0);

  /* if cooling rate is small, take the explicit solution */
  if (fabs(ratefact_cgs * LambdaNet_cgs * dt_cgs) <
      explicit_tolerance * u_0_cgs) {

    u_final_cgs = u_0_cgs + ratefact_cgs * LambdaNet_cgs * dt_cgs;

  } else {

    u_final_cgs =
        bisection_iter(u_0_cgs, n_H_cgs, cosmo->z, n_H_index, d_n_H, met_index,
                       d_met, red_index, d_red, Lambda_He_reion_cgs,
                       ratefact_cgs, cooling, abundance_ratio, dt_cgs, p->id);
  }

  /* Convert back to internal units */
  double u_final = u_final_cgs * cooling->internal_energy_from_cgs;

  /* We now need to check that we are not going to go below any of the limits */

  /* Absolute minimum */
  const double u_minimal = hydro_properties->minimal_internal_energy;
  u_final = max(u_final, u_minimal);

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho_physical = hydro_get_physical_density(p, cosmo);
  const double u_floor =
      gas_internal_energy_from_entropy(rho_physical, A_floor);
  u_final = max(u_final, u_floor);

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u = u_final - max(u_start, u_floor);

  /* Determine if we are in the slow- or rapid-cooling regime,
   * by comparing dt / t_cool to the rapid_cooling_threshold.
   * Note that dt / t_cool = fabs(delta_u) / u_start.
   * If rapid_cooling_threshold < 0, always use the slow-cooling
   * regime. */
  float cooling_du_dt;
  if ((cooling->rapid_cooling_threshold >= 0.0) &&
      (fabs(delta_u) / max(u_start, u_floor) >=
       cooling->rapid_cooling_threshold)) {
    /* Rapid-cooling regime. Update internal energy
     * to u_final and set du/dt = 0. */
    cooling_du_dt = 0.0;

    /* Update the particle's u and du/dt */
    hydro_set_physical_internal_energy(p, xp, cosmo, u_final);
    hydro_set_drifted_physical_internal_energy(p, cosmo, u_final);
    hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);

    /* Store the radiated energy */
    xp->cooling_data.radiated_energy -=
        hydro_get_mass(p) * (u_final - u_0);
  } else {
    /* Slow-cooling regime. Update du/dt so that
     * we can subsequently drift internal energy. */
    cooling_du_dt = delta_u / dt_therm;

    /* Update the internal energy time derivative */
    hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);

    /* Store the radiated energy */
    xp->cooling_data.radiated_energy -= hydro_get_mass(p) * dt * ((u_final - u_0) / dt_therm);
  }

  /* check if the particle is in an HII region and if yes, set the parameter
   * accordingly */
  if ((time <= xp->tracers_data.HIIregion_timer_gas) &&
      (xp->tracers_data.HIIregion_timer_gas > 0.)) {
    /*const float temp = cooling_get_temperature (phys_const, hydro_properties,
     * us, cosmo, cooling, p, xp); */

    const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
    /* HII region internal energy is the internal energy of a particle at a
     * temperature of cooling->HIIregion_temp */
    const float u_HII_cgs = cooling_get_internalenergy_for_temperature(
        phys_const, hydro_properties, us, cosmo, cooling, p, xp,
        cooling->HIIregion_temp);

    const float u_HII = u_HII_cgs / cooling->internal_energy_to_cgs;

    if (u_old < u_HII) {
      /* Inject energy into the particle */
      hydro_set_physical_internal_energy(p, xp, cosmo, u_HII);
      hydro_set_drifted_physical_internal_energy(p, cosmo, u_HII);

      /* internal energy should stay constant for the timestep */
      const float cooling_du_dt_HII = 0.;
      hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt_HII);
    }
  } else if ((time > xp->tracers_data.HIIregion_timer_gas) &&
             (xp->tracers_data.HIIregion_timer_gas > 0.)) {
    xp->tracers_data.HIIregion_timer_gas = -1.;
    xp->tracers_data.HIIregion_starid = -1;
  }

  /* set subgrid properties and hydrogen fractions */
  set_subgrid_part(phys_const, us, cosmo, hydro_properties, floor_props,
                   cooling, p, xp);
}

/**
 * @brief Computes the cooling time-step.
 *
 * The time-step is not set by the properties of cooling.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const #phys_const data struct.
 * @param us The internal system of units.
 * @param cosmo #cosmology struct.
 * @param hydro_props the properties of the hydro scheme.
 * @param p #part data.
 * @param xp extended particle data.
 */
__attribute__((always_inline)) INLINE float cooling_timestep(
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props, const struct part *restrict p,
    const struct xpart *restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const #phys_const data structure.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
__attribute__((always_inline)) INLINE void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
}

/**
 * @brief Compute the internal energy of a #part based on the cooling function
 * but for a given temperature. This is used e.g. for particles in HII regions
 * that are set to a constant temperature, but their internal energies should
 * reflect the particle composition
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @param desired gas temperature
 */
float cooling_get_internalenergy_for_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp,
    float temp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--temperature runtime flag?");
#endif

  /* Get the Hydrogen mass fraction */
  float const *metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float rho = hydro_get_physical_density(p, cosmo);
  const double n_H = rho * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  const float logZZsol =
      log10(chemistry_get_total_metal_mass_fraction_for_cooling(p) /
            cooling->Zsol[0]);

  /* compute hydrogen number density, metallicity and redshift indices and
   * offsets  */

  float d_red, d_met, d_n_H;
  int red_index, met_index, n_H_index;

  if (cosmo->z < cooling->H_reion_z) {
    get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
                 &red_index, &d_red);
  } else {
    red_index = colibre_cooling_N_redshifts - 2;
    d_red = 1.0;
  }

  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Compute the log10 of the temperature by interpolating the table */
  const double log_10_U =
      colibre_convert_temp_to_u(log10(temp), cosmo->z, n_H_index, d_n_H,
                                met_index, d_met, red_index, d_red, cooling);

  /* Undo the log! */
  return exp10(log_10_U);
}

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * The temperature returned is consistent with the cooling rates.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--temperature runtime flag?");
#endif

  /* Get physical internal energy */
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const double u_cgs = u * cooling->internal_energy_to_cgs;

  /* Get the Hydrogen mass fraction */
  float const *metal_fraction =
      chemistry_get_metal_mass_fraction_for_cooling(p);
  const float XH = metal_fraction[chemistry_element_H];

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float rho = hydro_get_physical_density(p, cosmo);
  const double n_H = rho * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  const float logZZsol =
      log10(chemistry_get_total_metal_mass_fraction_for_cooling(p) /
            cooling->Zsol[0]);

  /* compute hydrogen number density, metallicity and redshift indices and
   * offsets  */

  float d_red, d_met, d_n_H;
  int red_index, met_index, n_H_index;

  if (cosmo->z < cooling->H_reion_z) {
    get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
                 &red_index, &d_red);
  } else {
    red_index = colibre_cooling_N_redshifts - 2;
    d_red = 1.0;
  }

  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Compute the log10 of the temperature by interpolating the table */
  const double log_10_T =
      colibre_convert_u_to_temp(log10(u_cgs), cosmo->z, n_H_index, d_n_H,
                                met_index, d_met, red_index, d_red, cooling);

  /* Undo the log! */
  return exp10(log_10_T);
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp #xpart data struct
 */
__attribute__((always_inline)) INLINE float cooling_get_radiated_energy(
    const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Inject a fixed amount of energy to each particle in the simulation
 * to mimic Hydrogen reionization.
 *
 * @param cooling The properties of the cooling model.
 * @param cosmo The cosmological model.
 * @param s The #space containing the particles.
 */
void cooling_Hydrogen_reionization(const struct cooling_function_data *cooling,
                                   const struct cosmology *cosmo,
                                   struct space *s) {

  struct part *parts = s->parts;
  struct xpart *xparts = s->xparts;

  /* Energy to inject in internal units */
  const float extra_heat =
      cooling->H_reion_heat_cgs * cooling->internal_energy_from_cgs;

  message("Applying extra energy for H reionization!");

  /* Loop through particles and set new heat */
  for (size_t i = 0; i < s->nr_parts; i++) {

    struct part *p = &parts[i];
    struct xpart *xp = &xparts[i];

    const float old_u = hydro_get_physical_internal_energy(p, xp, cosmo);
    const float new_u = old_u + extra_heat;

    hydro_set_physical_internal_energy(p, xp, cosmo, new_u);
    hydro_set_drifted_physical_internal_energy(p, cosmo, new_u);
  }
}

/**
 * @brief Initialises properties stored in the cooling_function_data struct
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const #phys_const data structure
 * @param cooling #cooling_function_data struct to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          struct cooling_function_data *cooling) {

  /* read some parameters */

  /* Despite the names, the values of H_reion_heat_cgs and He_reion_heat_cgs
   * that are read in are actually in units of electron volts per proton mass.
   * We later convert to units just below */

  parser_get_param_string(parameter_file, "COLIBRECooling:dir_name",
                          cooling->cooling_table_path);

  cooling->H_reion_done = 0;
  cooling->H_reion_z =
      parser_get_param_float(parameter_file, "COLIBRECooling:H_reion_z");
  cooling->H_reion_heat_cgs =
      parser_get_param_float(parameter_file, "COLIBRECooling:H_reion_eV_p_H");
  cooling->He_reion_z_centre = parser_get_param_float(
      parameter_file, "COLIBRECooling:He_reion_z_centre");
  cooling->He_reion_z_sigma =
      parser_get_param_float(parameter_file, "COLIBRECooling:He_reion_z_sigma");
  cooling->He_reion_heat_cgs =
      parser_get_param_float(parameter_file, "COLIBRECooling:He_reion_eV_p_H");

  /* Properties of the HII region model ------------------------------------- */
  cooling->HIIregion_fion = parser_get_param_float(
      parameter_file, "COLIBRECooling:HIIregion_ionization_fraction");

  cooling->HIIregion_temp = parser_get_param_float(
      parameter_file, "COLIBRECooling:HIIregion_temperature");

  /* Properties for the subgrid properties model */
  cooling->dlogT_EOS = parser_get_param_float(
      parameter_file, "COLIBRECooling:delta_logTEOS_subgrid_properties");

  /* Check that it makes sense. */
  if (cooling->HIIregion_fion < 0.5 || cooling->HIIregion_fion > 1.0) {
    error("HIIregion_ionization_fraction has to be between 0.5 and 1.0");
  }

  /* Optional parameters to correct the abundances */
  cooling->Ca_over_Si_ratio_in_solar = parser_get_opt_param_float(
      parameter_file, "COLIBRECooling:Ca_over_Si_in_solar", 1.f);
  cooling->S_over_Si_ratio_in_solar = parser_get_opt_param_float(
      parameter_file, "COLIBRECooling:S_over_Si_in_solar", 1.f);

  /* Convert H_reion_heat_cgs and He_reion_heat_cgs to cgs
   * (units used internally by the cooling routines). This is done by
   * multiplying by 'eV/m_H' in internal units, then converting to cgs units.
   * Note that the dimensions of these quantities are energy/mass = velocity^2
   */

  cooling->H_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  cooling->He_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Compute conversion factors */
  cooling->internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  cooling->internal_energy_from_cgs = 1. / cooling->internal_energy_to_cgs;
  cooling->number_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Store some constants in CGS units */
  const float units_kB[5] = {1, 2, -2, 0, -1};
  const double kB_cgs = phys_const->const_boltzmann_k *
                        units_general_cgs_conversion_factor(us, units_kB);
  const double proton_mass_cgs =
      phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);

  cooling->log10_kB_cgs = log10(kB_cgs);
  cooling->inv_proton_mass_cgs = 1. / proton_mass_cgs;
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

#ifdef SWIFT_DEBUG_CHECKS
  /* Basic cross-check... */
  if (kB_cgs > 1.381e-16 || kB_cgs < 1.380e-16)
    error("Boltzmann's constant not initialised properly!");
#endif

  /* Compute the coefficient at the front of the Compton cooling expression */
  const double radiation_constant =
      4. * phys_const->const_stefan_boltzmann / phys_const->const_speed_light_c;
  const double compton_coefficient =
      4. * radiation_constant * phys_const->const_thomson_cross_section *
      phys_const->const_boltzmann_k /
      (phys_const->const_electron_mass * phys_const->const_speed_light_c);
  const float dimension_coefficient[5] = {1, 2, -3, 0, -5};

  /* This should be ~1.0178085e-37 g cm^2 s^-3 K^-5 */
  const double compton_coefficient_cgs =
      compton_coefficient *
      units_general_cgs_conversion_factor(us, dimension_coefficient);

#ifdef SWIFT_DEBUG_CHECKS
  const double expected_compton_coefficient_cgs = 1.0178085e-37;
  if (fabs(compton_coefficient_cgs - expected_compton_coefficient_cgs) /
          expected_compton_coefficient_cgs >
      0.01)
    error("compton coefficient incorrect.");
#endif

  /* And now the Compton rate */
  cooling->compton_rate_cgs = compton_coefficient_cgs * cooling->T_CMB_0 *
                              cooling->T_CMB_0 * cooling->T_CMB_0 *
                              cooling->T_CMB_0;

  /* Threshold in dt / t_cool above which we
   * are in the rapid cooling regime. If negative,
   * we never use this scheme (i.e. always drift
   * the internal energies). */
  cooling->rapid_cooling_threshold = parser_get_param_double(
      parameter_file, "COLIBRECooling:rapid_cooling_threshold");

  /* Finally, read the tables */
  read_cooling_header(cooling);
  read_cooling_tables(cooling);
}

/**
 * @brief Restore cooling tables (if applicable) after
 * restart
 *
 * @param cooling the #cooling_function_data structure
 * @param cosmo #cosmology structure
 */
void cooling_restore_tables(struct cooling_function_data *cooling,
                            const struct cosmology *cosmo) {

  read_cooling_header(cooling);
  read_cooling_tables(cooling);

  cooling_update(cosmo, cooling, /*space=*/NULL);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling #cooling_function_data struct.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  message("Cooling function is 'COLIBRE'.");
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * We simply free all the arrays.
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data *cooling) {

  /* Free the side arrays */
  free(cooling->Redshifts);
  free(cooling->nH);
  free(cooling->Temp);
  free(cooling->Metallicity);
  free(cooling->Therm);
  free(cooling->LogAbundances);
  free(cooling->Abundances);
  free(cooling->Abundances_inv);
  free(cooling->atomicmass);
  free(cooling->LogMassFractions);
  free(cooling->MassFractions);

  /* Free the tables */
  free(cooling->table.Tcooling);
  free(cooling->table.Ucooling);
  free(cooling->table.Theating);
  free(cooling->table.Uheating);
  free(cooling->table.Telectron_fraction);
  free(cooling->table.Uelectron_fraction);
  free(cooling->table.T_from_U);
  free(cooling->table.U_from_T);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the cooling routines. Helps debugging. */
  struct cooling_function_data cooling_copy = *cooling;
  cooling_copy.Redshifts = NULL;
  cooling_copy.nH = NULL;
  cooling_copy.Temp = NULL;
  cooling_copy.Metallicity = NULL;
  cooling_copy.Therm = NULL;
  cooling_copy.LogAbundances = NULL;
  cooling_copy.Abundances = NULL;
  cooling_copy.Abundances_inv = NULL;
  cooling_copy.atomicmass = NULL;
  cooling_copy.LogMassFractions = NULL;
  cooling_copy.MassFractions = NULL;

  cooling_copy.table.Tcooling = NULL;
  cooling_copy.table.Theating = NULL;
  cooling_copy.table.Telectron_fraction = NULL;
  cooling_copy.table.Ucooling = NULL;
  cooling_copy.table.Uheating = NULL;
  cooling_copy.table.Uelectron_fraction = NULL;
  cooling_copy.table.T_from_U = NULL;
  cooling_copy.table.U_from_T = NULL;

  restart_write_blocks((void *)&cooling_copy,
                       sizeof(struct cooling_function_data), 1, stream,
                       "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the cooling tables by
 * re-reading them.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo) {
  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  cooling_restore_tables(cooling, cosmo);
}
