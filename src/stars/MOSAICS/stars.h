/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2019 Joel Pfeffer (j.l.pfeffer@ljmu.ac.uk)
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
#ifndef SWIFT_MOSAICS_STARS_H
#define SWIFT_MOSAICS_STARS_H

#include <float.h>

/* Local includes */
#include "cosmology.h"
#include "hydro.h"
#include "cooling.h"
#include "star_formation.h"
#include "physical_constants.h"
#include "mosaics_clevo.h"
#include "mosaics_clform.h"

/**
 * @brief Computes the gravity time-step of a given star particle.
 *
 * @param sp Pointer to the s-particle data.
 */
__attribute__((always_inline)) INLINE static float stars_compute_timestep(
    const struct spart* const sp) {

  return FLT_MAX;
}

/**
 * @brief Prepares a s-particle for its interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_init_spart(
    struct spart* sp) {

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    sp->ids_ngbs_density[i] = -1;
  sp->num_ngb_density = 0;
#endif

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;

  if (sp->gpart) sp->gpart->calc_tensor = sp->calc_tensor;
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param stars_properties Properties of the stars model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param scale_factor The current scale-factor (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp, const struct stars_props* stars_properties,
    const int with_cosmology, const double scale_factor, const double time) {

  sp->time_bin = 0;
  sp->birth_density = 0.f;
  sp->SNII_f_E = -1.f;
  sp->count_since_last_enrichment = -1;

  if (stars_properties->overwrite_birth_time)
    sp->birth_time = stars_properties->spart_first_init_birth_time;

  if (with_cosmology)
    sp->last_enrichment_time = scale_factor;
  else
    sp->last_enrichment_time = time;

  sp->HIIregion_last_rebuild = -1.f;
  sp->HIIregion_mass_to_ionize = 0.f;
  sp->HIIregion_mass_in_kernel = -1.f;
  sp->star_timestep = 0.f;

  /* TODO temporary until we read gc props from ICs */
  /* i.e. does this particle have any clusters with M>0? */
  sp->gcflag = 0;

  stars_init_spart(sp);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void stars_predict_extra(
    struct spart* restrict sp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The particle.
 */
__attribute__((always_inline)) INLINE static void stars_reset_predicted_values(
    struct spart* restrict sp) {}

/**
 * @brief Finishes the calculation of (non-gravity) forces acting on stars
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_end_feedback(
    struct spart* sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void stars_kick_extra(
    struct spart* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_density(
    struct spart* sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->density.wcount *= h_inv_dim;
  sp->density.wcount_dh *= h_inv_dim_plus_one;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_spart_has_no_neighbours(
    struct spart* restrict sp, const struct cosmology* cosmo) {

  /* Re-set problematic values */
  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_acceleration(
    struct spart* restrict p) {

#ifdef DEBUG_INTERACTIONS_STARS
  p->num_ngb_force = 0;
#endif
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_feedback(
    struct spart* restrict p) {

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    p->ids_ngbs_force[i] = -1;
  p->num_ngb_force = 0;
#endif
}

/**
 * @brief Do the mosaics subgrid star cluster model
 *
 * @param sp The particle to act upon
 * @param stars_properties Properties of the stars model.
 * @param sf_props the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology Are we running with cosmological time integration.
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static void stars_do_mosaics(
    struct spart* restrict sp, const struct stars_props* stars_properties,
    const struct star_formation* sf_props, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const int with_cosmology, const float time) {

  /* Have we already been here this timestep? */
  /* This funcation can be called twice: once at SF, then again in star loop */
  /* Now fixed so not necessary
  if (!sp->new_star) {
    if (with_cosmology) {
      if (sp->birth_scale_factor == (float)cosmo->a) {
        return;
      }
    } else {
      if (sp->birth_time == time) {
        return;
      }
    }
  }
  */

  /* shift old tensors along, regardless if have clusters */
  for (int i = 0; i < 2; i++) {
    sp->tidal_tensor[i][0] = sp->tidal_tensor[i + 1][0];
    sp->tidal_tensor[i][1] = sp->tidal_tensor[i + 1][1];
    sp->tidal_tensor[i][2] = sp->tidal_tensor[i + 1][2];
    sp->tidal_tensor[i][3] = sp->tidal_tensor[i + 1][3];
    sp->tidal_tensor[i][4] = sp->tidal_tensor[i + 1][4];
    sp->tidal_tensor[i][5] = sp->tidal_tensor[i + 1][5];
  }

  /* Now retrieve the new ones from the gpart */
  sp->tidal_tensor[2][0] = sp->gpart->tidal_tensor[0];
  sp->tidal_tensor[2][1] = sp->gpart->tidal_tensor[1];
  sp->tidal_tensor[2][2] = sp->gpart->tidal_tensor[2];
  sp->tidal_tensor[2][3] = sp->gpart->tidal_tensor[3];
  sp->tidal_tensor[2][4] = sp->gpart->tidal_tensor[4];
  sp->tidal_tensor[2][5] = sp->gpart->tidal_tensor[5];

  /* TODO just temporary */
  sp->potential = sp->gpart->potential;

  /* Formation or evolution? */
  if (sp->new_star) {
    /* Do cluster formation */
    sp->gcflag = 1;

    /* Go make clusters */
    mosaics_clform(sp, stars_properties, sf_props, phys_const, cosmo);

    /* We're done with cluster formation */
    sp->new_star = 0;

  } else if (sp->gcflag) {
    /* Do cluster evolution, if the particle has clusters */
    mosaics_clevo(sp, stars_properties);
  }

}

/**
 * @brief Gather some extra props needed at star formation time for mosaics
 *
 * @param p The gas particle
 * @param xp The additional properties of the gas particles.
 * @param sp The star particle
 * @param cosmo the cosmological parameters and properties.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * cosmology).
 */
__attribute__((always_inline)) INLINE static void
stars_mosaics_copy_extra_properties(
    const struct part* p, const struct xpart* xp, struct spart* restrict sp,
    const struct cosmology* cosmo, const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* Flag it for cluster formation */
  sp->new_star = 1;

  /* Store the birth properties in the star particle */
  sp->hbirth = p->h;
  sp->gas_vel_disp = sqrt(p->sf_data.sigma_v2);
  sp->gas_mass_unweighted = p->sf_data.gas_mass_unweighted;

  /* Hydro pressure */
  sp->birth_pressure = hydro_get_physical_pressure(p, cosmo);

  /* Calculate the hydro sound speed */
  sp->sound_speed_subgrid = hydro_get_physical_soundspeed(p, cosmo);

  /* Set up for the stellar neighbour search */
  sp->scount = 0;
  sp->stars_rho = 0.f;
  sp->stars_sigma_v2 = 0.f;
  sp->stars_mass_unweighted = 0.f;

#if defined(COOLING_COLIBRE) || defined(COOLING_CHIMES) || \
    defined(COOLING_CHIMES_HYBRID)

  /* Correction for subgrid ISM */

  /* Calculate the temperature */
  const double temperature = cooling_get_temperature(phys_const, hydro_props,
                                                     us, cosmo, cooling, p, xp);

  /* Get the subgrid temperature from the tracers */
  const double subgrid_temperature = xp->tracers_data.subgrid_temp;

  /* Get the subgrid sound speed */
  sp->sound_speed_subgrid *= sqrt(subgrid_temperature / temperature);

  /* Get the subgrid density */
  sp->birth_subgrid_dens = xp->tracers_data.subgrid_dens;

#else

  /* Get the hydro density */
  sp->birth_subgrid_dens = sp->birth_density;

#endif
}

#endif /* SWIFT_MOSAICS_STARS_H */
