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
#include "engine.h"
#include "cosmology.h"
#include "hydro.h"
#include "cooling.h"
#include "star_formation.h"
#include "physical_constants.h"
#include "mosaics_clevo.h"
#include "mosaics_clform.h"
#include "dsyevj3.h"

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
#if !defined(STAR_FORMATION_NONE)
  sp->sf_data.birth_density = 0.f;
#endif
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

  /* Only MOSAICS task for this particle is stellar evo. for field */
  sp->num_clusters = 0;
  sp->initial_num_clusters = 0;
  sp->initial_num_clusters_evo = 0;
  sp->field_mass = sp->mass;

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
 * @param e The #engine.
 * @param with_cosmology Are we running with cosmological time integration.
 */
__attribute__((always_inline)) INLINE static void stars_do_mosaics(
    struct spart* restrict sp, const struct engine* e,
    const int with_cosmology) {

  /* Shift the old tensors along */
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

    /* For applying stellar evolutionary mass loss */
    sp->mass_prev_timestep = sp->mass;

    /* Go make clusters */
    mosaics_clform(sp, e, with_cosmology);

    /* We're done with cluster formation */
    sp->new_star = 0;

  } else {

    /* Do cluster evolution, and apply stellar evo. to field props */
    mosaics_clevo(sp, e, with_cosmology);
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

#if !defined(STAR_FORMATION_NONE)

  /* Flag it for cluster formation and the stellar neighbour search */
  sp->new_star = 1;

  /* Set up for the stellar neighbour search */
  sp->scount = 0;
  sp->stars_rho = 0.f;
  sp->stars_sigma_v2 = 0.f;
  sp->stars_mass_unweighted = 0.f;
  sp->Omega_birth = 0.f;
  sp->kappa_birth = 0.f;

  /* Store the birth properties in the star particle */
  sp->hbirth = p->h;
  sp->gas_mass_unweighted = p->sf_data.gas_mass_unweighted;
  sp->birth_toomre_length = p->sf_data.toomre_length * cosmo->a;

  /* Hydro pressure */
  sp->birth_pressure = hydro_get_physical_pressure(p, cosmo);

  /* Calculate the hydro sound speed */
  sp->sound_speed_subgrid = hydro_get_physical_soundspeed(p, cosmo);

#if defined(COOLING_COLIBRE) || defined(COOLING_CHIMES) || \
    defined(COOLING_CHIMES_HYBRID)

  /* Correction for subgrid ISM */
  sp->sound_speed_subgrid *= 
      sqrt(sp->sf_data.birth_subgrid_temperature / sp->sf_data.birth_temperature);

#endif /* cooling model */

#endif /* !defined(STAR_FORMATION_NONE) */
}

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param e The #engine
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void mosaics_calc_toomre_length(
    struct part* p, const struct engine* e, const struct cosmology* cosmo) {

  const struct stars_props* star_props = e->stars_properties;
  const double const_G = e->physical_constants->const_newton_G;

  /* Did it already become a star? */
  if (!p->gpart)
    return;

  /* Temporary array for eigvec/val calculation
   * Note the lower elements are never referenced by dsyevj3 */
  double tensor[3][3];
  tensor[0][0] = p->gpart->tidal_tensor[0];
  tensor[0][1] = p->gpart->tidal_tensor[1];
  tensor[0][2] = p->gpart->tidal_tensor[2];
  tensor[1][1] = p->gpart->tidal_tensor[3];
  tensor[1][2] = p->gpart->tidal_tensor[4];
  tensor[2][2] = p->gpart->tidal_tensor[5];

  double sum_tensor = 0.f;
  for (int i = 0; i < 6; i++)
    sum_tensor += p->gpart->tidal_tensor[i];
  /*TODO I think gravity is calculated by this point, so must be something
   * going on? */
  //if (sum_tensor == 0)
  //  error("Tidal tensor hasn't been calculated yet for this particle!");

  /* Can't do anything in this case, so keep the length we have */
  if (sum_tensor == 0)
    return;

  /* Get the eigenvalues */
  double tideval[3];
  dsyevj3(tensor, tideval);

  /* sort eigenvalues, tideval[2] is largest eigenvalue */
  sort3(tideval);

  /* Gas properties */
  const double rho = hydro_get_physical_density(p, cosmo);
  const double pressure = rho * p->sf_data.sigma_v2 / 3.f;

  /* Gas surface density, assuming phi_P = 1 */
  const double SigmaG = sqrt(2. * pressure / (M_PI * const_G));

  /* Circular and epicyclic frequencies */
  double Omega2;
  if (star_props->Omega_is_lambda2) {
    /* Correct version (in principle) */
    Omega2 = fabs(-tideval[1]);
  } else {
    /* Version in E-MOSAICS */
    Omega2 = fabs(-tideval[0] - tideval[1] - tideval[2]) / 3.f;
  }
  double kappa2 = fabs(3.f * Omega2 - tideval[2]);
  kappa2 *= cosmo->a3_inv;

  p->sf_data.toomre_length = 4.f * M_PI * M_PI * const_G * SigmaG / kappa2;

  /* Maximum length */
  if (p->sf_data.toomre_length > star_props->max_toomre_length)
    p->sf_data.toomre_length = star_props->max_toomre_length;

  /* Back in comoving units */
  p->sf_data.toomre_length /= cosmo->a;

  /* Copy to gpart for next gravity calculation */
  p->gpart->htoomre = p->sf_data.toomre_length;
}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 *
 * @param p Pointer to the particle data.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 */
__attribute__((always_inline)) INLINE static void
mosaics_first_init_part(struct part* restrict p, 
                        const struct gravity_props* grav_props) {

  p->sf_data.toomre_length = grav_props->epsilon_baryon_cur * 
      grav_props->init_toomre_length_softening_factor;
}

#endif /* SWIFT_MOSAICS_STARS_H */
