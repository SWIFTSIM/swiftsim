/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_STAR_H
#define SWIFT_DEFAULT_STAR_H

#include <float.h>
#include "minmax.h"

/**
 * @brief Computes the gravity time-step of a given star particle.
 *
 * @param sp Pointer to the s-particle data.
 */
__attribute__((always_inline)) INLINE static float star_compute_timestep(
    const struct spart* const sp) {

  return FLT_MAX;
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void star_first_init_spart(
    struct spart* sp) {

  sp->time_bin = 0;
}

/**
 * @brief Prepares a s-particle for its interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void star_init_spart(
    struct spart* sp) {

#ifdef DEBUG_INTERACTIONS_SPH
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS; ++i) sp->ids_ngbs_density[i] = -1;
  sp->num_ngb_density = 0;
#endif

  sp->wcount = 0.f;
  sp->wcount_dh = 0.f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The particle.
 */
__attribute__((always_inline)) INLINE static void star_reset_predicted_values(
    struct spart* restrict sp) {}

/**
 * @brief Finishes the calculation of (non-gravity) forces acting on stars
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void star_end_force(
    struct spart* sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void star_kick_extra(
    struct spart* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on stars
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void star_end_density(
    struct spart* sp) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  sp->wcount += kernel_root;
  sp->wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  sp->wcount *= h_inv_dim;
  sp->wcount_dh *= h_inv_dim_plus_one;
}


/**
 * @brief Sets all particle fields to sensible values when the #spart has 0 ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_spart_has_no_neighbours(
    struct spart *restrict sp, const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  sp->wcount = kernel_root * h_inv_dim;
  sp->wcount_dh = 0.f;
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_prepare_force(
    struct part *restrict sp, const struct cosmology *cosmo) {}


/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all star acceleration and time derivative fields in preparation
 * for the sums taking place in the variaous force tasks
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void star_reset_acceleration(
    struct spart *restrict p) {

#ifdef DEBUG_INTERACTIONS_SPH
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS; ++i) sp->ids_ngbs_force[i] = -1;
  sp->num_ngb_force = 0;
#endif

}

#endif /* SWIFT_DEFAULT_STAR_H */
