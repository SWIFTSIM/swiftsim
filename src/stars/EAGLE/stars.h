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
#ifndef SWIFT_EAGLE_STARS_H
#define SWIFT_EAGLE_STARS_H

#include <float.h>

#include "dimension.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "inline.h"
#include "minmax.h"
#include "part.h"
#include "stars_part.h"
#include "stars_properties.h"

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

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
  sp->density.neighbour_mass = 0.f;
  sp->rho_gas = 0.f;
  sp->feedback.probability = 0.f;
  sp->feedback.delta_u = 0.f;
  sp->feedback.ti_current = -1;
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp) {

  sp->time_bin = 0;
  sp->birth_density = -1.f;
  sp->birth_time = -1.f;

  stars_init_spart(sp);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * Nothing to do in the EAGLE model. No additional field needs drifting.
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
 * Multiplies the forces and accelerations by the appropiate constants.
 * And reset properties.
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_end_feedback(
    struct spart* sp) {

  /* Make sure this star won't do feedback again */
  sp->feedback.probability = 0.f;
}

/**
 * @brief Kick the additional variables
 *
 * Nothing to do here in the EAGLE model. Stars do not carry additional
 * variables that are integrated forward in time.
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void stars_kick_extra(
    struct spart* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on stars.
 *
 * This loop is used to multiply in missing factors.
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
  sp->rho_gas *= h_inv_dim;
  sp->density.wcount *= h_inv_dim;
  sp->density.wcount_dh *= h_inv_dim_plus_one;

  /* Note: Nothing to do for the total neighbour gas mass */
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * This should only happen in pathological cases and is here to be used as
 * a damage mitigation exercise.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_spart_has_no_neighbours(
    struct spart* restrict sp, const struct cosmology* cosmo) {

  /* Re-set problematic values */
  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
  sp->density.neighbour_mass = 0.f;
  sp->rho_gas = 0.f;
}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 * @param stars_properties The #stars_props
 */
void stars_evolve_spart(struct spart* restrict sp,
                        const struct stars_props* stars_properties,
                        const struct cosmology* cosmo);

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
 * @param with_cosmology Are we running with cosmology on?
 * @param ti_current The current time on the time-line (for random numbers).
 */
void stars_prepare_feedback(
    struct spart* restrict sp, const struct stars_props* stars_properties,
    const struct hydro_props* hydro_properties, const struct unit_system* us,
    const struct phys_const* pyhs_consts, const struct cosmology* cosmo,
    const int with_cosmology, const integertime_t ti_current,
    const double time_base);

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * Nothing to do there in the EAGLE model as we do not accumulate quantities
 * in the feedback loop.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_feedback(
    struct spart* restrict p) {}

/**
 * @brief Reset the energy received by a particle
 */
__attribute__((always_inline)) INLINE static void feedback_part_reset(
    struct part* p) {

#ifdef SWIFT_DEBUG_CHECKS
  if (p->feedback_data.delta_u != 0.f)
    error("Feedback energy was not injected for particle %lld", p->id);
#endif

  // p->feedback_data.delta_u = 0.f;
}

__attribute__((always_inline)) INLINE static void feedback_apply(
    struct part* p, struct xpart* xp, const struct cosmology* cosmo,
    const integertime_t ti_current) {

  /* Anything to do here? */
  if (p->feedback_data.delta_u == 0.f) return;

  /* Apply the total change in energy from all the SN events */
  const float delta_u = p->feedback_data.delta_u;
  const float old_u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float new_u = old_u + delta_u;

  /* Update the internal energy and the drifted one. */
  hydro_set_physical_internal_energy(p, xp, cosmo, new_u);
  hydro_set_drifted_physical_internal_energy(p, cosmo, new_u);

  /* Be safe, make sure the energy injection is reset */
  p->feedback_data.delta_u = 0.f;
}

#endif /* SWIFT_EAGLE_STARS_H */
