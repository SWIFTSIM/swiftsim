/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_GEAR_H
#define SWIFT_FEEDBACK_GEAR_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param cosmo The #cosmology.
 */
__attribute__((always_inline)) INLINE static void feedback_update_part(
    struct part *restrict p, struct xpart *restrict xp, const struct cosmology *restrict cosmo) {

  /* Update mass */
  const float mass = hydro_get_mass(p);

  hydro_set_mass(p, mass + xp->feedback_data.delta_mass);

  xp->feedback_data.delta_mass = 0;

  /* Update internal energy */
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float u_new = u + xp->feedback_data.delta_u;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, u_new);

  xp->feedback_data.delta_u = 0.;

  /* // TODO activate particle */


}

/**
 * @brief Should we do feedback for this star?
 *
 * @param sp The star to consider.
 */
__attribute__((always_inline)) INLINE static int feedback_do_feedback(
    const struct spart* sp) {

  return (sp->birth_time != -1.);
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param time The current simulation time (Non-cosmological runs).
 * @param cosmo The cosmological model (cosmological runs).
 * @param with_cosmology Are we doing a cosmological run?
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const float time, const struct cosmology* cosmo,
    const int with_cosmology) {

  if (sp->birth_time == -1.) return 0;

  if (with_cosmology) {
    return ((float)cosmo->a) > sp->birth_scale_factor;
  } else {
    return time > sp->birth_time;
  }
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.enrichment_weight = 0.f;

}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.to_distribute.mass = 0.f;

  /* Zero the energy to inject */
  sp->feedback_data.to_distribute.energy = 0.f;
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {

  sp->feedback_data.explosion_time = -1;
  feedback_init_spart(sp);
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 */
__attribute__((always_inline)) INLINE static void feedback_evolve_spart(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const double star_age_beg_step, const double dt) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv);       /* 1/h^d */
  
  sp->feedback_data.enrichment_weight *= hi_inv_dim;

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  /* Does the star need to explode? */
  if (sp->feedback_data.explosion_time > 0)
    return;

  if (sp->mass < feedback_props->mass_ejected)
    error("Stars cannot have negative mass. (%g < %g)",
	  sp->mass, feedback_props->mass_ejected);
  sp->mass -= feedback_props->mass_ejected;

}

#endif /* SWIFT_FEEDBACK_GEAR_H */
