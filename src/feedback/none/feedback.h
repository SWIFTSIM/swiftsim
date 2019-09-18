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
#ifndef SWIFT_FEEDBACK_NONE_H
#define SWIFT_FEEDBACK_NONE_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {}

/**
 * @brief Should we do feedback for this star?
 *
 * @param sp The star to consider.
 */
__attribute__((always_inline)) INLINE static int feedback_do_feedback(
    const struct spart* sp) {

  return 0;
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * Note: Since this 'none' feedback mode is used for testing the neighbour
 * loops only, we want to always do feedback irrespective of the particle
 * or of the system's state.
 *
 * @param sp The #spart.
 * @param time The current simulation time (Non-cosmological runs).
 * @param cosmo The cosmological model (cosmological runs).
 * @param with_cosmology Are we doing a cosmological run?
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const float time, const struct cosmology* cosmo,
    const int with_cosmology) {

  return 1;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {}

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
    struct spart* sp, const struct feedback_props* feedback_props) {}

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
 * @brief changes the particle properties of gas particles flagged as HII region
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param hydro_properties the hydro_props struct
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param feedback_props feedback_props data structure
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param time Time since Big Bang
 *
 */
INLINE static void heating_HII_part(const struct phys_const* phys_const,
                                    const struct unit_system* us,
                                    const struct hydro_props* hydro_properties,
                                    const struct cosmology* cosmo,
                                    const struct cooling_function_data* cooling,
                                    const struct feedback_props* feedback_props,
                                    struct part* restrict p,
                                    struct xpart* restrict xp, double time) {}

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
 * @param time_beg_of_step Time at the beginning of the time step
 */
__attribute__((always_inline)) INLINE static void feedback_evolve_spart(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const double star_age_beg_step, const double dt,
    const double time_beg_of_step) {}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
static INLINE void feedback_struct_dump(const struct feedback_props* feedback,
                                        FILE* stream) {}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
static INLINE void feedback_struct_restore(struct feedback_props* feedback,
                                           FILE* stream) {}

#endif /* SWIFT_FEEDBACK_NONE_H */
