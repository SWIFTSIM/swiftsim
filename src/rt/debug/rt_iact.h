/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IACT_DEBUG_H
#define SWIFT_RT_IACT_DEBUG_H

#include "rt_gradients.h"

/**
 * @file src/rt/debug/rt_iact.h
 * @brief Main header file for the debug radiative transfer scheme particle
 * interactions.
 */

/**
 * @brief Injection step interaction between star and hydro particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si Star particle.
 * @param pj Hydro particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_inject(
    const float r2, float *dx, const float hi, const float hj,
    struct spart *restrict si, struct part *restrict pj, float a, float H) {

  struct rt_spart_data *restrict sd = &(si->rt_data);
  struct rt_part_data *restrict pd = &(pj->rt_data);

  sd->iact_hydro_inject += 1;
  sd->calls_tot += 1;
  sd->calls_per_step += 1;

  pd->iact_stars_inject += 1;
  pd->calls_tot += 1;
  pd->calls_per_step += 1;
}

/**
 * @brief Check the time step sizes of the star and hydro particle
 * are compatible. Debugging-only function.
 *
 * @param si Star particle.
 * @param pj Hydro particle.
 * @param e pointer to the engine.
 */
__attribute__((always_inline)) INLINE static void
rt_injection_timestep_debugging_check(struct spart *restrict sp,
                                      struct part *restrict p,
                                      const struct engine *e) {

  const integertime_t ti_current = e->ti_current;
  const integertime_t pti_end = get_integer_time_end(ti_current, p->time_bin);
  const integertime_t sti_end = get_integer_time_end(ti_current, sp->time_bin);

  if (sti_end < ti_current)
    error(
        "s-particle in an impossible time-zone! sp->ti_end=%lld "
        "e->ti_current=%lld",
        sti_end, ti_current);
  if (pti_end < ti_current)
    error(
        "particle in an impossible time-zone! p->ti_end=%lld "
        "e->ti_current=%lld",
        pti_end, ti_current);
  if (pti_end > sti_end)
    message(
        "WARNING: Got star that whose time step ends before the interacting "
        "hydro particle's time step ends. This needs to be dealt with.");
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_rt_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param mode 0 if non-symmetric interaction, 1 if symmetric
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_flux_common(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, int mode) {

  if (pi->rt_data.injection_done != 1)
    error(
        "Trying to do iact transport when "
        "finalise injection count is %d",
        pi->rt_data.injection_done);

  if (pi->rt_data.calls_iact_gradient == 0)
    error(
        "Called iact transport on particle "
        "with iact gradient count 0");

  if (pi->rt_data.gradients_done != 1)
    error(
        "Trying to do iact transport when "
        "rt_finalise_gradient count is %d",
        pi->rt_data.gradients_done);

  if (mode == 1) {

    if (pj->rt_data.injection_done != 1)
      error(
          "Trying to do iact transport when "
          "finalise injection count is %d",
          pj->rt_data.injection_done);

    if (pj->rt_data.calls_iact_gradient == 0)
      error(
          "Called iact transport on particle "
          "with iact gradient count 0");

    if (pj->rt_data.gradients_done != 1)
      error(
          "Trying to do iact transport when "
          "rt_finalise_gradient count is %d",
          pj->rt_data.gradients_done);

    pi->rt_data.calls_tot += 1;
    pi->rt_data.calls_per_step += 1;
    pi->rt_data.calls_iact_transport += 1;
    pj->rt_data.calls_tot += 1;
    pj->rt_data.calls_per_step += 1;
    pj->rt_data.calls_iact_transport += 1;
  } else {

    if (pj->rt_data.injection_done != 1)
      message(
          "Trying to do iact transport when finalise injection "
          "count is %d in nonsym flux. You should look into this",
          pj->rt_data.injection_done);

    if (pj->rt_data.calls_iact_gradient == 0)
      message(
          "Called iact transport on particle with iact gradient "
          "count 0 in nonsym flux. You should look into this");

    if (pj->rt_data.gradients_done != 1)
      message(
          "Trying to do iact transport when rt_finalise_gradient "
          "count is %d in nonsym flux. You should look into this",
          pj->rt_data.gradients_done);

    pi->rt_data.calls_tot += 1;
    pi->rt_data.calls_per_step += 1;
    pi->rt_data.calls_iact_transport += 1;
  }
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_rt_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_transport(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 1);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * This method calls runner_iact_rt_fluxes_common with mode 0.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_transport(float r2, const float *dx, float hi, float hj,
                                struct part *restrict pi,
                                struct part *restrict pj, float a, float H) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 0);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around rt_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  rt_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around rt_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_gradient(float r2, const float *dx, float hi, float hj,
                               struct part *restrict pi,
                               struct part *restrict pj, float a, float H,
                               unsigned long long cellID) {

  rt_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj, cellID);
}

#endif /* SWIFT_RT_IACT_DEBUG_H */
