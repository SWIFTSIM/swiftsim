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
 * @param xpj Hydro particle extra data.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_inject(
    const float r2, float *dx, const float hi, const float hj,
    struct spart *restrict si, struct part *restrict pj) {

  struct rt_spart_data *restrict sd = &(si->rt_data);
  struct rt_part_data *restrict pd = &(pj->rt_data);

  sd->calls_tot += 1;
  sd->calls_per_step += 1;
  sd->iact_hydro_inject += 1;

  pd->calls_tot += 1;
  pd->calls_per_step += 1;
  pd->iact_stars_inject += 1;

  if (r2 > 0.f) {
    sd->calls_self_inject += 1;
    pd->calls_self_inject += 1;
  } else {
    sd->calls_pair_inject += 1;
    pd->calls_pair_inject += 1;
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
 * @param mode 0 if non-symmetric interaction, 1 if symmetric
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_flux_common(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, int mode, long long ciID, long long cjID) {

  if (mode == 1) {
    pi->rt_data.calls_iact_transport += 1;
    pi->rt_data.calls_iact_transport_sym += 1;

    pj->rt_data.calls_iact_transport += 1;
    pj->rt_data.calls_iact_transport_sym += 1;

    int f; 
    f = pi->rt_data.neigh_iact_transp_free;
    if (f == 400) error("Reached 400 neighbours for transport particle debugging. Raise limit");
    pi->rt_data.neigh_iact_transp[f] = pj->id;
    pi->rt_data.neigh_cell_iact_transp[f] = cjID;
    pi->rt_data.neigh_iact_transp_free++;
    pi->rt_data.this_cell = ciID;

    f = pj->rt_data.neigh_iact_transp_free;
    if (f == 400) error("Reached 400 neighbours for transport particle debugging. Raise limit");
    pj->rt_data.neigh_iact_transp[f] = pi->id;
    pj->rt_data.neigh_cell_iact_transp[f] = ciID;
    pj->rt_data.neigh_iact_transp_free++;
    pj->rt_data.this_cell = cjID;


  } else {
    pi->rt_data.calls_iact_transport += 1;
    pi->rt_data.calls_iact_transport_nonsym += 1;

    int f = pi->rt_data.neigh_iact_transp_free;
    if (f == 400) error("Reached 400 neighbours for transport particle debugging. Raise limit");
    pi->rt_data.neigh_iact_transp[f] = pj->id;
    pi->rt_data.neigh_cell_iact_transp[f] = cjID;
    pi->rt_data.neigh_iact_transp_free++;
    pi->rt_data.this_cell = ciID;
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
    struct part *restrict pj, float a, float H, long long ciID, long long cjID) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 1, ciID, cjID);
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
                                struct part *restrict pj, float a, float H, long long ciID, long long cjID) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 0, ciID, cjID);
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
    struct part *restrict pj, float a, float H, long long ciID, long long cjID) {

  rt_gradients_collect(r2, dx, hi, hj, pi, pj, ciID, cjID);
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
                               struct part *restrict pj, float a, float H, long long ciID, long long cjID) {

  rt_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj, ciID, cjID);
}

#endif /* SWIFT_RT_IACT_DEBUG_H */
