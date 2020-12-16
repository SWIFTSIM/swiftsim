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
#ifndef SWIFT_RT_GRADIENTS_DEBUG_H
#define SWIFT_RT_GRADIENTS_DEBUG_H

#include "atomic.h"

/**
 * @file src/rt/debug/rt_gradients.h
 * @brief Main header file for the debug radiative transfer scheme gradients
 */

/**
 * @brief symmetric gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, long long ciID, long long cjID) {

  if (pi->rt_data.ghost_finished != 1) 
    printf("--- RT grad sym: Particle %6lld has ghost_finished = %i\n", pi->id, pi->rt_data.ghost_finished);
  if (pj->rt_data.ghost_finished != 1) 
    printf("--- RT grad sym: Particle %6lld has ghost_finished = %i\n", pj->id, pj->rt_data.ghost_finished);

  atomic_inc(&pi->rt_data.calls_iact_gradient);
  atomic_inc(&pi->rt_data.calls_iact_gradient_sym);

  atomic_inc(&pj->rt_data.calls_iact_gradient);
  atomic_inc(&pj->rt_data.calls_iact_gradient_sym);

  int f;
  f = pi->rt_data.neigh_iact_grad_free;
  if (f == 400) error("Reached 400 neighbours for grad particle debugging. Raise limit");
  pi->rt_data.neigh_iact_grad[f] = pj->id;
  pi->rt_data.neigh_cell_iact_grad[f] = cjID;
  pi->rt_data.neigh_iact_grad_free++;
  if (llabs(pi->rt_data.this_cell_grad) < llabs(ciID))
    pi->rt_data.this_cell_grad = ciID;
  pi->rt_data.h_grad = hi;

  f = pj->rt_data.neigh_iact_grad_free;
  if (f == 400) error("Reached 400 neighbours for grad particle debugging. Raise limit");
  pj->rt_data.neigh_iact_grad[f] = pi->id;
  pj->rt_data.neigh_cell_iact_grad[f] = ciID;
  pj->rt_data.neigh_iact_grad_free++;
  if (llabs(pj->rt_data.this_cell_grad) < llabs(cjID))
    pj->rt_data.this_cell_grad = cjID;
  pj->rt_data.h_grad = hj;
}

/**
 * @brief Non-symmetric gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_nonsym_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, long long ciID, long long cjID) {

  if (pi->rt_data.ghost_finished != 1) 
    printf("--- RT grad nonsym: Particle %6lld has ghost_finished = %i\n", pi->id, pi->rt_data.ghost_finished);
  if (pj->rt_data.ghost_finished != 1) 
    printf("--- RT grad nonsym: Particle %6lld has ghost_finished = %i\n", pj->id, pj->rt_data.ghost_finished);

  atomic_inc(&pi->rt_data.calls_iact_gradient);
  atomic_inc(&pi->rt_data.calls_iact_gradient_nonsym);

  int f = pi->rt_data.neigh_iact_grad_free;
  if (f == 400) error("Reached 400 neighbours for grad particle debugging. Raise limit");
  pi->rt_data.neigh_iact_grad[f] = pj->id;
  pi->rt_data.neigh_cell_iact_grad[f] = cjID;
  pi->rt_data.neigh_iact_grad_free++;
  if (llabs(pi->rt_data.this_cell_grad) < llabs(ciID))
    pi->rt_data.this_cell_grad = ciID;
  pi->rt_data.h_grad = hi;

}

#endif /* SWIFT_RT_GRADIENT_DEBUG_H */