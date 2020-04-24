/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2020 Joel Pfeffer (j.l.pfeffer@ljmu.ac.uk)
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
#ifndef SWIFT_MOSAICS_STARS_IACT_H
#define SWIFT_MOSAICS_STARS_IACT_H

#include "random.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_density(const float r2, const float *dx,
                                 const float hi, const float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, const float a,
                                 const float H) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_density[si->num_ngb_density] = pj->id;
  ++si->num_ngb_density;
#endif
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_feedback(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  const struct spart *restrict si,
                                  struct part *restrict pj, const float a,
                                  const float H) {

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_feedback < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_feedback[si->num_ngb_feedback] = pj->id;
  ++si->num_ngb_feedback;
#endif
}

/**
 * @brief Stellar velocity dispersion/mass fraction computation (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param pi First particle.
 * @param si Second sparticle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_star_veldisp(const float r2, const float *dx,
                                const float hi, struct spart *restrict si, 
                                const struct spart *restrict sj, const float a,
                                const float H) {

  float wi;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Number of particles which contribute to the calculation */
  si->scount += 1;

  /* Compute contribution to stellar density */
  si->stars_rho += wi * sj->mass;
  si->stars_mass_unweighted += sj->mass;

  /* Calculation of the velocity dispersion */

  /* Calculate the velocity difference */
  const float vi_min_vj[3] = {si->v[0]-sj->v[0], si->v[1]-sj->v[1], 
      si->v[2]-sj->v[2]};

  /* Calculate the constant for the position */
  const float a2H = a*a*H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + vi_min_vj[0], 
      a2H * dx[1] + vi_min_vj[1], a2H * dx[2] + vi_min_vj[2]};

  /* Calculate the velocity dispersion */
  const float norm_v2 = v_plus_H_flow[0] * v_plus_H_flow[0] + 
      v_plus_H_flow[1] * v_plus_H_flow[1] + 
      v_plus_H_flow[2] * v_plus_H_flow[2];

  si->stars_sigma_v2 += norm_v2 * wi * sj->mass;
}

#endif /* SWIFT_MOSAICS_STARS_IACT_H */
