/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_COLIBRE_STAR_FORMATION_IACT_H
#define SWIFT_COLIBRE_STAR_FORMATION_IACT_H

/**
 * @file COLIBRE/star_formation_iact.h
 * @brief Density computation
 */

/**
 * @brief do star_formation computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_star_formation(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  /* Calculation of the velocity dispersion */

  /* Calculating a handfull of often used variables*/
  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  /* Calculate the fraction of the kernel */
  const float ui = r * hi_inv;
  const float uj = r * hj_inv;

  /* Calculate the velocity difference vector */
  const float vi_min_vj[3] = {pi->v[0] - pj->v[0], pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};

  /* Calculate the constant for the position */
  const float a2H = a * a * H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + vi_min_vj[0],
                                  a2H * dx[1] + vi_min_vj[1],
                                  a2H * dx[2] + vi_min_vj[2]};

  /* Calculate the velocity dispersion */
  const float norm_v2 = v_plus_H_flow[0] * v_plus_H_flow[0] +
                        v_plus_H_flow[1] * v_plus_H_flow[1] +
                        v_plus_H_flow[2] * v_plus_H_flow[2];

  /* Calculate the kernel */
  float wi;
  float wj;

  kernel_eval(ui, &wi);
  kernel_eval(uj, &wj);

  /* Finalize the calculation of the velocity dispersion squared */
  pi->sf_data.sigma_v2 += norm_v2 * wi * pj->mass;
  pj->sf_data.sigma_v2 += norm_v2 * wj * pi->mass;

  pi->sf_data.gas_mass_unweighted += pj->mass;
  pj->sf_data.gas_mass_unweighted += pi->mass;
}

/**
 * @brief do star_formation computation after the runner_iact_density (non
 * symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_star_formation(float r2, const float *dx, float hi, float hj,
                                  struct part *restrict pi,
                                  const struct part *restrict pj, float a,
                                  float H) {

  /* Calculation of the velocity dispersion */

  /* Calculating a handfull of often used variables */
  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;

  /* Calculate the kernel fraction*/
  const float ui = r * hi_inv;

  /* Calculate the velocity difference */
  const float vi_min_vj[3] = {pi->v[0] - pj->v[0], pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};

  /* Calculate the constant for the position */
  const float a2H = a * a * H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + vi_min_vj[0],
                                  a2H * dx[1] + vi_min_vj[1],
                                  a2H * dx[2] + vi_min_vj[2]};

  /* Calculate the velocity dispersion */
  const float norm_v2 = v_plus_H_flow[0] * v_plus_H_flow[0] +
                        v_plus_H_flow[1] * v_plus_H_flow[1] +
                        v_plus_H_flow[2] * v_plus_H_flow[2];

  float wi;
  kernel_eval(ui, &wi);

  pi->sf_data.sigma_v2 += norm_v2 * wi * pj->mass;

  pi->sf_data.gas_mass_unweighted += pj->mass;
}

#endif /* SWIFT_COLIBRE_STAR_FORMATION_IACT_H */
