/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_HYDRO_IACT_H
#define SWIFT_DEFAULT_HYDRO_IACT_H

/**
 * @file Default/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added diffusive physics (Cullen & Denhen 2011 AV,
 *        Price 2017 (PHANTOM) diffusion) (interaction routines)
 */

#include "adiabatic_index.h"
#include "minmax.h"

#include "./hydro_parameters.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    struct part* pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;
  float dv[3], curlvr[3];

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  /* Now we need to compute the div terms */
  const float r_inv = 1.f / r;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;
  pj->viscosity.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  /* Negative because of the change in sign of dx & dv. */
  pj->density.rot_v[0] += facj * curlvr[0];
  pj->density.rot_v[1] += facj * curlvr[1];
  pj->density.rot_v[2] += facj * curlvr[2];
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    const struct part* pj, float a, float H) {

  float wi, wi_dx;
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  const float r_inv = 1.f / r;
  const float faci = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
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
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);
  pj->viscosity.v_sig = max(pj->viscosity.v_sig, new_v_sig);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    struct part* pj, float a, float H) {

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float v_sig = pi->force.soundspeed + pj->force.soundspeed -
                      const_viscosity_beta * mu_ij;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  /* Balsara includes the alphas */
  const float visc = -0.125f * v_sig * mu_ij * (balsara_i + balsara_j);

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * pi->force.f / rhoi + wj_dr * pj->force.f / rhoj) *
      r_inv;

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * pi->force.f;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * pj->force.f;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = P_over_rho2_j * dvdr * r_inv * wj_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */

  const float v_diff =
      sqrtf(2.0 * fabsf(pressurei - pressurej) / (rhoi + rhoj)) +
      dvdr_Hubble * r_inv;

  const float alpha_diff = 0.5f * (pi->diffusion.alpha + pj->diffusion.alpha);
  /* wi_dx + wj_dx / 2 is F_ij */
  const float diff_du_term =
      alpha_diff * v_diff * (pi->u - pj->u) * 0.5f *
      (wi_dr * pi->force.f / pi->rho + wj_dr * pi->force.f / pi->rho);

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term - diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * pi->force.f * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * pj->force.f * r_inv / rhoi * wj_dr;
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    const struct part* pj, float a, float H) {

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  // const float mi = pi->mass;
  const float mj = pj->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float v_sig = pi->force.soundspeed + pj->force.soundspeed -
                      const_viscosity_beta * mu_ij;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  /* Balsara includes the alphas */
  const float visc = -0.125f * v_sig * mu_ij * (balsara_i + balsara_j);

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * pi->force.f / rhoi + wj_dr * pj->force.f / rhoj) *
      r_inv;

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * pi->force.f;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * pj->force.f;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  const float v_diff =
      sqrtf(2.0 * fabsf(pressurei - pressurej) / (rhoi + rhoj)) +
      dvdr_Hubble * r_inv;

  const float alpha_diff = 0.5f * (pi->diffusion.alpha + pj->diffusion.alpha);
  /* wi_dx + wj_dx / 2 is F_ij */
  const float diff_du_term =
      alpha_diff * v_diff * (pi->u - pj->u) * 0.5f *
      (wi_dr * pi->force.f / pi->rho + wj_dr * pi->force.f / pi->rho);

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * pi->force.f * r_inv / rhoj * wi_dr;
}

/**
 * @brief Timestep limiter loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_limiter(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* Nothing to do here if both particles are active */
}

/**
 * @brief Timestep limiter loop (non-symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_limiter(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* Wake up the neighbour? */
  if (pi->viscosity.v_sig >
      const_limiter_max_v_sig_ratio * pj->viscosity.v_sig) {

    pj->wakeup = time_bin_awake;
  }
}

#endif /* SWIFT_DEFAULT_HYDRO_IACT_H */
