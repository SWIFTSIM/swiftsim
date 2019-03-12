/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COLIBRE_CHEMISTRY_IACT_H
#define SWIFT_COLIBRE_CHEMISTRY_IACT_H

/**
 * @file COLIBRE/chemistry_iact.h
 * @brief Smooth metal interaction functions following the COLIBRE model.
 */

/**
 * @brief do chemistry computation after the runner_iact_density (symmetric
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
 *
 * @brief after chemistry computation comes shear tensor computation
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, dwi_dx;
  float wj, dwj_dx;

  /* Get the masses. */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_eval(ui, &wi, &dwi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_eval(uj, &wj, &dwj_dx);

  /* Compute contribution to the smooth metallicity */
  for (int i = 0; i < chemistry_element_count; i++) {
    chi->smoothed_metal_mass_fraction[i] +=
        mj * chj->metal_mass_fraction[i] * wi;
    chj->smoothed_metal_mass_fraction[i] +=
        mi * chi->metal_mass_fraction[i] * wj;
  }

  // Smooth metal mass fraction of all metals
  chi->smoothed_metal_mass_fraction_total +=
      mj * chj->metal_mass_fraction_total * wi;
  chj->smoothed_metal_mass_fraction_total +=
      mi * chi->metal_mass_fraction_total * wj;

  // Smooth iron mass fraction from SNIa
  chi->smoothed_iron_mass_fraction_from_SNIa +=
      mj * chj->iron_mass_fraction_from_SNIa * wi;
  chj->smoothed_iron_mass_fraction_from_SNIa +=
      mi * chi->iron_mass_fraction_from_SNIa * wj;
    
  // Let's now calculate the shear tensor
  struct diffusion_part_data *di = &pi->diffusion_data;
  struct diffusion_part_data *dj = &pj->diffusion_data;
    
  float dwj_r = dwj_dx / r;
  float mi_dwj_r = mi * dwj_r;

  float dwi_r = dwi_dx / r;
  float mj_dwi_r = mj * dwi_r;
    
  /* Compute shear tensor */
  for (int k = 0; k < 3; k++){
        di->shear_tensor[k][0] += (pj->v[0] - pi->v[0]) * dx[k] * mj_dwi_r;
        di->shear_tensor[k][1] += (pj->v[1] - pi->v[1]) * dx[k] * mj_dwi_r;
        di->shear_tensor[k][2] += (pj->v[2] - pi->v[2]) * dx[k] * mj_dwi_r;
        
        dj->shear_tensor[k][0] += (pi->v[0] - pj->v[0]) * dx[k] * mi_dwj_r;
        dj->shear_tensor[k][1] += (pi->v[1] - pj->v[1]) * dx[k] * mi_dwj_r;
        dj->shear_tensor[k][2] += (pi->v[2] - pj->v[2]) * dx[k] * mi_dwj_r;
  }
}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 * @brief after chemistry computation comes shear tensor computation
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, dwi_dx;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_eval(ui, &wi, &dwi_dx);

  /* Compute contribution to the smooth metallicity */
  for (int i = 0; i < chemistry_element_count; i++) {
    chi->smoothed_metal_mass_fraction[i] +=
        mj * chj->metal_mass_fraction[i] * wi;
  }

  // Smooth metal mass fraction of all metals
  chi->smoothed_metal_mass_fraction_total +=
      mj * chj->metal_mass_fraction_total * wi;

  // Smooth iron mass fraction from SNIa
  chi->smoothed_iron_mass_fraction_from_SNIa +=
      mj * chj->iron_mass_fraction_from_SNIa * wi;
    
  /* Compute shear tensor */
   float dwi_r = dwi_dx / r;
   float mj_dwi_r = mj * dwi_r;
    
   for (int k = 0; k < 3; k++){
        di->shear_tensor[k][0] += (pj->v[0] - pi->v[0]) * dx[k] * mj_dwi_r;
        di->shear_tensor[k][1] += (pj->v[1] - pi->v[1]) * dx[k] * mj_dwi_r;
        di->shear_tensor[k][2] += (pj->v[2] - pi->v[2]) * dx[k] * mj_dwi_r;
   }
}


/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
                                                                        float r2, float hi, float hj, struct part *restrict pi,
                                                                        struct part *restrict pj, float a, float H, float time_base, integertime_t t_current,
                                                                        const struct cosmology *cosmo, const int with_cosmology) {
    
    struct diffusion_part_data *di = &pi->diffusion_data;
    struct diffusion_part_data *dj = &pj->diffusion_data;
    
    struct chemistry_part_data *chi = &pi->chemistry_data;
    struct chemistry_part_data *chj = &pj->chemistry_data;
    
    /* Get mass */
    float mj = pj->mass;
    float mi = pi->mass;
    float wi, wj, dwi_dx, dwj_dx;
    
    /* Get r */
    float r = sqrtf(r2);
    
    /* part j */
    /* Get the kernel for hj */
    float hj_inv = 1.0f / hj;
    const float hj_inv_dim = pow_dimension(hj_inv);       /* 1/h^d */
    const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */
    const float rho_j_inv = 1.0f / pj->rho;
    
    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);
    
    /* part i */
    /* Get the kernel for hi */
    float hi_inv = 1.0f / hi;
    const float hi_inv_dim = pow_dimension(hi_inv);       /* 1/h^d */
    const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
    const float rho_i_inv = 1.0f / pi->rho;
    
    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);
    
    float dw_r = 0.5f * (dwi_dx * hi_inv_dim_plus_one + dwj_dx * hj_inv_dim_plus_one) / r;
    float mj_dw_r = mj * dw_r;
    
    /* Compute K_ij coefficient (see Correa et al., in prep.) */
    float K_ij;
    K_ij = 4.0f * dj->diffusion_coefficient * di->diffusion_coefficient;
    K_ij /= (dj->diffusion_coefficient + di->diffusion_coefficient);
    K_ij *= rho_i_inv * rho_j_inv * mj_dw_r;
    float K_ji = K_ij * mi / mj;
    
    /* Manage time interval of particles i & j to be the smallest */
    double dt;
    if (with_cosmology) {
        const integertime_t ti_step = get_integer_timestep(pi->time_bin);
        const integertime_t ti_begin =
        get_integer_time_begin(t_current - 1, pi->time_bin);
        dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
    } else {
        dt = get_timestep(pi->time_bin, time_base);
    }
    /* Same for particle j */
    double dt_j;
    if (with_cosmology) {
        const integertime_t tj_step = get_integer_timestep(pj->time_bin);
        const integertime_t tj_begin =
        get_integer_time_begin(t_current - 1, pj->time_bin);
        dt_j = cosmology_get_delta_time(cosmo, tj_begin, tj_begin + tj_step);
    } else {
        dt_j = get_timestep(pj->time_bin, time_base);
    }
    if (dt_j < dt) dt = dt_j;
    
    /* Compute contribution to the metal abundance */
    for (int i = 0; i < chemistry_element_count; i++) {
        di->dmetal_mass_fraction[i] += K_ij * (chi->metal_mass_fraction[i] - chj->metal_mass_fraction[i]) * dt;
        dj->dmetal_mass_fraction[i] += K_ji * (chj->metal_mass_fraction[i] - chi->metal_mass_fraction[i]) * dt;
    }
    
}


/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
                                                                                      float r2, float hi, float hj, struct part *restrict pi,
                                                                                      struct part *restrict pj, float a, float H, float time_base, integertime_t t_current,
                                                                                      const struct cosmology *cosmo, const int with_cosmology) {
    struct diffusion_part_data *di = &pi->diffusion_data;
    struct diffusion_part_data *dj = &pj->diffusion_data;
    
    struct chemistry_part_data *chi = &pi->chemistry_data;
    struct chemistry_part_data *chj = &pj->chemistry_data;
    
    /* Get mass */
    float mj = pj->mass;
    float wi, wj, dwi_dx, dwj_dx;
    
    /* Get r */
    float r = sqrtf(r2);
    
    /* part j */
    /* Get the kernel for hj */
    float hj_inv = 1.0f / hj;
    const float hj_inv_dim = pow_dimension(hj_inv);       /* 1/h^d */
    const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */
    const float rho_j_inv = 1.0f / pj->rho;
    
    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);
    
    /* part i */
    /* Get the kernel for hi */
    float hi_inv = 1.0f / hi;
    const float hi_inv_dim = pow_dimension(hi_inv);       /* 1/h^d */
    const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
    const float rho_i_inv = 1.0f / pi->rho;
    
    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);
    
    float dw_r = 0.5f * (dwi_dx * hi_inv_dim_plus_one + dwj_dx * hj_inv_dim_plus_one) / r;
    float mj_dw_r = mj * dw_r;
    
    /* Compute K_ij coefficient (see Correa et al., in prep.) */
    float K_ij;
    K_ij = 4.0f * dj->diffusion_coefficient * di->diffusion_coefficient;
    K_ij /= (dj->diffusion_coefficient + di->diffusion_coefficient);
    K_ij *= rho_i_inv * rho_j_inv * mj_dw_r;
    
    /* Manage time interval of particles i & j to be the smallest */
    double dt;
    if (with_cosmology) {
        const integertime_t ti_step = get_integer_timestep(pi->time_bin);
        const integertime_t ti_begin =
        get_integer_time_begin(t_current - 1, pi->time_bin);
        dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
    } else {
        dt = get_timestep(pi->time_bin, time_base);
    }
    /* Same for particle j */
    double dt_j;
    if (with_cosmology) {
        const integertime_t tj_step = get_integer_timestep(pj->time_bin);
        const integertime_t tj_begin =
        get_integer_time_begin(t_current - 1, pj->time_bin);
        dt_j = cosmology_get_delta_time(cosmo, tj_begin, tj_begin + tj_step);
    } else {
        dt_j = get_timestep(pj->time_bin, time_base);
    }
    if (dt_j < dt) dt = dt_j;
    
    /* Compute contribution to the metal abundance */
    for (int i = 0; i < chemistry_element_count; i++) {
        di->dmetal_mass_fraction[i] += K_ij * (chi->metal_mass_fraction[i] - chj->metal_mass_fraction[i]) * dt;
    }
    
}

#endif /* SWIFT_COLIBRE_CHEMISTRY_IACT_H */
