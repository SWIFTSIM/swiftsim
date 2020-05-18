/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenunuiv.nl)
 *               2020 Camila Correa (c.a.correa@uva.nl)
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

  struct chemistry_part_data *di = &pi->chemistry_data;
  struct chemistry_part_data *dj = &pj->chemistry_data;

  float dwi_dx, dwj_dx;

  /* Get the masses */
  const float mj = hydro_get_mass(pj);
  const float mi = hydro_get_mass(pi);

  /* Get r */
  const float r_inv = 1.f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_eval(ui, &dwi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_eval(uj, &dwj_dx);

  const float dwj_r = dwj_dx * r_inv;
  const float mi_dwj_r = mi * dwj_r;

  const float dwi_r = dwi_dx * r_inv;
  const float mj_dwi_r = mj * dwi_r;

  /* Compute velocity shear tensor */
  for (int k = 0; k < 3; k++) {
    di->shear_tensor[k][0] += (pj->v[0] - pi->v[0]) * dx[k] * mj_dwi_r;
    di->shear_tensor[k][1] += (pj->v[1] - pi->v[1]) * dx[k] * mj_dwi_r;
    di->shear_tensor[k][2] += (pj->v[2] - pi->v[2]) * dx[k] * mj_dwi_r;

    dj->shear_tensor[k][0] -= (pi->v[0] - pj->v[0]) * dx[k] * mi_dwj_r;
    dj->shear_tensor[k][1] -= (pi->v[1] - pj->v[1]) * dx[k] * mi_dwj_r;
    dj->shear_tensor[k][2] -= (pi->v[2] - pj->v[2]) * dx[k] * mi_dwj_r;
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

  struct chemistry_part_data *di = &pi->chemistry_data;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  float dwi_dx;

  /* Get r */
  const float r_inv = 1.f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_eval(ui, &dwi_dx);

  const float dwi_r = dwi_dx * r_inv;
  const float mj_dwi_r = mj * dwi_r;

  /* Compute velocity shear tensor */
  for (int k = 0; k < 3; k++) {
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, float time_base,
    integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  if (chj->diffusion_coefficient > 0 || chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float mi = hydro_get_mass(pi);
    const float rhoj = hydro_get_comoving_density(pj);
    const float rhoi = hydro_get_comoving_density(pi);

    float wi, wj, dwi_dx, dwj_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part j */
    /* Get the kernel for hj */
    const float hj_inv = 1.0f / hj;
    const float hj_inv_dim = pow_dimension(hj_inv);        /* 1/h^d */
    const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */
    const float rho_j_inv = 1.0f / rhoj;

    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;
    const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
    const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
    const float rho_i_inv = 1.0f / rhoi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);

    float dw_r = 0.5f *
                 (dwi_dx * hi_inv_dim_plus_one + dwj_dx * hj_inv_dim_plus_one) *
                 r_inv;

    const float mj_dw_r = mj * dw_r;
    const float mi_dw_r = mi * dw_r;

    /* Compute K_ij coefficient (see Correa et al., in prep.) */
    /* K_ij in physical coordinates */
    float K;
    K = 4.0f * chj->diffusion_coefficient * chi->diffusion_coefficient;
    K /= (chj->diffusion_coefficient + chi->diffusion_coefficient);
    K *= rho_i_inv * rho_j_inv;
    K *= a;

    const float K_ij = K * mj_dw_r;
    const float K_ji = K * mi_dw_r;

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
      chi->dmetal_mass_fraction[i] +=
          K_ij * (chi->metal_mass_fraction[i] - chj->metal_mass_fraction[i]) *
          dt;
      chj->dmetal_mass_fraction[i] +=
          K_ji * (chj->metal_mass_fraction[i] - chi->metal_mass_fraction[i]) *
          dt;

      chi->diffusion_rate[i] +=
          K_ij * (chi->metal_mass_fraction[i] - chj->metal_mass_fraction[i]);
      chj->diffusion_rate[i] +=
          K_ji * (chj->metal_mass_fraction[i] - chi->metal_mass_fraction[i]);
    }

    /* Diffusing the total metal mass */
    chi->dmetal_mass_fraction_total +=
        K_ij *
        (chi->metal_mass_fraction_total - chj->metal_mass_fraction_total) * dt;

    chj->dmetal_mass_fraction_total +=
        K_ji *
        (chj->metal_mass_fraction_total - chi->metal_mass_fraction_total) * dt;

    /* And diffuse the individual feedback tracers */
    chi->dmetal_mass_fraction_from_SNIa +=
        K_ij *
        (chi->metal_mass_fraction_from_SNIa -
         chj->metal_mass_fraction_from_SNIa) *
        dt;

    chj->dmetal_mass_fraction_from_SNIa +=
        K_ji *
        (chj->metal_mass_fraction_from_SNIa -
         chi->metal_mass_fraction_from_SNIa) *
        dt;

    chi->dmetal_mass_fraction_from_AGB += K_ij *
                                          (chi->metal_mass_fraction_from_AGB -
                                           chj->metal_mass_fraction_from_AGB) *
                                          dt;

    chj->dmetal_mass_fraction_from_AGB += K_ji *
                                          (chj->metal_mass_fraction_from_AGB -
                                           chi->metal_mass_fraction_from_AGB) *
                                          dt;

    chi->dmetal_mass_fraction_from_SNII +=
        K_ij *
        (chi->metal_mass_fraction_from_SNII -
         chj->metal_mass_fraction_from_SNII) *
        dt;

    chj->dmetal_mass_fraction_from_SNII +=
        K_ji *
        (chj->metal_mass_fraction_from_SNII -
         chi->metal_mass_fraction_from_SNII) *
        dt;

    chi->diron_mass_fraction_from_SNIa += K_ij *
                                          (chi->iron_mass_fraction_from_SNIa -
                                           chj->iron_mass_fraction_from_SNIa) *
                                          dt;

    chj->diron_mass_fraction_from_SNIa += K_ji *
                                          (chj->iron_mass_fraction_from_SNIa -
                                           chi->iron_mass_fraction_from_SNIa) *
                                          dt;

    float chi_Eu_fraction_from_NSM = chi->mass_from_NSM / mi;
    float chj_Eu_fraction_from_NSM = chj->mass_from_NSM / mj;
      chi->dEu_mass_fraction_from_NSM += K_ij * (chi_Eu_fraction_from_NSM - chj_Eu_fraction_from_NSM) * dt;

      chj->dEu_mass_fraction_from_NSM += K_ji * (chj_Eu_fraction_from_NSM - chi_Eu_fraction_from_NSM) * dt;

    float chi_Eu_fraction_from_CEJSN = chi->mass_from_CEJSN / mi;
    float chj_Eu_fraction_from_CEJSN = chj->mass_from_CEJSN / mj;
      chi->dEu_mass_fraction_from_CEJSN += K_ij * (chi_Eu_fraction_from_CEJSN - chj_Eu_fraction_from_CEJSN) * dt;
      
      chj->dEu_mass_fraction_from_CEJSN += K_ji * (chj_Eu_fraction_from_CEJSN - chi_Eu_fraction_from_CEJSN) * dt;

    float chi_Eu_fraction_from_collapsar = chi->mass_from_collapsar / mi;
    float chj_Eu_fraction_from_collapsar = chj->mass_from_collapsar / mj;
      chi->dEu_mass_fraction_from_collapsar += K_ij * (chi_Eu_fraction_from_collapsar - chj_Eu_fraction_from_collapsar) * dt;
      
      chj->dEu_mass_fraction_from_collapsar += K_ji * (chj_Eu_fraction_from_collapsar - chi_Eu_fraction_from_collapsar) * dt;
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, float time_base,
    integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  if (chj->diffusion_coefficient > 0 || chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mi = hydro_get_mass(pi);
    const float mj = hydro_get_mass(pj);
    const float rhoj = hydro_get_comoving_density(pj);
    const float rhoi = hydro_get_comoving_density(pi);

    float wi, wj, dwi_dx, dwj_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part j */
    /* Get the kernel for hj */
    const float hj_inv = 1.0f / hj;
    const float hj_inv_dim = pow_dimension(hj_inv);        /* 1/h^d */
    const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */
    const float rho_j_inv = 1.0f / rhoj;

    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);

    /* part i */
    /* Get the kernel for hi */
    float hi_inv = 1.0f / hi;
    const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
    const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
    const float rho_i_inv = 1.0f / rhoi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);
    float dw_r = 0.5f *
                 (dwi_dx * hi_inv_dim_plus_one + dwj_dx * hj_inv_dim_plus_one) *
                 r_inv;
    float mj_dw_r = mj * dw_r;

    /* Compute K_ij coefficient (see Correa et al., in prep.) */
    float K;
    K = 4.0f * chj->diffusion_coefficient * chi->diffusion_coefficient;
    K /= (chj->diffusion_coefficient + chi->diffusion_coefficient);
    K *= rho_i_inv * rho_j_inv;
    K *= a;

    const float K_ij = K * mj_dw_r;

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
      chi->dmetal_mass_fraction[i] +=
          K_ij * (chi->metal_mass_fraction[i] - chj->metal_mass_fraction[i]) *
          dt;
      chi->diffusion_rate[i] +=
          K_ij * (chi->metal_mass_fraction[i] - chj->metal_mass_fraction[i]);
    }

    /* Diffusing the total metal mass */
    chi->dmetal_mass_fraction_total +=
        K_ij *
        (chi->metal_mass_fraction_total - chj->metal_mass_fraction_total) * dt;

    /* And diffuse the individual feedback tracers */
    chi->dmetal_mass_fraction_from_SNIa +=
        K_ij *
        (chi->metal_mass_fraction_from_SNIa -
         chj->metal_mass_fraction_from_SNIa) *
        dt;

    chi->dmetal_mass_fraction_from_AGB += K_ij *
                                          (chi->metal_mass_fraction_from_AGB -
                                           chj->metal_mass_fraction_from_AGB) *
                                          dt;

    chi->dmetal_mass_fraction_from_SNII +=
        K_ij *
        (chi->metal_mass_fraction_from_SNII -
         chj->metal_mass_fraction_from_SNII) *
        dt;

    chi->diron_mass_fraction_from_SNIa += K_ij *
                                          (chi->iron_mass_fraction_from_SNIa -
                                           chj->iron_mass_fraction_from_SNIa) *
                                          dt;
      
    float chi_Eu_fraction_from_NSM = chi->mass_from_NSM / mi;
    float chj_Eu_fraction_from_NSM = chj->mass_from_NSM / mj;
    chi->dEu_mass_fraction_from_NSM += K_ij * (chi_Eu_fraction_from_NSM - chj_Eu_fraction_from_NSM) * dt;
      
    float chi_Eu_fraction_from_CEJSN = chi->mass_from_CEJSN / mi;
    float chj_Eu_fraction_from_CEJSN = chj->mass_from_CEJSN / mj;
    chi->dEu_mass_fraction_from_CEJSN += K_ij * (chi_Eu_fraction_from_CEJSN - chj_Eu_fraction_from_CEJSN) * dt;

    float chi_Eu_fraction_from_collapsar = chi->mass_from_collapsar / mi;
    float chj_Eu_fraction_from_collapsar = chj->mass_from_collapsar / mj;
    chi->dEu_mass_fraction_from_collapsar += K_ij * (chi_Eu_fraction_from_collapsar - chj_Eu_fraction_from_collapsar) * dt;
  }
}

#endif /* SWIFT_COLIBRE_CHEMISTRY_IACT_H */
