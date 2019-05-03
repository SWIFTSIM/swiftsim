/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_GEAR_FEEDBACK_IACT_H
#define SWIFT_GEAR_FEEDBACK_IACT_H

/* Local includes */
#include "random.h"
#include "hydro.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float *dx,
                                    const float hi, const float hj,
                                    struct spart *restrict si,
                                    const struct part *restrict pj,
                                    const struct xpart *restrict xpj,
                                    const struct cosmology *restrict cosmo,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  /* The normalization by 1 / h^d is done in feedback.h */
  si->feedback_data.enrichment_weight += mj * wi;

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
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  struct spart *restrict si,
                                  struct part *restrict pj,
                                  struct xpart *restrict xpj,
				  const struct feedback_props* fp,
                                  const struct cosmology *restrict cosmo,
                                  const integertime_t ti_current) {

  const float mj = hydro_get_mass(pj);
  const float r = sqrtf(r2);

  /* Get the kernel for hi. */
  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv);       /* 1/h^d */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;

  /* Does the star need to explode? */
  if (si->feedback_data.explosion_time > 0 && si->feedback_data.explosion_time != ti_current)
    return;

  si->feedback_data.explosion_time = ti_current;
  message("exploding");

  /* Mass received */
  const double m_ej = fp->mass_ejected;
  // TODO compute inverse before feedback loop
  const double weight = mj * wi / si->feedback_data.enrichment_weight;
  const double dm = m_ej * weight;
  const double new_mass = mj + dm;

  /* Energy received */
  const double e_sn = fp->energy_per_supernovae;
  // TODO compute inverse before feedback loop
  const int n_sn = 1;
  const double u_sn = n_sn * e_sn / m_ej;
  const double du = dm * u_sn * weight / new_mass;

  xpj->feedback_data.delta_mass = dm;
  xpj->feedback_data.delta_u = du;

  /* Compute the norm of the speed. */
  float vi2 = 0.;
  float vj2 = 0.;
  for(int i = 0; i < 3; i++) {
    vi2 += si->v[i] * si->v[i];
    vj2 += pj->v[i] * pj->v[i];
  }

  /* Update velocity in order to conserve momentum */
  // TODO Update value in feedback_update_part
  float fac = (mj + weight * m_ej * vi2 / vj2) / new_mass;
  fac = sqrt(fac);
  for(int i = 0; i < 3; i++) {
    float dv = (fac - 1.) * pj->v[i];
    xpj->v_full[i] += dv;
    pj->v[i] += dv;
  }

}

#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
