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
#ifndef SWIFT_EAGLE_STARS_IACT_H
#define SWIFT_EAGLE_STARS_IACT_H

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
runner_iact_nonsym_stars_density(float r2, const float *dx, float hi, float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, float a,
                                 float H) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

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

  /* Compute contribution to the density */
  si->rho_gas += mj * wi;

  /* Mass of the gas neighbours */
  si->density.neighbour_mass += mj;
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
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
runner_iact_nonsym_stars_feedback(float r2, const float *dx, float hi, float hj,
                                  const struct spart *restrict spi,
                                  struct part *restrict pj, float a, float H) {

  /* Skip stars that don't do feedback */
  if (spi->feedback.probability <= 0.) return;

  /* Are we lucky? */
  const float prob_j = random_unit_interval(pj->id, spi->feedback.ti_current,
                                            random_number_stellar_feedback);

  /* Should we inject energy in this particle? */
  if (prob_j < spi->feedback.probability) {

    /* Inject energy into the gas particle */
    pj->feedback_data.delta_u += spi->feedback.delta_u;

    message("star %lld injecting energy in gas particle %lld. prob=%f r=%f",
            spi->id, pj->id, spi->feedback.probability, prob_j);
  }
}

#endif /* SWIFT_EAGLE_STARS_IACT_H */
