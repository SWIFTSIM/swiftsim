/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Camila Correa (correa@strw.leidenuniv.nl)
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
#ifndef SWIFT_DIFFUSION_IACT_H
#define SWIFT_DIFFUSION_IACT_H

/**
 * @file DiffusionOfScalar/diffusion_iact.h
 * @brief Functions that calculate the tensor of the gas velocity shear,
 * diffusion coefficient, diffusion rate and diffused scalar
 * following the Smagorinsky model described in Correa et al (2019).
 */

/**
 * @brief do shear tensor computation after the runner_iact_density in tools.c
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static void runner_iact_shear_tensor(
        float r2, float hi, float hj, struct part *restrict pi,
        struct part *restrict pj) {
}

/**
 * @brief do shear tensor computation after the runner_iact_nonsym_density in tools.c
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_shear_tensor(
    float r2, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {
}
#endif /* SWIFT_DIFFUSION_IACT_H */
