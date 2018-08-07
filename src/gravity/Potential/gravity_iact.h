/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_POTENTIAL_GRAVITY_IACT_H
#define SWIFT_POTENTIAL_GRAVITY_IACT_H

/* Includes. */
#include "kernel_gravity.h"
#include "kernel_long_gravity.h"
#include "multipole.h"

/* Standard headers */
#include <float.h>

/**
 * @brief Computes the intensity of the force at a point generated by a
 * point-mass.
 *
 * The returned quantity needs to be multiplied by the distance vector to obtain
 * the force vector.
 *
 * @param r2 Square of the distance to the point-mass.
 * @param h2 Square of the softening length.
 * @param h_inv Inverse of the softening length.
 * @param h_inv3 Cube of the inverse of the softening length.
 * @param mass Mass of the point-mass.
 * @param f_ij (return) The force intensity.
 * @param pot_ij (return) The potential.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_full(
    const float r2, const float h2, const float h_inv, const float h_inv3,
    const float mass, float *f_ij, float *pot_ij) {

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2 + FLT_MIN);

  /* Should we soften ? */
  if (r2 >= h2) {

    /* Get Newtonian gravity */
    *f_ij = mass * r_inv * r_inv * r_inv;
    *pot_ij = -mass * r_inv;

  } else {

    const float r = r2 * r_inv;
    const float ui = r * h_inv;

    float W_f_ij, W_pot_ij;
    kernel_grav_force_eval(ui, &W_f_ij);
    kernel_grav_pot_eval(ui, &W_pot_ij);

    /* Get softened gravity */
    *f_ij = mass * h_inv3 * W_f_ij;
    *pot_ij = mass * h_inv * W_pot_ij;
  }
}

/**
 * @brief Computes the intensity of the force at a point generated by a
 * point-mass truncated for long-distance periodicity.
 *
 * The returned quantity needs to be multiplied by the distance vector to obtain
 * the force vector.
 *
 * @param r2 Square of the distance to the point-mass.
 * @param h2 Square of the softening length.
 * @param h_inv Inverse of the softening length.
 * @param h_inv3 Cube of the inverse of the softening length.
 * @param mass Mass of the point-mass.
 * @param r_s_inv Inverse of the mesh smoothing scale.
 * @param f_ij (return) The force intensity.
 * @param pot_ij (return) The potential.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_truncated(
    const float r2, const float h2, const float h_inv, const float h_inv3,
    const float mass, const float r_s_inv, float *f_ij, float *pot_ij) {

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2 + FLT_MIN);
  const float r = r2 * r_inv;

  /* Should we soften ? */
  if (r2 >= h2) {

    /* Get Newtonian gravity */
    *f_ij = mass * r_inv * r_inv * r_inv;
    *pot_ij = -mass * r_inv;

  } else {

    const float ui = r * h_inv;
    float W_f_ij, W_pot_ij;

    kernel_grav_force_eval(ui, &W_f_ij);
    kernel_grav_pot_eval(ui, &W_pot_ij);

    /* Get softened gravity */
    *f_ij = mass * h_inv3 * W_f_ij;
    *pot_ij = mass * h_inv * W_pot_ij;
  }

  /* Get long-range correction */
  const float u_lr = r * r_s_inv;
  float corr_f_lr, corr_pot_lr;
  kernel_long_grav_force_eval(u_lr, &corr_f_lr);
  kernel_long_grav_pot_eval(u_lr, &corr_pot_lr);
  *f_ij *= corr_f_lr;
  *pot_ij *= corr_pot_lr;
}

/**
 * @brief Computes the forces at a point generated by a multipole.
 *
 * This assumes M_100 == M_010 == M_001 == 0.
 * This uses the quadrupole and trace of the octupole terms only and defaults to
 * the monopole if the code is compiled with low-order gravity only.
 *
 * @param r_x x-component of the distance vector to the multipole.
 * @param r_y y-component of the distance vector to the multipole.
 * @param r_z z-component of the distance vector to the multipole.
 * @param r2 Square of the distance vector to the multipole.
 * @param h The softening length.
 * @param h_inv Inverse of the softening length.
 * @param m The multipole.
 * @param f_x (return) The x-component of the acceleration.
 * @param f_y (return) The y-component of the acceleration.
 * @param f_z (return) The z-component of the acceleration.
 * @param pot (return) The potential.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pm_full(
    const float r_x, const float r_y, const float r_z, const float r2,
    const float h, const float h_inv, const struct multipole *m, float *f_x,
    float *f_y, float *f_z, float *pot) {

/* In the case where the order is < 3, then there is only a monopole term left.
 * We can default to the normal P-P interaction with the mass of the multipole
 * and its CoM as the "particle" property */
#if SELF_GRAVITY_MULTIPOLE_ORDER < 3

  float f_ij, pot_ij;
  runner_iact_grav_pp_full(r2, h * h, h_inv, h_inv * h_inv * h_inv, m->M_000,
                           &f_ij, &pot_ij);
  *f_x = f_ij * r_x;
  *f_y = f_ij * r_y;
  *f_z = f_ij * r_z;
  *pot = pot_ij;

#else

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);

  /* Compute the derivatives of the potential */
  struct potential_derivatives_M2P d;
  compute_potential_derivatives_M2P(r_x, r_y, r_z, r2, r_inv, h, h_inv, 0, 0.f,
                                    &d);

  /* 0th order contributions */
  *f_x = m->M_000 * d.D_100;
  *f_y = m->M_000 * d.D_010;
  *f_z = m->M_000 * d.D_001;
  *pot = m->M_000 * d.D_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

/* 1st order contributions */

/* 1st order contributions are all 0 since the dipole is 0 */

/* *f_x = m->M_001 * d.D_101 + m->M_010 * d.D_110 + m->M_100 * d.D_200 ; */
/* *f_y = m->M_001 * d.D_011 + m->M_010 * d.D_020 + m->M_100 * d.D_110 ; */
/* *f_z = m->M_001 * d.D_002 + m->M_010 * d.D_011 + m->M_100 * d.D_101 ; */
/* *pot = m->M_001 * d.D_001 + m->M_010 * d.D_010 + m->M_100 * d.D_100 ; */

#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 2nd order contributions */
  *f_x += m->M_002 * d.D_102 + m->M_011 * d.D_111 + m->M_020 * d.D_120 +
          m->M_101 * d.D_201 + m->M_110 * d.D_210 + m->M_200 * d.D_300;
  *f_y += m->M_002 * d.D_012 + m->M_011 * d.D_021 + m->M_020 * d.D_030 +
          m->M_101 * d.D_111 + m->M_110 * d.D_120 + m->M_200 * d.D_210;
  *f_z += m->M_002 * d.D_003 + m->M_011 * d.D_012 + m->M_020 * d.D_021 +
          m->M_101 * d.D_102 + m->M_110 * d.D_111 + m->M_200 * d.D_201;
  *pot += m->M_002 * d.D_002 + m->M_011 * d.D_011 + m->M_020 * d.D_020 +
          m->M_101 * d.D_101 + m->M_110 * d.D_110 + m->M_200 * d.D_200;

#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 3rd order contributions */
  *f_x += m->M_003 * d.D_103 + m->M_012 * d.D_112 + m->M_021 * d.D_121 +
          m->M_030 * d.D_130 + m->M_102 * d.D_202 + m->M_111 * d.D_211 +
          m->M_120 * d.D_220 + m->M_201 * d.D_301 + m->M_210 * d.D_310 +
          m->M_300 * d.D_400;
  *f_y += m->M_003 * d.D_013 + m->M_012 * d.D_022 + m->M_021 * d.D_031 +
          m->M_030 * d.D_040 + m->M_102 * d.D_112 + m->M_111 * d.D_121 +
          m->M_120 * d.D_130 + m->M_201 * d.D_211 + m->M_210 * d.D_220 +
          m->M_300 * d.D_310;
  *f_z += m->M_003 * d.D_004 + m->M_012 * d.D_013 + m->M_021 * d.D_022 +
          m->M_030 * d.D_031 + m->M_102 * d.D_103 + m->M_111 * d.D_112 +
          m->M_120 * d.D_121 + m->M_201 * d.D_202 + m->M_210 * d.D_211 +
          m->M_300 * d.D_301;
  *pot += m->M_003 * d.D_003 + m->M_012 * d.D_012 + m->M_021 * d.D_021 +
          m->M_030 * d.D_030 + m->M_102 * d.D_102 + m->M_111 * d.D_111 +
          m->M_120 * d.D_120 + m->M_201 * d.D_201 + m->M_210 * d.D_210 +
          m->M_300 * d.D_300;

#endif

  /* Take care of the the sign convention */
  *f_x *= -1.f;
  *f_y *= -1.f;
  *f_z *= -1.f;
  *pot *= -1.f;
#endif
}

/**
 * @brief Computes the forces at a point generated by a multipole, truncated for
 * long-range periodicity.
 *
 * This assumes M_100 == M_010 == M_001 == 0.
 * This uses the quadrupole term and trace of the octupole terms only and
 * defaults to the monopole if the code is compiled with low-order gravity only.
 *
 * @param r_x x-component of the distance vector to the multipole.
 * @param r_y y-component of the distance vector to the multipole.
 * @param r_z z-component of the distance vector to the multipole.
 * @param r2 Square of the distance vector to the multipole.
 * @param h The softening length.
 * @param h_inv Inverse of the softening length.
 * @param r_s_inv The inverse of the gravity mesh-smoothing scale.
 * @param m The multipole.
 * @param f_x (return) The x-component of the acceleration.
 * @param f_y (return) The y-component of the acceleration.
 * @param f_z (return) The z-component of the acceleration.
 * @param pot (return) The potential.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pm_truncated(
    const float r_x, const float r_y, const float r_z, const float r2,
    const float h, const float h_inv, const float r_s_inv,
    const struct multipole *m, float *f_x, float *f_y, float *f_z, float *pot) {

/* In the case where the order is < 3, then there is only a monopole term left.
 * We can default to the normal P-P interaction with the mass of the multipole
 * and its CoM as the "particle" property */
#if SELF_GRAVITY_MULTIPOLE_ORDER < 3

  float f_ij, pot_ij;
  runner_iact_grav_pp_truncated(r2, h * h, h_inv, h_inv * h_inv * h_inv,
                                m->M_000, r_s_inv, &f_ij, &pot_ij);
  *f_x = f_ij * r_x;
  *f_y = f_ij * r_y;
  *f_z = f_ij * r_z;
  *pot = -pot_ij;

#else

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);

  /* Compute the derivatives of the potential */
  struct potential_derivatives_M2P d;
  compute_potential_derivatives_M2P(r_x, r_y, r_z, r2, r_inv, h, h_inv, 1,
                                    r_s_inv, &d);

  /* 0th order contributions */
  *f_x = m->M_000 * d.D_100;
  *f_y = m->M_000 * d.D_010;
  *f_z = m->M_000 * d.D_001;
  *pot = m->M_000 * d.D_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

/* 1st order contributions */

/* 1st order contributions are all 0 since the dipole is 0 */

/* *f_x = m->M_001 * d.D_101 + m->M_010 * d.D_110 + m->M_100 * d.D_200 ; */
/* *f_y = m->M_001 * d.D_011 + m->M_010 * d.D_020 + m->M_100 * d.D_110 ; */
/* *f_z = m->M_001 * d.D_002 + m->M_010 * d.D_011 + m->M_100 * d.D_101 ; */
/* *pot = m->M_001 * d.D_001 + m->M_010 * d.D_010 + m->M_100 * d.D_100 ; */

#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 2nd order contributions */
  *f_x += m->M_002 * d.D_102 + m->M_011 * d.D_111 + m->M_020 * d.D_120 +
          m->M_101 * d.D_201 + m->M_110 * d.D_210 + m->M_200 * d.D_300;
  *f_y += m->M_002 * d.D_012 + m->M_011 * d.D_021 + m->M_020 * d.D_030 +
          m->M_101 * d.D_111 + m->M_110 * d.D_120 + m->M_200 * d.D_210;
  *f_z += m->M_002 * d.D_003 + m->M_011 * d.D_012 + m->M_020 * d.D_021 +
          m->M_101 * d.D_102 + m->M_110 * d.D_111 + m->M_200 * d.D_201;
  *pot += m->M_002 * d.D_002 + m->M_011 * d.D_011 + m->M_020 * d.D_020 +
          m->M_101 * d.D_101 + m->M_110 * d.D_110 + m->M_200 * d.D_200;

#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 3rd order contributions */
  *f_x += m->M_003 * d.D_103 + m->M_012 * d.D_112 + m->M_021 * d.D_121 +
          m->M_030 * d.D_130 + m->M_102 * d.D_202 + m->M_111 * d.D_211 +
          m->M_120 * d.D_220 + m->M_201 * d.D_301 + m->M_210 * d.D_310 +
          m->M_300 * d.D_400;
  *f_y += m->M_003 * d.D_013 + m->M_012 * d.D_022 + m->M_021 * d.D_031 +
          m->M_030 * d.D_040 + m->M_102 * d.D_112 + m->M_111 * d.D_121 +
          m->M_120 * d.D_130 + m->M_201 * d.D_211 + m->M_210 * d.D_220 +
          m->M_300 * d.D_310;
  *f_z += m->M_003 * d.D_004 + m->M_012 * d.D_013 + m->M_021 * d.D_022 +
          m->M_030 * d.D_031 + m->M_102 * d.D_103 + m->M_111 * d.D_112 +
          m->M_120 * d.D_121 + m->M_201 * d.D_202 + m->M_210 * d.D_211 +
          m->M_300 * d.D_301;
  *pot += m->M_003 * d.D_003 + m->M_012 * d.D_012 + m->M_021 * d.D_021 +
          m->M_030 * d.D_030 + m->M_102 * d.D_102 + m->M_111 * d.D_111 +
          m->M_120 * d.D_120 + m->M_201 * d.D_201 + m->M_210 * d.D_210 +
          m->M_300 * d.D_300;

#endif

  /* Take care of the the sign convention */
  *f_x *= -1.f;
  *f_y *= -1.f;
  *f_z *= -1.f;
  *pot *= -1.f;
#endif
}

#endif /* SWIFT_POTENTIAL_GRAVITY_IACT_H */
