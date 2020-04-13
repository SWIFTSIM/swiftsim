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
#ifndef SWIFT_MULTI_SOFTENING_GRAVITY_IACT_H
#define SWIFT_MULTI_SOFTENING_GRAVITY_IACT_H

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
    const float mass, float *restrict f_ij, float *restrict pot_ij) {

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
    const float mass, const float r_s_inv, float *restrict f_ij,
    float *restrict pot_ij) {

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
    const float h, const float h_inv, const struct multipole *m,
    float *restrict f_x, float *restrict f_y, float *restrict f_z,
    float *restrict pot) {

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);

  /* Compute the derivatives of the potential */
  struct potential_derivatives_M2P d;
  potential_derivatives_compute_M2P(r_x, r_y, r_z, r2, r_inv, h,
                                    /*periodic=*/0, /*r_s_inv=*/0.f, &d);

  float F_000 = 0.f;
  float F_100 = 0.f;
  float F_010 = 0.f;
  float F_001 = 0.f;

  const float M_000 = m->M_000;
  const float D_000 = d.D_000;

  const float D_100 = d.D_100;
  const float D_010 = d.D_010;
  const float D_001 = d.D_001;

  /*  0th order term */
  F_000 -= M_000 * D_000;

  /*  1st order multipole term (addition to rank 1) */
  F_100 -= M_000 * D_100;
  F_010 -= M_000 * D_010;
  F_001 -= M_000 * D_001;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* The dipole term is zero when using the CoM */
  /* The compiler will optimize out the terms in the equations */
  /* below. We keep them written to maintain the logical structure. */
  const float M_100 = 0.f;
  const float M_010 = 0.f;
  const float M_001 = 0.f;

  const float D_200 = d.D_200;
  const float D_020 = d.D_020;
  const float D_002 = d.D_002;
  const float D_110 = d.D_110;
  const float D_101 = d.D_101;
  const float D_011 = d.D_011;

  /*  1st order multipole term (addition to rank 0) */
  F_000 += M_100 * D_100 + M_010 * D_010 + M_001 * D_001;

  /*  2nd order multipole term (addition to rank 1)*/
  F_100 += M_100 * D_200 + M_010 * D_110 + M_001 * D_101;
  F_010 += M_100 * D_110 + M_010 * D_020 + M_001 * D_011;
  F_001 += M_100 * D_101 + M_010 * D_011 + M_001 * D_002;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  const float M_200 = m->M_200;
  const float M_020 = m->M_020;
  const float M_002 = m->M_002;
  const float M_110 = m->M_110;
  const float M_101 = m->M_101;
  const float M_011 = m->M_011;

  const float D_300 = d.D_300;
  const float D_030 = d.D_030;
  const float D_003 = d.D_003;
  const float D_210 = d.D_210;
  const float D_201 = d.D_201;
  const float D_021 = d.D_021;
  const float D_120 = d.D_120;
  const float D_012 = d.D_012;
  const float D_102 = d.D_102;
  const float D_111 = d.D_111;

  /*  2nd order multipole term (addition to rank 0)*/
  F_000 -= M_200 * D_200 + M_020 * D_020 + M_002 * D_002;
  F_000 -= M_110 * D_110 + M_101 * D_101 + M_011 * D_011;

  /*  3rd order multipole term (addition to rank 1)*/
  F_100 -= M_200 * D_300 + M_020 * D_120 + M_002 * D_102;
  F_100 -= M_110 * D_210 + M_101 * D_201 + M_011 * D_111;
  F_010 -= M_200 * D_210 + M_020 * D_030 + M_002 * D_012;
  F_010 -= M_110 * D_120 + M_101 * D_111 + M_011 * D_021;
  F_001 -= M_200 * D_201 + M_020 * D_021 + M_002 * D_003;
  F_001 -= M_110 * D_111 + M_101 * D_102 + M_011 * D_012;

#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  const float M_300 = m->M_300;
  const float M_030 = m->M_030;
  const float M_003 = m->M_003;
  const float M_210 = m->M_210;
  const float M_201 = m->M_201;
  const float M_021 = m->M_021;
  const float M_120 = m->M_120;
  const float M_012 = m->M_012;
  const float M_102 = m->M_102;
  const float M_111 = m->M_111;

  const float D_400 = d.D_400;
  const float D_040 = d.D_040;
  const float D_004 = d.D_004;
  const float D_310 = d.D_310;
  const float D_301 = d.D_301;
  const float D_031 = d.D_031;
  const float D_130 = d.D_130;
  const float D_013 = d.D_013;
  const float D_103 = d.D_103;
  const float D_220 = d.D_220;
  const float D_202 = d.D_202;
  const float D_022 = d.D_022;
  const float D_211 = d.D_211;
  const float D_121 = d.D_121;
  const float D_112 = d.D_112;

  /*  3rd order multipole term (addition to rank 0)*/
  F_000 += M_300 * D_300 + M_030 * D_030 + M_003 * D_003;
  F_000 += M_210 * D_210 + M_201 * D_201 + M_120 * D_120;
  F_000 += M_021 * D_021 + M_102 * D_102 + M_012 * D_012;
  F_000 += M_111 * D_111;

  /* Compute 4th order field tensor terms (addition to rank 1) */
  F_001 += M_003 * D_004 + M_012 * D_013 + M_021 * D_022 + M_030 * D_031 +
           M_102 * D_103 + M_111 * D_112 + M_120 * D_121 + M_201 * D_202 +
           M_210 * D_211 + M_300 * D_301;
  F_010 += M_003 * D_013 + M_012 * D_022 + M_021 * D_031 + M_030 * D_040 +
           M_102 * D_112 + M_111 * D_121 + M_120 * D_130 + M_201 * D_211 +
           M_210 * D_220 + M_300 * D_310;
  F_100 += M_003 * D_103 + M_012 * D_112 + M_021 * D_121 + M_030 * D_130 +
           M_102 * D_202 + M_111 * D_211 + M_120 * D_220 + M_201 * D_301 +
           M_210 * D_310 + M_300 * D_400;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  const float M_400 = m->M_400;
  const float M_040 = m->M_040;
  const float M_004 = m->M_004;
  const float M_310 = m->M_310;
  const float M_301 = m->M_301;
  const float M_031 = m->M_031;
  const float M_130 = m->M_130;
  const float M_013 = m->M_013;
  const float M_103 = m->M_103;
  const float M_220 = m->M_220;
  const float M_202 = m->M_202;
  const float M_022 = m->M_022;
  const float M_211 = m->M_211;
  const float M_121 = m->M_121;
  const float M_112 = m->M_112;

  const float D_500 = d.D_500;
  const float D_050 = d.D_050;
  const float D_005 = d.D_005;
  const float D_410 = d.D_410;
  const float D_401 = d.D_401;
  const float D_041 = d.D_041;
  const float D_140 = d.D_140;
  const float D_014 = d.D_014;
  const float D_104 = d.D_104;
  const float D_320 = d.D_320;
  const float D_302 = d.D_302;
  const float D_230 = d.D_230;
  const float D_032 = d.D_032;
  const float D_203 = d.D_203;
  const float D_023 = d.D_023;
  const float D_122 = d.D_122;
  const float D_212 = d.D_212;
  const float D_221 = d.D_221;
  const float D_311 = d.D_311;
  const float D_131 = d.D_131;
  const float D_113 = d.D_113;

  /* Compute 4th order field tensor terms (addition to rank 0) */
  F_000 -= M_004 * D_004 + M_013 * D_013 + M_022 * D_022 + M_031 * D_031 +
           M_040 * D_040 + M_103 * D_103 + M_112 * D_112 + M_121 * D_121 +
           M_130 * D_130 + M_202 * D_202 + M_211 * D_211 + M_220 * D_220 +
           M_301 * D_301 + M_310 * D_310 + M_400 * D_400;

  /* Compute 5th order field tensor terms (addition to rank 1) */
  F_001 -= M_004 * D_005 + M_013 * D_014 + M_022 * D_023 + M_031 * D_032 +
           M_040 * D_041 + M_103 * D_104 + M_112 * D_113 + M_121 * D_122 +
           M_130 * D_131 + M_202 * D_203 + M_211 * D_212 + M_220 * D_221 +
           M_301 * D_302 + M_310 * D_311 + M_400 * D_401;
  F_010 -= M_004 * D_014 + M_013 * D_023 + M_022 * D_032 + M_031 * D_041 +
           M_040 * D_050 + M_103 * D_113 + M_112 * D_122 + M_121 * D_131 +
           M_130 * D_140 + M_202 * D_212 + M_211 * D_221 + M_220 * D_230 +
           M_301 * D_311 + M_310 * D_320 + M_400 * D_410;
  F_100 -= M_004 * D_104 + M_013 * D_113 + M_022 * D_122 + M_031 * D_131 +
           M_040 * D_140 + M_103 * D_203 + M_112 * D_212 + M_121 * D_221 +
           M_130 * D_230 + M_202 * D_302 + M_211 * D_311 + M_220 * D_320 +
           M_301 * D_401 + M_310 * D_410 + M_400 * D_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  const float M_500 = m->M_500;
  const float M_050 = m->M_050;
  const float M_005 = m->M_005;
  const float M_410 = m->M_410;
  const float M_401 = m->M_401;
  const float M_041 = m->M_041;
  const float M_140 = m->M_140;
  const float M_014 = m->M_014;
  const float M_104 = m->M_104;
  const float M_320 = m->M_320;
  const float M_302 = m->M_302;
  const float M_230 = m->M_230;
  const float M_032 = m->M_032;
  const float M_203 = m->M_203;
  const float M_023 = m->M_023;
  const float M_122 = m->M_122;
  const float M_212 = m->M_212;
  const float M_221 = m->M_221;
  const float M_311 = m->M_311;
  const float M_131 = m->M_131;
  const float M_113 = m->M_113;

  /* Compute 5th order field tensor terms (addition to rank 0) */
  F_000 += M_005 * D_005 + M_014 * D_014 + M_023 * D_023 + M_032 * D_032 +
           M_041 * D_041 + M_050 * D_050 + M_104 * D_104 + M_113 * D_113 +
           M_122 * D_122 + M_131 * D_131 + M_140 * D_140 + M_203 * D_203 +
           M_212 * D_212 + M_221 * D_221 + M_230 * D_230 + M_302 * D_302 +
           M_311 * D_311 + M_320 * D_320 + M_401 * D_401 + M_410 * D_410 +
           M_500 * D_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for orders >5"
#endif

  *f_x = F_100;
  *f_y = F_010;
  *f_z = F_001;
  *pot = F_000;
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
    const struct multipole *m, float *restrict f_x, float *restrict f_y,
    float *restrict f_z, float *restrict pot) {

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);

  /* Compute the derivatives of the potential */
  struct potential_derivatives_M2P d;
  potential_derivatives_compute_M2P(r_x, r_y, r_z, r2, r_inv, h, /*periodic=*/1,
                                    r_s_inv, &d);

  float F_000 = 0.f;
  float F_100 = 0.f;
  float F_010 = 0.f;
  float F_001 = 0.f;

  const float M_000 = m->M_000;
  const float D_000 = d.D_000;

  const float D_100 = d.D_100;
  const float D_010 = d.D_010;
  const float D_001 = d.D_001;

  /*  0th order term */
  F_000 -= M_000 * D_000;

  /*  1st order multipole term (addition to rank 1) */
  F_100 -= M_000 * D_100;
  F_010 -= M_000 * D_010;
  F_001 -= M_000 * D_001;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* The dipole term is zero when using the CoM */
  /* The compiler will optimize out the terms in the equations */
  /* below. We keep them written to maintain the logical structure. */
  const float M_100 = 0.f;
  const float M_010 = 0.f;
  const float M_001 = 0.f;

  const float D_200 = d.D_200;
  const float D_020 = d.D_020;
  const float D_002 = d.D_002;
  const float D_110 = d.D_110;
  const float D_101 = d.D_101;
  const float D_011 = d.D_011;

  /*  1st order multipole term (addition to rank 0) */
  F_000 += M_100 * D_100 + M_010 * D_010 + M_001 * D_001;

  /*  2nd order multipole term (addition to rank 1)*/
  F_100 += M_100 * D_200 + M_010 * D_110 + M_001 * D_101;
  F_010 += M_100 * D_110 + M_010 * D_020 + M_001 * D_011;
  F_001 += M_100 * D_101 + M_010 * D_011 + M_001 * D_002;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  const float M_200 = m->M_200;
  const float M_020 = m->M_020;
  const float M_002 = m->M_002;
  const float M_110 = m->M_110;
  const float M_101 = m->M_101;
  const float M_011 = m->M_011;

  const float D_300 = d.D_300;
  const float D_030 = d.D_030;
  const float D_003 = d.D_003;
  const float D_210 = d.D_210;
  const float D_201 = d.D_201;
  const float D_021 = d.D_021;
  const float D_120 = d.D_120;
  const float D_012 = d.D_012;
  const float D_102 = d.D_102;
  const float D_111 = d.D_111;

  /*  2nd order multipole term (addition to rank 0)*/
  F_000 -= M_200 * D_200 + M_020 * D_020 + M_002 * D_002;
  F_000 -= M_110 * D_110 + M_101 * D_101 + M_011 * D_011;

  /*  3rd order multipole term (addition to rank 1)*/
  F_100 -= M_200 * D_300 + M_020 * D_120 + M_002 * D_102;
  F_100 -= M_110 * D_210 + M_101 * D_201 + M_011 * D_111;
  F_010 -= M_200 * D_210 + M_020 * D_030 + M_002 * D_012;
  F_010 -= M_110 * D_120 + M_101 * D_111 + M_011 * D_021;
  F_001 -= M_200 * D_201 + M_020 * D_021 + M_002 * D_003;
  F_001 -= M_110 * D_111 + M_101 * D_102 + M_011 * D_012;

#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  const float M_300 = m->M_300;
  const float M_030 = m->M_030;
  const float M_003 = m->M_003;
  const float M_210 = m->M_210;
  const float M_201 = m->M_201;
  const float M_021 = m->M_021;
  const float M_120 = m->M_120;
  const float M_012 = m->M_012;
  const float M_102 = m->M_102;
  const float M_111 = m->M_111;

  const float D_400 = d.D_400;
  const float D_040 = d.D_040;
  const float D_004 = d.D_004;
  const float D_310 = d.D_310;
  const float D_301 = d.D_301;
  const float D_031 = d.D_031;
  const float D_130 = d.D_130;
  const float D_013 = d.D_013;
  const float D_103 = d.D_103;
  const float D_220 = d.D_220;
  const float D_202 = d.D_202;
  const float D_022 = d.D_022;
  const float D_211 = d.D_211;
  const float D_121 = d.D_121;
  const float D_112 = d.D_112;

  /*  3rd order multipole term (addition to rank 0)*/
  F_000 += M_300 * D_300 + M_030 * D_030 + M_003 * D_003;
  F_000 += M_210 * D_210 + M_201 * D_201 + M_120 * D_120;
  F_000 += M_021 * D_021 + M_102 * D_102 + M_012 * D_012;
  F_000 += M_111 * D_111;

  /* Compute 4th order field tensor terms (addition to rank 1) */
  F_001 += M_003 * D_004 + M_012 * D_013 + M_021 * D_022 + M_030 * D_031 +
           M_102 * D_103 + M_111 * D_112 + M_120 * D_121 + M_201 * D_202 +
           M_210 * D_211 + M_300 * D_301;
  F_010 += M_003 * D_013 + M_012 * D_022 + M_021 * D_031 + M_030 * D_040 +
           M_102 * D_112 + M_111 * D_121 + M_120 * D_130 + M_201 * D_211 +
           M_210 * D_220 + M_300 * D_310;
  F_100 += M_003 * D_103 + M_012 * D_112 + M_021 * D_121 + M_030 * D_130 +
           M_102 * D_202 + M_111 * D_211 + M_120 * D_220 + M_201 * D_301 +
           M_210 * D_310 + M_300 * D_400;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  const float M_400 = m->M_400;
  const float M_040 = m->M_040;
  const float M_004 = m->M_004;
  const float M_310 = m->M_310;
  const float M_301 = m->M_301;
  const float M_031 = m->M_031;
  const float M_130 = m->M_130;
  const float M_013 = m->M_013;
  const float M_103 = m->M_103;
  const float M_220 = m->M_220;
  const float M_202 = m->M_202;
  const float M_022 = m->M_022;
  const float M_211 = m->M_211;
  const float M_121 = m->M_121;
  const float M_112 = m->M_112;

  const float D_500 = d.D_500;
  const float D_050 = d.D_050;
  const float D_005 = d.D_005;
  const float D_410 = d.D_410;
  const float D_401 = d.D_401;
  const float D_041 = d.D_041;
  const float D_140 = d.D_140;
  const float D_014 = d.D_014;
  const float D_104 = d.D_104;
  const float D_320 = d.D_320;
  const float D_302 = d.D_302;
  const float D_230 = d.D_230;
  const float D_032 = d.D_032;
  const float D_203 = d.D_203;
  const float D_023 = d.D_023;
  const float D_122 = d.D_122;
  const float D_212 = d.D_212;
  const float D_221 = d.D_221;
  const float D_311 = d.D_311;
  const float D_131 = d.D_131;
  const float D_113 = d.D_113;

  /* Compute 4th order field tensor terms (addition to rank 0) */
  F_000 -= M_004 * D_004 + M_013 * D_013 + M_022 * D_022 + M_031 * D_031 +
           M_040 * D_040 + M_103 * D_103 + M_112 * D_112 + M_121 * D_121 +
           M_130 * D_130 + M_202 * D_202 + M_211 * D_211 + M_220 * D_220 +
           M_301 * D_301 + M_310 * D_310 + M_400 * D_400;

  /* Compute 5th order field tensor terms (addition to rank 1) */
  F_001 -= M_004 * D_005 + M_013 * D_014 + M_022 * D_023 + M_031 * D_032 +
           M_040 * D_041 + M_103 * D_104 + M_112 * D_113 + M_121 * D_122 +
           M_130 * D_131 + M_202 * D_203 + M_211 * D_212 + M_220 * D_221 +
           M_301 * D_302 + M_310 * D_311 + M_400 * D_401;
  F_010 -= M_004 * D_014 + M_013 * D_023 + M_022 * D_032 + M_031 * D_041 +
           M_040 * D_050 + M_103 * D_113 + M_112 * D_122 + M_121 * D_131 +
           M_130 * D_140 + M_202 * D_212 + M_211 * D_221 + M_220 * D_230 +
           M_301 * D_311 + M_310 * D_320 + M_400 * D_410;
  F_100 -= M_004 * D_104 + M_013 * D_113 + M_022 * D_122 + M_031 * D_131 +
           M_040 * D_140 + M_103 * D_203 + M_112 * D_212 + M_121 * D_221 +
           M_130 * D_230 + M_202 * D_302 + M_211 * D_311 + M_220 * D_320 +
           M_301 * D_401 + M_310 * D_410 + M_400 * D_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  const float M_500 = m->M_500;
  const float M_050 = m->M_050;
  const float M_005 = m->M_005;
  const float M_410 = m->M_410;
  const float M_401 = m->M_401;
  const float M_041 = m->M_041;
  const float M_140 = m->M_140;
  const float M_014 = m->M_014;
  const float M_104 = m->M_104;
  const float M_320 = m->M_320;
  const float M_302 = m->M_302;
  const float M_230 = m->M_230;
  const float M_032 = m->M_032;
  const float M_203 = m->M_203;
  const float M_023 = m->M_023;
  const float M_122 = m->M_122;
  const float M_212 = m->M_212;
  const float M_221 = m->M_221;
  const float M_311 = m->M_311;
  const float M_131 = m->M_131;
  const float M_113 = m->M_113;

  /* Compute 5th order field tensor terms (addition to rank 0) */
  F_000 += M_005 * D_005 + M_014 * D_014 + M_023 * D_023 + M_032 * D_032 +
           M_041 * D_041 + M_050 * D_050 + M_104 * D_104 + M_113 * D_113 +
           M_122 * D_122 + M_131 * D_131 + M_140 * D_140 + M_203 * D_203 +
           M_212 * D_212 + M_221 * D_221 + M_230 * D_230 + M_302 * D_302 +
           M_311 * D_311 + M_320 * D_320 + M_401 * D_401 + M_410 * D_410 +
           M_500 * D_500;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for orders >5"
#endif

  *f_x = F_100;
  *f_y = F_010;
  *f_z = F_001;
  *pot = F_000;
}

#endif /* SWIFT_MULTI_SOFTENING_GRAVITY_IACT_H */
