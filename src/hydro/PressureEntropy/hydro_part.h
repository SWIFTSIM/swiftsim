/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENTROPY_HYDRO_PART_H
#define SWIFT_PRESSURE_ENTROPY_HYDRO_PART_H

/**
 * @file PressureEntropy/hydro_part.h
 * @brief Pressure-Entropy implementation of SPH (Particle definition)
 *
 * The thermal variable is the entropy (S) and the entropy is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows eqautions (19), (21) and (22) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

#include "../config.h"
#include "cooling_struct.h"

/* Extra particle data not needed during the SPH loops over neighbours. */
struct xpart {

  /*! Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Velocity at the last full step. */
  float v_full[3];

  /*! Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

#ifdef SWIFT_STREAMING_IO
  /*! For logging, last values of these quantities recorded. */
  struct {
    double x[3];
    float v[3];
    float h;
    float u;
    float rho;
    int t_x, t_v, t_h, t_u, t_rho;
  } last_logged;
#endif

} SWIFT_STRUCT_ALIGN;

/* Data of a single particle. */
struct part {

  /*! Particle position. */
  double x[3];

  /*! Particle predicted velocity. */
  float v[3];

  /*! Particle acceleration. */
  float a_hydro[3];

  /*! Particle cutoff radius. */
  float h;

  /*! Particle mass. */
  float mass;

  /*! Particle time of beginning of time-step. */
  int ti_begin;

  /*! Particle time of end of time-step. */
  int ti_end;

  /*! Particle density. */
  float rho;

  /*! Particle weighted density */
  float rho_bar;

  /*! Particle entropy. */
  float entropy;

  /*! Entropy time derivative */
  float entropy_dt;

  /*! Particle entropy to the power 1/gamma. */
  float entropy_one_over_gamma;

  union {

    struct {

      /*! Number of neighbours. */
      float wcount;

      /*! Number of neighbours spatial derivative. */
      float wcount_dh;

      /*! Derivative of density with respect to h */
      float rho_dh;

      /*! Derivative of pressure with respect to h */
      float pressure_dh;

      /*! Particle velocity curl. */
      float rot_v[3];

      /*! Particle velocity divergence. */
      float div_v;

    } density;

    struct {

      /*! Balsara switch */
      float balsara;

      /*! "Grad h" term */
      float f;

      /*! Pressure term */
      float P_over_rho2;

      /*! Particle sound speed. */
      float soundspeed;

      /*! Signal velocity. */
      float v_sig;

      /*! Time derivative of the smoothing length */
      float h_dt;

    } force;
  };

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_PART_H */
