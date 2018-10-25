/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl) 
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
#ifndef SWIFT_POTENTIAL_HERNQUIST_H
#define SWIFT_POTENTIAL_ISOTHERMAL_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - Hernquist potential
 */
struct external_potential {

  /*! Position of the centre of potential */
  double x[3];

  /*! Mass of the potential */
  double M;

  /*! Scale length (often as a, to prevent confusion we use al) */
  double al;

  /*! Time-step condition pre-factor */
  double timestep_mult;
};

/**
 * @brief Computes the time-step due to the acceleration from an Hernquist 
 * potential
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float external_gravity_timestep(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const,
    const struct gpart* restrict g) {

  const float G_newton = phys_const->const_newton_G;

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  const float r2_sqrt = sqrtf(dx * dx + dy * dy + dz * dz);
  const float r_plus_a_inv = 1.f / ( r2_sqrt + potential->al );

  const float r2_plus_epsilon2_inv =
      1.f / (dx * dx + dy * dy + dz * dz + potential->epsilon2);
  const float drdv =
      dx * (g->v_full[0]) + dy * (g->v_full[1]) + dz * (g->v_full[2]);
  const double vrot = potential->vrot;

  const float dota_x = vrot * vrot * r2_plus_epsilon2_inv *
                       (g->v_full[0] - 2.f * drdv * dx * r2_plus_epsilon2_inv);
  const float dota_y = vrot * vrot * r2_plus_epsilon2_inv *
                       (g->v_full[1] - 2.f * drdv * dy * r2_plus_epsilon2_inv);
  const float dota_z = vrot * vrot * r2_plus_epsilon2_inv *
                       (g->v_full[2] - 2.f * drdv * dz * r2_plus_epsilon2_inv);
  const float dota_2 = dota_x * dota_x + dota_y * dota_y + dota_z * dota_z;
  const float a_2 = g->a_grav[0] * g->a_grav[0] + g->a_grav[1] * g->a_grav[1] +
                    g->a_grav[2] * g->a_grav[2];

  return potential->timestep_mult * sqrtf(a_2 / dota_2);
}

/**
 * @brief Computes the gravitational acceleration from an isothermal potential.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * a_x = -(v_rot^2 / G) * x / (r^2 + epsilon^2)
 * a_y = -(v_rot^2 / G) * y / (r^2 + epsilon^2)
 * a_z = -(v_rot^2 / G) * z / (r^2 + epsilon^2)
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, struct gpart* g) {

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];
  const float r2_plus_epsilon2_inv =
      1.f / (dx * dx + dy * dy + dz * dz + potential->epsilon2);

  const float term = -potential->vrot2_over_G * r2_plus_epsilon2_inv;

  g->a_grav[0] += term * dx;
  g->a_grav[1] += term * dy;
  g->a_grav[2] += term * dz;
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * isothermal potential.
 *
 * phi = -0.5 * vrot^2 * ln(r^2 + epsilon^2)
 *
 * @param time The current time (unused here).
 * @param potential The #external_potential used in the run.
 * @param phys_const Physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
external_gravity_get_potential_energy(
    double time, const struct external_potential* potential,
    const struct phys_const* const phys_const, const struct gpart* g) {

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  return 0.5f * potential->vrot * potential->vrot *
         logf(dx * dx + dy * dy + dz * dz + potential->epsilon2);
}

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
static INLINE void potential_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct space* s,
    struct external_potential* potential) {

  const int useabspos = 
      parser_get_param_int(parameter_file, "HernquistPotential:useabspos");
  parser_get_param_double_array(parameter_file, "HernquistPotential:position",
                                3, potential->x);
  if ( useabspos == 0) {
    potential->x[0] += s->dim[0] / 2.;
    potential->x[1] += s->dim[1] / 2.;
    potential->x[2] += s->dim[2] / 2.;
  } 

  potential->M = 
      parser_get_param_double(parameter_file, "HernquistPotential:Mass");
  potential->al = parser_get_param_double(parameter_file, 
                                          "Hernquistpotential:scalelength");
  potential->timestep_mult = parser_get_param_float(
      parameter_file, "HernquistPotential:timestep_mult");
}

/**
 * @brief prints the properties of the external potential to stdout.
 *
 * @param  potential the external potential properties.
 */
static inline void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "external potential is 'hernquist' with properties are (x,y,z) = (%e, "
      "%e, %e), mass = %e "
      "scale length = %e"
      "timestep multiplier = %e",
      potential->x[0], potential->x[1], potential->x[2], potential->M,
      potential->al, potential->timestep_mult;
}

#endif /* SWIFT_ISOTHERMAL_H */
