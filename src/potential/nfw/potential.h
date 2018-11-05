/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018   Ashley Kelly ()
 *                      Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_POTENTIAL_NFW_H
#define SWIFT_POTENTIAL_NFW_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief External Potential Properties - NFW Potential
		rho(r) = rho_0 / ( (r/R_s)*(1+r/R_s)^2 )

	We however parameterise this in terms of c and virial_mass
 */
struct external_potential {

  /*! Position of the centre of potential */
  double x[3];

  /*! The density parameter for the potential rho_0 */
  double rho_0;

  /*! The scale radius of the NFW potential */
  double r_s;

  /*! The pre-factor 4*pi*G*rho_0*r_s*r_s*r_s */
  double pre_factor;

  /*! The critical density of the universe */
  double rho_c;

  /*! The concentration parameter */
  double c;

  /*! The virial mass */
  double virial_mass;

  /*! Time-step condition pre_factor, this factor is used to multiply times the 
   * orbital time, so in the case of 0.01 we take 1% of the orbital time as 
   * the time integration steps */
  double timestep_mult;

  /*! Minimum time step based on the orbital time at the softening times
   * the timestep_mult */
  double mintime;

  /*! the NFW lnthing */
  /* lnthing = ln(1+c) - c/(1+c) */
  double lnthing;

  double eps;

  double vel_20;
};

/**
 * @brief Computes the time-step due to the acceleration from the NFW potential
 *        as a fraction (timestep_mult) of the dynamical time of that particle.
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

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];
  
  const float r = sqrtf(dx*dx + dy*dy + dz*dz + potential->eps*potential->eps);

  const float mr = potential->virial_mass * ( log(1.f + r/potential->r_s) - 
       r/(r+potential->r_s) ) / potential->lnthing;

  const float period = 2 * M_PI * r * sqrtf(r) / sqrtf(phys_const->const_newton_G * mr);
  
  return max(potential->timestep_mult*period, potential->mintime);
}

/**
 * @brief Computes the gravitational acceleration from an NFW Halo potential.
 *
 * Note that the accelerations are multiplied by Newton's G constant
 * later on.
 *
 * a_x = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * x
 * a_y = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * y
 * a_z = 4 pi \rho_0 r_s^3 ( 1/((r+rs)*r^2) - log(1+r/rs)/r^3) * z
 *
 * @param time The current time.
 * @param potential The #external_potential used in the run.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static void external_gravity_acceleration(
    double time, const struct external_potential* restrict potential,
    const struct phys_const* restrict phys_const, struct gpart* restrict g) {

  const float dx = g->x[0] - potential->x[0];
  const float dy = g->x[1] - potential->x[1];
  const float dz = g->x[2] - potential->x[2];

  const float r = sqrtf(dx*dx + dy*dy + dz*dz + potential->eps*potential->eps);
  const float term1 = 1.0f*potential->pre_factor;
  const float term2 = (1.0f/((r + potential->r_s)*r*r) -
                       logf( 1.0f + r/potential->r_s )/(r*r*r));

  g->a_grav[0] += term1*term2*dx;
  g->a_grav[1] += term1*term2*dy;
  g->a_grav[2] += term1*term2*dz;
}

/**
 * @brief Computes the gravitational potential energy of a particle in an
 * NFW potential.
 *
 * phi = -4 * pi * G * rho_0 * r_s^3 * ln(1+r/r_s)
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

  const float r = sqrtf(dx*dx + dy*dy +dz*dz + potential->eps*potential->eps);
  const float term1 = -potential->pre_factor / r;
  const float term2 = logf(1.0f + r/potential->r_s);

  return term1*term2;
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
    struct swift_params* parameter_file,
    const struct phys_const* phys_const, const struct unit_system* us,
    const struct space* s, struct external_potential* potential) {

  const int useabspos = 
      parser_get_param_int(parameter_file, "NFWPotential:useabspos");
  parser_get_param_double_array(parameter_file, "NFWPotential:position",
                                 3, potential->x);   
  if ( useabspos == 0) {     
    potential->x[0] += s->dim[0] / 2.;     
    potential->x[1] += s->dim[1] / 2.;     
    potential->x[2] += s->dim[2] / 2.;   
  } 
  potential->timestep_mult = parser_get_param_float(
      parameter_file, "NFWPotential:timestep_mult");
  potential->c = parser_get_param_float(
      parameter_file, "NFWPotential:concentration");
	potential->virial_mass = parser_get_param_float(
	  	parameter_file, "NFWPotential:virial_mass");
	potential->rho_c = parser_get_param_float(
	  	parameter_file, "NFWPotential:critical_density");
	potential->eps = 0.05;

	double virial_radius = cbrtf( 3.0*potential->virial_mass/(800.0f*potential->rho_c*M_PI) );
	potential->r_s = virial_radius / potential->c;

	potential->rho_0 = (200.0f/3.0f) * potential->c * potential->c * potential->c * potential->rho_c /
			( logf(1.0f + potential->c) - potential->c/(1.0f + potential->c) );

	potential->pre_factor = 4.0f * M_PI * potential->rho_0
 												* potential->r_s * potential->r_s * potential->r_s;

 	potential->vel_20 = sqrtf(-potential->pre_factor * phys_const->const_newton_G *
      ( 1.0f/(20.0f + potential->r_s) - logf( 1.0f + 20.0f/potential->r_s )/20.0f ));
  potential->lnthing = log(1.f + potential->c) - potential->c/(1 + potential->c);
  const float sqrtgm = sqrtf( phys_const->const_newton_G * potential->virial_mass);
  const float epslnthing = log(1.f + potential->eps/potential->r_s) 
                            - potential->eps/(potential->eps + potential->r_s);
  
  potential->mintime = 2.f * M_PI * potential->eps * sqrtf(potential->eps) * 
              sqrtf(potential->lnthing / epslnthing) / sqrtgm * potential->timestep_mult;
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void potential_print_backend(
    const struct external_potential* potential) {

  message(
      "External potential is 'NFW' with properties are (x,y,z) = (%e, "
      "%e, %e), scale radius = %e "
      "timestep multiplier = %e, scale density = %e, vel at 20kpc = %e,"
      "mintime = %e",
      potential->x[0], potential->x[1], potential->x[2], potential->r_s,
      potential->timestep_mult, potential->rho_0, potential->vel_20, 
      potential->mintime);
}

#endif /* SWIFT_POTENTIAL_NFW_H */
