/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Richard Bower (r.g.bower@durham.ac.uk)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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

#ifndef SWIFT_COOLING_CONST_LAMBDA_H
#define SWIFT_COOLING_CONST_LAMBDA_H

/* Some standard headers. */
#include <math.h>
#include <float.h>

/* Local includes. */
#include "const.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"


/**
 * @brief Calculates du/dt in code units for a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data..
 */
__attribute__((always_inline)) INLINE static float cooling_rate(
    const struct phys_const* const phys_const, const struct UnitSystem* us,
    const struct cooling_function_data* cooling, const struct part* p) {

  /* Get particle density */
  const float rho = p->rho;
  
  /* Get cooling function properties */
  const float X_H = cooling->hydrogen_mass_abundance;

  /* Calculate du_dt */
  const float du_dt = -cooling->lambda * ( X_H * rho / phys_const->const_proton_mass) * 
                                         ( X_H * rho / phys_const->const_proton_mass) / rho;
  return du_dt;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct UnitSystem* restrict us,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {

  /* Get current internal energy (dt=0) */
  const float u_old = hydro_get_internal_energy(p, 0.f);
  
  /* Internal energy floor */
  const float u_floor = cooling->min_energy;

  /* Calculate du_dt */
  const float du_dt = cooling_rate(phys_const, us, cooling, p);

  /* Intergrate cooling equation, but enforce energy floor */
  float u_new;
  if (u_old + du_dt * dt > u_floor) {
    u_new = u_old + du_dt * dt;
  } else {
    u_new = u_floor;
  }
  /* Update the internal energy */
  hydro_set_internal_energy(p, u_new);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy += hydro_get_mass(p) * (u_old - u_new);
}

/**
 * @brief Computes the time-step due to cooling
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct UnitSystem* restrict us, const struct part* restrict p) {

  /* Get current internal energy (dt=0) */
  const float u = hydro_get_internal_energy(p, 0.f);

  /* If we are close to the energy floor, ignore cooling timestep condition*/

  /* if (abs((u - cooling->min_energy)/cooling->min_energy) < 1.e-4){ */
  /*   return FLT_MAX; */
  /* } */
  /* else{ */
  /*   /\* Get du_dt *\/ */
  /*   const float du_dt = cooling_rate(phys_const, us, cooling, p); */
  /*   return u / abs(du_dt); */
  /* } */
  const float du_dt = cooling_rate(phys_const, us, cooling, p);

  return u/abs(du_dt);
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void cooling_init_part(
    const struct part* restrict p, struct xpart* restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cooling The cooling properties to initialize
 */
static INLINE void cooling_init_backend(
    const struct swift_params* parameter_file, const struct UnitSystem* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {

  const double lambda_cgs =
      parser_get_param_double(parameter_file, "LambdaCooling:lambda_cgs");
  const float min_temperature = parser_get_param_double(
      parameter_file, "LambdaCooling:minimum_temperature");
  cooling->hydrogen_mass_abundance = parser_get_param_double(
      parameter_file, "LambdaCooling:hydrogen_mass_abundance");
  cooling->mean_molecular_weight = parser_get_param_double(
      parameter_file, "LambdaCooling:mean_molecular_weight");
  cooling->cooling_tstep_mult = parser_get_param_double(
      parameter_file, "LambdaCooling:cooling_tstep_mult");

  /* convert minimum temperature into minimum internal energy */
  const float u_floor =
      phys_const->const_boltzmann_k * min_temperature /
      (hydro_gamma_minus_one * cooling->mean_molecular_weight *
       phys_const->const_proton_mass);

  cooling->min_energy = u_floor;

  /* convert lambda to code units */
  cooling->lambda = lambda_cgs *
                    units_cgs_conversion_factor(us,UNIT_CONV_TIME) / 
                    (units_cgs_conversion_factor(us,UNIT_CONV_ENERGY) * 
		    units_cgs_conversion_factor(us,UNIT_CONV_VOLUME));
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message(
      "Cooling function is 'Constant lambda' with "
      "(lambda,min_energy,hydrogen_mass_abundance,mean_molecular_weight) "
      "=  (%g,%g,%g,%g)",
      cooling->lambda, cooling->min_energy,
      cooling->hydrogen_mass_abundance, cooling->mean_molecular_weight);
}

#endif /* SWIFT_COOLING_CONST_LAMBDA_H */
