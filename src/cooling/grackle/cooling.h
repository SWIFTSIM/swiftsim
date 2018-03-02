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
#ifndef SWIFT_COOLING_GRACKLE_H
#define SWIFT_COOLING_GRACKLE_H

/**
 * @file src/cooling/none/cooling.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/* Some standard headers. */
#include <float.h>
#include <grackle.h>
#include <math.h>

/* Local includes. */
#include "../config.h"
#include "cooling_io.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3


/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct cooling_function_data* cooling) {

  xp->cooling_data.radiated_energy = 0.f;

  /* metal cooling = 1 */
  xp->cooling_data.metal_frac = cooling->chemistry.SolarMetalFractionByMass;

#if COOLING_GRACKLE_MODE >= 1
  gr_float zero = 1.e-20;

  /* primordial chemistry >= 1 */
  xp->cooling_data.HI_frac = grackle_data->HydrogenFractionByMass;
  xp->cooling_data.HII_frac = zero;
  xp->cooling_data.HeI_frac = 1. - grackle_data->HydrogenFractionByMass;
  xp->cooling_data.HeII_frac = zero;
  xp->cooling_data.HeIII_frac = zero;
  xp->cooling_data.e_frac = zero;

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  xp->cooling_data.HM_frac = zero;
  xp->cooling_data.H2I_frac = zero;
  xp->cooling_data.H2II_frac = zero;

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  xp->cooling_data.DI_frac = grackle_data->HydrogenFractionByMass *
                             grackle_data->DeuteriumToHydrogenRatio;
  xp->cooling_data.DII_frac = zero;
  xp->cooling_data.HDI_frac = zero;
#endif  // MODE >= 3

#endif  // MODE >= 2

#endif  // MODE >= 1
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
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'Grackle'.");
  message("Using Grackle    = %i", cooling->chemistry.use_grackle);
  message("Chemical network = %i",
          cooling->chemistry.primordial_chemistry);
  message("CloudyTable      = %s", cooling->cloudy_table);
  message("Redshift         = %g", cooling->redshift);
  message("UV background    = %d", cooling->with_uv_background);
  message("Metal cooling    = %i", cooling->chemistry.metal_cooling);
  message("Self Shielding   = %i", cooling->self_shielding_method);
  message("Specific Heating Rates   = %i",
	  cooling->provide_specific_heating_rates);
  message("Volumetric Heating Rates = %i",
	  cooling->provide_volumetric_heating_rates);
  message("Units:");
  message("\tComoving     = %i", cooling->units.comoving_coordinates);
  message("\tLength       = %g", cooling->units.length_units);
  message("\tDensity      = %g", cooling->units.density_units);
  message("\tTime         = %g", cooling->units.time_units);
  message("\tScale Factor = %g", cooling->units.a_units);
}

/**
 * @brief allocate required field
 *
 * @param data the #grackle_field_data
 */
__attribute__((always_inline)) INLINE static void cooling_malloc_data(
    grackle_field_data* data) {

#if COOLING_GRACKLE_MODE >= 1
  /* primordial chemistry >= 1 */
  data->HI_density = malloc(sizeof(gr_float));
  data->HII_density = malloc(sizeof(gr_float));
  data->HeI_density = malloc(sizeof(gr_float));
  data->HeII_density = malloc(sizeof(gr_float));
  data->HeIII_density = malloc(sizeof(gr_float));
  data->e_density = malloc(sizeof(gr_float));
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  data->HM_density = malloc(sizeof(gr_float));
  data->H2I_density = malloc(sizeof(gr_float));
  data->H2II_density = malloc(sizeof(gr_float));
#endif  // MODE 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  data->DI_density = malloc(sizeof(gr_float));
  data->DII_density = malloc(sizeof(gr_float));
  data->HDI_density = malloc(sizeof(gr_float));
#endif  // MODE >= 3

  /* metal cooling = 1 */
  data->metal_density = malloc(sizeof(gr_float));

  /* /\* volumetric heating rate *\/ */
  /* data->volumetric_heating_rate = NULL; */

  /* /\* specific heating rate *\/ */
  /* data->specific_heating_rate = NULL; */
}

/**
 * @brief free the allocated memory
 *
 * @param data the #grackle_field_data
 */

__attribute__((always_inline)) INLINE static void cooling_free_data(
    grackle_field_data* data) {

#if COOLING_GRACKLE_MODE >= 1
  /* primordial chemistry >= 1 */
  free(data->HI_density);
  free(data->HII_density);
  free(data->HeI_density);
  free(data->HeII_density);
  free(data->HeIII_density);
  free(data->e_density);
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  free(data->HM_density);
  free(data->H2I_density);
  free(data->H2II_density);
#endif  // MODE 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  free(data->DI_density);
  free(data->DII_density);
  free(data->HDI_density);
#endif  // MODE >= 3

  /* metal cooling = 1 */
  free(data->metal_density);

  /* /\* volumetric heating rate *\/ */
  /* data->volumetric_heating_rate = NULL; */

  /* /\* specific heating rate *\/ */
  /* data->specific_heating_rate = NULL; */
}

/**
 * @brief copy xp to data
 *
 * requires the particle to have been transformed into density
 *
 * @param data the #grackle_field_data
 * @param xp the #xpart
 */
__attribute__((always_inline)) INLINE static void cooling_copy_to_data(
    grackle_field_data* data, const struct xpart* xp, const gr_float rho) {

#if COOLING_GRACKLE_MODE >= 1
  /* primordial chemistry >= 1 */
  data->HI_density[0] = xp->cooling_data.HI_frac * rho;
  data->HII_density[0] = xp->cooling_data.HII_frac * rho;
  data->HeI_density[0] = xp->cooling_data.HeI_frac * rho;
  data->HeII_density[0] = xp->cooling_data.HeII_frac * rho;
  data->HeIII_density[0] = xp->cooling_data.HeIII_frac * rho;
  data->e_density[0] = xp->cooling_data.e_frac * rho;
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  data->HM_density[0] = xp->cooling_data.HM_frac * rho;
  data->H2I_density[0] = xp->cooling_data.H2I_frac * rho;
  data->H2II_density[0] = xp->cooling_data.H2II_frac * rho;
#endif  // MODE 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  data->DI_density[0] = xp->cooling_data.DI_frac * rho;
  data->DII_density[0] = xp->cooling_data.DII_frac * rho;
  data->HDI_density[0] = xp->cooling_data.HDI_frac * rho;
#endif  // MODE >= 3

  /* metal cooling = 1 */
  data->metal_density[0] = xp->cooling_data.metal_frac * rho;

  /* volumetric heating rate */
  data->volumetric_heating_rate = NULL;

  /* specific heating rate */
  data->specific_heating_rate = NULL;
}

/**
 * @brief copy data to xp
 *
 * @param data the #grackle_field_data
 * @param xp the #xpart
 */
__attribute__((always_inline)) INLINE static void cooling_copy_to_particle(
    const grackle_field_data* data, struct xpart* xp, const gr_float rho) {

#if COOLING_GRACKLE_MODE >= 1
  /* primordial chemistry >= 1 */
  xp->cooling_data.HI_frac = data->HI_density[0] / rho;
  xp->cooling_data.HII_frac = data->HII_density[0] / rho;
  xp->cooling_data.HeI_frac = data->HeI_density[0] / rho;
  xp->cooling_data.HeII_frac = data->HeII_density[0] / rho;
  xp->cooling_data.HeIII_frac = data->HeIII_density[0] / rho;
  xp->cooling_data.e_frac = data->e_density[0] / rho;
#endif  // MODE >= 1

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  xp->cooling_data.HM_frac = data->HM_density[0] / rho;
  xp->cooling_data.H2I_frac = data->H2I_density[0] / rho;
  xp->cooling_data.H2II_frac = data->H2II_density[0] / rho;
#endif  // MODE 2

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  xp->cooling_data.DI_frac = data->DI_density[0] / rho;
  xp->cooling_data.DII_frac = data->DII_density[0] / rho;
  xp->cooling_data.HDI_frac = data->HDI_density[0] / rho;
#endif  // MODE >= 3

  /* metal cooling = 1 */
  xp->cooling_data.metal_frac = data->metal_density[0] / rho;

  /* /\* volumetric heating rate *\/ */
  /* data->volumetric_heating_rate = NULL; */

  /* /\* specific heating rate *\/ */
  /* data->specific_heating_rate = NULL; */
}

/**
 * @brief Compute the cooling rate and update the particle chemistry data
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 *
 * @return du / dt
 */
__attribute__((always_inline)) INLINE static double cooling_rate(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, double dt) {

  /* set current time */
  code_units units = cooling->units;
  if (cooling->redshift == -1)
    error("TODO time dependant redshift");
  else
    units.a_value = 1. / (1. + cooling->redshift);

  /* initialize data */
  grackle_field_data data;

  /* set values */
  /* grid */
  int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  int grid_start[GRACKLE_RANK] = {0, 0, 0};
  int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};

  data.grid_rank = GRACKLE_RANK;
  data.grid_dimension = grid_dimension;
  data.grid_start = grid_start;
  data.grid_end = grid_end;

  /* general particle data */
  const gr_float energy_before = hydro_get_internal_energy(p);
  gr_float density = hydro_get_density(p);
  gr_float energy = energy_before;

  /* initialize density */
  data.density = &density;

  /* initialize energy */
  data.internal_energy = &energy;

  /* grackle 3.0 doc: "Currently not used" */
  data.x_velocity = NULL;
  data.y_velocity = NULL;
  data.z_velocity = NULL;

  /* allocate grackle data */
  cooling_malloc_data(&data);

  /* copy data from particle to grackle data */
  cooling_copy_to_data(&data, xp, density);

  /* solve chemistry with table */
  if (solve_chemistry(&units, &data, dt) == 0) {
    error("Error in solve_chemistry.");
  }

  /* copy from grackle data to particle */
  cooling_copy_to_particle(&data, xp, density);

  /* free allocated memory */
  cooling_free_data(&data);

  /* compute rate */
  return (energy - energy_before) / dt;
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
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, double dt) {

  if (dt == 0.) return;

  /* Current du_dt */
  const float hydro_du_dt = hydro_get_internal_energy_dt(p);

  /* compute cooling rate */
  const float du_dt = cooling_rate(phys_const, us, cooling, p, xp, dt);

  /* record energy lost */
  xp->cooling_data.radiated_energy += -du_dt * dt * hydro_get_mass(p);

  /* Update the internal energy */
  hydro_set_internal_energy_dt(p, hydro_du_dt + du_dt);
}

/**
 * @brief Computes the cooling time-step.
 *
 * We return FLT_MAX so as to impose no limit on the time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const struct part* restrict p) {

  return FLT_MAX;
}



/**
 * @brief Initialises the cooling unit system.
 *
 * @param us The current internal system of units.
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_units(
    const struct unit_system* us,
    struct cooling_function_data* cooling) {

  /* These are conversions from code units to cgs. */

  /* first cosmo */
  cooling->units.a_units = 1.0;  // units for the expansion factor (1/1+zi)
  cooling->units.a_value = 1.0;

  /* We assume here all physical quantities to
     be in proper coordinate (not comobile)  */
  cooling->units.comoving_coordinates = 0;

  /* then units */
  cooling->units.density_units =
      us->UnitMass_in_cgs / pow(us->UnitLength_in_cgs, 3);
  cooling->units.length_units = us->UnitLength_in_cgs;
  cooling->units.time_units = us->UnitTime_in_cgs;
  cooling->units.velocity_units = cooling->units.a_units *
                                  cooling->units.length_units /
                                  cooling->units.time_units;
}

/**
 * @brief Initialises Grackle.
 *
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_grackle(
    struct cooling_function_data* cooling) {
  
#ifdef SWIFT_DEBUG_CHECKS
  /* enable verbose for grackle */
  grackle_verbose = 1;
#endif

  chemistry_data* chemistry = &cooling->chemistry;

  /* Create a chemistry object for parameters and rate data. */
  if (set_default_chemistry_parameters(chemistry) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  // Set parameter values for chemistry.
  chemistry->use_grackle = 1;
  chemistry->with_radiative_cooling = 1;

  /* molecular network with H, He, D
   From Cloudy table */
  chemistry->primordial_chemistry = COOLING_GRACKLE_MODE;
  chemistry->metal_cooling = cooling->with_metal_cooling;
  chemistry->UVbackground = cooling->with_uv_background;
  chemistry->grackle_data_file = cooling->cloudy_table;

  /* radiative transfer */
  chemistry->use_radiative_transfer =
    cooling->provide_specific_heating_rates ||
    cooling->provide_volumetric_heating_rates;
  chemistry->use_volumetric_heating_rate =
    cooling->provide_volumetric_heating_rates;
  chemistry->use_specific_heating_rate =
    cooling->provide_specific_heating_rates;

  if (cooling->provide_specific_heating_rates &&
      cooling->provide_volumetric_heating_rates)
    message("WARNING: You should specified either the specific or the volumetric heating rates, not both");

  /* self shielding */
  chemistry->self_shielding_method =
    cooling->self_shielding_method;

  /* Initialize the chemistry object. */
  if (initialize_chemistry_data(&cooling->units) == 0) {
    error("Error in initialize_chemistry_data.");
  }

#ifdef SWIFT_DEBUG_CHECKS
  message("***************************************");
  message("initializing grackle cooling function");
  message("");
  cooling_print_backend(cooling);
  message("");
  message("***************************************");
#endif
  
}
  
/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_init_backend(
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {

  if (GRACKLE_NPART != 1)
    error("Grackle with multiple particles not implemented");

  /* read parameters */
  cooling_parse_arguments(parameter_file, cooling);

  /* Set up the units system. */
  cooling_init_units(us, cooling);

  cooling_init_grackle(cooling);

}

#endif /* SWIFT_COOLING_GRACKLE_H */
