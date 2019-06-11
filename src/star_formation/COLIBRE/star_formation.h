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
 *******************************************************************************/
#ifndef SWIFT_COLIBRE_STAR_FORMATION_H
#define SWIFT_COLIBRE_STAR_FORMATION_H

/* Local includes */
#include "adiabatic_index.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "exp10.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "random.h"
#include "stars.h"
#include "units.h"

/**
 * @file src/star_formation/COLIBRE/star_formation.h
 * @brief Star formation model used in the COLIBRE model
 */

/**
 * @brief Properties of the COLIBRE star formation model.
 */
struct star_formation {

  /*! Critical overdensity */
  double min_over_den;

  /*! Star formation efficiency */
  double sfe;

  /*! free fall time constant sqrt(3 pi/ 32G) (cm^-1.5 g^.5 s^1) */
  double ff_const;

  /*! star formation efficiency over free fall time constant (cm^1.5 g^-.5 s^-1)
   */
  double mdot_const;

  /*! Dalla Vecchia & Schaye temperature criteria */
  double temperature_margin_threshold_dex;

  /*! 10^Tdex of Dalla Vecchia & SChaye temperature criteria */
  double ten_to_temperature_margin_threshold_dex;

  /*! Absolute temperature criteria */
  double temperature_threshold;

  /*! Maximal physical density, particles with a higher density instantaneously
   * turn into stars */
  double maximal_density_HpCM3;

  /*! Maximal physical density (internal units), particles with a higher density
   * instantaneously turn into stars */
  double maximal_density;
};

/**
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 *
 */
INLINE static int star_formation_is_star_forming(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const struct entropy_floor_properties* restrict entropy_props) {

  /* Minimal density (converted from critical density) for star formation */
  const double rho_crit_times_min_over_den =
      cosmo->critical_density * starform->min_over_den;

  /* Physical density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Deside whether we should form stars or not,
   * first we deterime if we have the correct over density
   * if that is true we calculate if either the maximum density
   * threshold is reached or if the metallicity dependent
   * threshold is reached, after this we calculate if the
   * temperature is appropriate */
  if (physical_density < rho_crit_times_min_over_den) return 0;

  /* Did we reach the maximal physical density in which we need to be turned
   * instantaneously into a star?*/
  if (physical_density > starform->maximal_density) return 1;

  /* Calculate the temperature */
  const double temperature = cooling_get_temperature(phys_const, hydro_props,
                                                     us, cosmo, cooling, p, xp);

  /* Check if the temperature criteria is satisfied */
  if (temperature < starform->temperature_threshold) return 1;

  /* Get the EOS temperature from the entropy floor */
  const double temperature_eos =
      entropy_floor_temperature(p, cosmo, entropy_props);

  /* Check if the constant temperature is satisfied */
  return (temperature <
          temperature_eos * starform->ten_to_temperature_margin_threshold_dex);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart. The star formation is calculated as a simple
 * Schmidt law with an efficiency per free-fall time that can be specified,
 * the free-fall time is based on the total SPH density.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const double dt_star) {

  /* Abort early if time-step size is 0 */
  if (dt_star == 0.) {

    xp->sf_data.SFR = 0.f;
    return;
  }

  /* Mass density of this particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* In the case the particle is immediatly converted to a star particle */
  if (physical_density > starform->maximal_density) {
    xp->sf_data.SFR = hydro_get_mass(p) / dt_star;
  }

  /* Calculate the SFR per gas mass */
  double SFRpergasmass = starform->mdot_const * sqrt(physical_density);

  /* Store the SFR */
  xp->sf_data.SFR = SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.
 *
 * Equation 21 of Schaye & Dalla Vecchia 2008.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct engine* e,
    const double dt_star) {

  /* Calculate the propability of forming a star */
  const double prob = xp->sf_data.SFR * dt_star / hydro_get_mass(p);

  /* Get a unique random number between 0 and 1 for star formation */
  const double random_number =
      random_unit_interval(p->id, e->ti_current, random_number_star_formation);

  /* Have we been lucky and need to form a star? */
  return (prob > random_number);
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct star_formation* starform, const int with_cosmology) {

  /* Check if it is the first time steps after star formation */
  if (xp->sf_data.SFR > 0.f) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      xp->sf_data.SFR = -e->cosmology->a;
    } else {
      xp->sf_data.SFR = -e->time;
    }
  }
}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 * @param starform the star formation law properties to use.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
INLINE static void star_formation_copy_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* Store the current mass */
  sp->mass = hydro_get_mass(p);

  /* Store the current mass as the initial mass */
  sp->mass_init = hydro_get_mass(p);

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Store the chemistry struct in the star particle */
  sp->chemistry_data = p->chemistry_data;

  /* Store the tracers data */
  sp->tracers_data = xp->tracers_data;

  /* Store the birth density in the star particle */
  sp->birth_density = hydro_get_physical_density(p, cosmo);

  /* Store the birth temperature in the star particle */
  sp->birth_temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                  cosmo, cooling, p, xp);

  /* Flag that this particle has not done feedback yet */
  sp->f_E = -1.f;
}

/**
 * @brief initialization of the star formation law
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units.
 * @param hydro_props The propertis of the hydro model.
 * @param starform the star formation law properties to initialize
 */
INLINE static void starformation_init_backend(
    struct swift_params* parameter_file, const struct phys_const* phys_const,
    const struct unit_system* us, const struct hydro_props* hydro_props,
    struct star_formation* starform) {

  /* Get the Gravitational constant */
  const double G_newton = phys_const->const_newton_G;

  /* Read the critical density contrast from the parameter file*/
  starform->min_over_den = parser_get_param_double(
      parameter_file, "COLIBREStarFormation:min_over_density");

  /* Get the star formation efficiency */
  starform->sfe = parser_get_param_double(
      parameter_file, "COLIBREStarFormation:star_formation_efficiency");

  /* Calculate the ff constant */
  starform->ff_const = sqrt(3.0 * M_PI / (32.0 * G_newton));

  /* Calculate the constant */
  starform->mdot_const = starform->sfe / starform->ff_const;

  /* Get the temperature threshold */
  starform->temperature_margin_threshold_dex = parser_get_opt_param_double(
      parameter_file, "COLIBREStarFormation:EOS_temperature_margin_dex", 10.f);

  /* Write it as a constant temperature multiplication */
  starform->ten_to_temperature_margin_threshold_dex =
      exp10(starform->temperature_margin_threshold_dex);

  /* Get the temperature threshold */
  starform->temperature_threshold = parser_get_opt_param_double(
      parameter_file, "COLIBREStarFormation:temperature_threshold_K", FLT_MAX);

  starform->maximal_density_HpCM3 = parser_get_opt_param_double(
      parameter_file, "COLIBREStarFormation:threshold_max_density_H_p_cm3",
      FLT_MAX);

  /* Conversion of number density from cgs */
  const double number_density_from_cgs =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Calculate the maximal physical density to a mass density using the mean
   * molecular weight of primordial gas (mu=1.22)*/
  starform->maximal_density = starform->maximal_density_HpCM3 *
                              phys_const->const_proton_mass *
                              number_density_from_cgs * hydro_props->mu_neutral;
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  /* Print the star formation properties */
  message("Star formation law is COLIBRE");
  message(
      "With properties: Star formation efficiency = %e, temperature "
      "threshold = %e, minimum over density = %e T dex = %e maximal "
      "density = %e",
      starform->sfe, starform->temperature_threshold, starform->min_over_den,
      starform->temperature_margin_threshold_dex,
      starform->maximal_density_HpCM3);
}

/**
 * @brief Finishes the density calculation.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the COLIBRE star formation model.
 *
 * @param p The particle to act upon
 * @param cd The global star_formation information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_formation_end_density(
    struct part* restrict p, const struct star_formation* cd,
    const struct cosmology* cosmo) {}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the COLIBRE star formation model.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(struct part* restrict p,
                                      struct xpart* restrict xp,
                                      const struct star_formation* cd,
                                      const struct cosmology* cosmo) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_part(const struct phys_const* restrict phys_const,
                               const struct unit_system* restrict us,
                               const struct cosmology* restrict cosmo,
                               const struct star_formation* data,
                               const struct part* restrict p,
                               struct xpart* restrict xp) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state.
 *
 * Nothing to do here. We do not need to compute any quantity in the hydro
 * density loop for the COLIBRE star formation model.
 *
 * @param p Pointer to the particle data.
 * @param data The global star_formation information.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    struct part* restrict p, const struct star_formation* data) {}

#endif /* SWIFT_COLIBRE_STAR_FORMATION_H */
