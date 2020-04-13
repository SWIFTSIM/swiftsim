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
//#include "stars.h"
#include "units.h"

/**
 * @file src/star_formation/COLIBRE/star_formation.h
 * @brief Star formation model used in the COLIBRE model
 */

#define star_formation_need_update_dx_max 0

/**
 * @file src/star_formation/COLIBRE/star_formation.h
 * @brief Star formation model used in the COLIBRE model
 */

/**
 * @brief Functional form of the star formation law
 */
enum star_formation_law {
  colibre_star_formation_schmidt_law, /*< Schmidt law */
  colibre_star_formation_pressure_law /*< Pressure law */
};

/**
 * @brief Properties of the COLIBRE star formation model.
 */
struct star_formation {

  /* Starting with the star forming gas criteria */

  /*! Critical overdensity */
  double min_over_den;

  /*! Maximal physical density, particles with a higher density instantaneously
   * turn into stars */
  double maximal_density_HpCM3;

  /*! Maximal physical density (internal units), particles with a higher density
   * instantaneously turn into stars */
  double maximal_density;

  /*! Virial constant used in the calculation (no units)*/
  double alpha_virial;

  /*! Inverse of the virial constant including the constants (G) used in the
   * calculation (no units)*/
  double alpha_virial_inv;

  /*! subgrid density threshold in user units */
  double subgrid_density_threshold_HpCM3;

  /*! subgrid density threshold in internal units */
  double subgrid_density_threshold;

  /*! Absolute subgrid temperature threshold */
  double temperature_threshold;

  /* Model dependent parameters */

  /* Which model are we using */
  enum star_formation_law SF_law;

  /* The Schmidt model parameters */
  struct {

    /*! Star formation efficiency */
    double sfe;

    /*! star formation efficiency over free fall time constant (cm^1.5 g^-.5
     * s^-1) */
    double mdot_const;

  } schmidt_law;

  /* The pressure model parameters */
  struct {

    /*! gas fraction */
    double fgas;

    /*! Slope of the KS law */
    double KS_power_law;

    /*! Star formation law slope */
    double SF_power_law;

    /*! Normalization of the KS star formation law (Msun / kpc^2 / yr) */
    double KS_normalization_MSUNpYRpKPC2;

    /*! Normalization of the KS star formation law (internal units) */
    double KS_normalization;

    /*! star formation normalization (internal units) */
    double SF_normalization;

  } pressure_law;
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
    const struct part* p, const struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* entropy_props) {

  /* Physical density of the particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Did we reach the maximal physical density in which we need to be turned
   * instantaneously into a star?*/
  if (physical_density > starform->maximal_density) return 1;

  /* Minimal density (converted from critical density) for star formation */
  const double rho_crit_times_min_over_den =
      cosmo->critical_density * starform->min_over_den;

  /* Deside whether we should form stars or not,
   * first we deterime if we have the correct over density
   * if that is true we calculate if either the maximum density
   * threshold is reached after this we calculate if the
   * temperature is appropriate */
  if (physical_density < rho_crit_times_min_over_den) return 0;

  /* Get the subgrid density */
  const double subgrid_density = p->cooling_data.subgrid_dens;

  /* Check if the subgrid density criterion is satisfied */
  if (subgrid_density > starform->subgrid_density_threshold) return 1;

  /* Calculate the temperature */
  const double temperature = cooling_get_temperature(phys_const, hydro_props,
                                                     us, cosmo, cooling, p, xp);

  /* Get the subgrid temperature from the tracers */
  const double subgrid_temperature = p->cooling_data.subgrid_temp;

  /* Check if we satisfy the subgrid temperature criterion */
  if (subgrid_temperature < starform->temperature_threshold) return 1;

  /* Calculate the thermal velocity dispersion */
  const double sigma_thermal_squared =
      3. * hydro_get_physical_pressure(p, cosmo) / physical_density;

  /* Get the subgrid thermal velocity dispersion */
  const double sigma_thermal2_subgrid =
      sigma_thermal_squared * subgrid_temperature / temperature;

  /* Get the gas mass */
  const double gas_mass = hydro_get_mass(p);

  /* Get the gas_mass to the power 2/3 */
  const double gas_mass_to_1_3th = cbrt(gas_mass);
  const double gas_mass_to_2_3th = gas_mass_to_1_3th * gas_mass_to_1_3th;

  /* Get the subgrid density to the power 1/3 */
  const double subgrid_density_to_1_3th = cbrt(subgrid_density);

  /* Check if we satisfy the subgrid virial criterion by calculation
   * alpha/alpha_virial, in which starform->alpha_virial_inv already includes
   * the gravitational constant */
  const double alpha = starform->alpha_virial_inv *
                       (p->sf_data.sigma_v2 + sigma_thermal2_subgrid) /
                       (subgrid_density_to_1_3th * gas_mass_to_2_3th);

  /* Check if the virial critiria is below one */
  return (alpha < 1.);
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
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR_schmidt_law(
    const struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

  /* Mass density of this particle */
  const double physical_density = hydro_get_physical_density(p, cosmo);

  /* Calculate the SFR per gas mass */
  const double SFRpergasmass =
      starform->schmidt_law.mdot_const * sqrt(physical_density);

  /* Store the SFR */
  xp->sf_data.SFR = SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Compute the star-formation rate of a given particle and store
 * it into the #xpart. The star formation is calculated using a pressure
 * law based on Schaye and Dalla Vecchia (2008), the star formation
 * rate is calculated as:
 *
 * \dot{m}_\star = A (1 Msun / pc^-2)^-n m_gas (\gamma/G * f_g *
 *                 pressure)**((n-1)/2)
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR_pressure_law(
    const struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

  /* Get the pressure used for the star formation */
  const double pressure = hydro_get_physical_pressure(p, cosmo);

  /* Calculate the specific star formation rate */
  const double SFRpergasmass =
      starform->pressure_law.SF_normalization *
      pow(pressure, starform->pressure_law.SF_power_law);

  /* Store the SFR */
  xp->sf_data.SFR = SFRpergasmass * hydro_get_mass(p);
}

/**
 * @brief Determine which star formation model is used and call the star
 * formation function either a pressure law or a Schmidt law to store the star
 * formation rate in the xpart. If the star exceeded the maximal allowed density
 * the gas particle is immediately converted to a star
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    const struct part* p, struct xpart* xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {

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
    return;
  }

  /* Determine which star formation model to use */
  switch (starform->SF_law) {

    case colibre_star_formation_schmidt_law:
      star_formation_compute_SFR_schmidt_law(p, xp, starform, phys_const,
                                             hydro_props, cosmo, dt_star);
      break;
    case colibre_star_formation_pressure_law:
      star_formation_compute_SFR_pressure_law(p, xp, starform, phys_const,
                                              hydro_props, cosmo, dt_star);
      break;
    default:
      error("Invalid star formation model!!!");
  }
}

/**
 * @brief Decides whether a new particle should be created or if the hydro
 * particle needs to be transformed.
 *
 * In COLIBRE, we always convert full particles --> return 0
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 *
 * @return 1 if a new spart needs to be created.
 */
INLINE static int star_formation_should_spawn_spart(
    struct part* p, struct xpart* xp, const struct star_formation* starform) {
  return 0;
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

  /* Calculate the probability of forming a star */
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
 * @param convert_part Did we convert a part (or spawned one)?
 */
INLINE static void star_formation_copy_properties(
    const struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cooling_function_data* cooling,
    const int convert_part) {

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
  sp->sf_data.birth_density = hydro_get_physical_density(p, cosmo);

  /* Store the birth temperature in the star particle */
  sp->sf_data.birth_temperature = cooling_get_temperature(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);

  /* Store the birth subgrid density in the star particle */
  sp->sf_data.birth_subgrid_density = p->cooling_data.subgrid_dens;

  /* Store the birth subgrid temperature in the star particle */
  sp->sf_data.birth_subgrid_temperature = p->cooling_data.subgrid_temp;

  /* Store the birth velocity dispersion in the star particle */
  sp->sf_data.birth_velocity_dispersion = p->sf_data.sigma_v2;

  /* Flag that this particle has not done feedback yet */
  sp->SNII_f_E = -1.f;

  /* And also no enrichment */
  sp->last_enrichment_time = sp->birth_time;
  sp->count_since_last_enrichment = -1;

  /* Initialize HII region */
  sp->HIIregion_last_rebuild = -1.f;
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

  /* Read the maximal density that we instantly turn into stars */
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

  /* Get the target number of neighbours and put this as a constant in the
   * virial criterion, scales as target number of neighbours to the power
   * 2/3 */
  const double target_neighbours_to_1_over_3 =
      cbrt(hydro_props->target_neighbours);
  const double target_neighbours_to_2_over_3 =
      target_neighbours_to_1_over_3 * target_neighbours_to_1_over_3;

  /* Store a variable for the virial criterion */
  starform->alpha_virial = parser_get_param_double(
      parameter_file, "COLIBREStarFormation:alpha_virial");

  starform->alpha_virial_inv =
      1. / (starform->alpha_virial * G_newton * target_neighbours_to_2_over_3);

  /* Get the subgrid density threshold */
  starform->subgrid_density_threshold_HpCM3 = parser_get_opt_param_double(
      parameter_file, "COLIBREStarFormation:subgrid_density_threshold_H_p_CM3",
      FLT_MAX);

  /* Calculate the subgrid density threshold using primordial gas mean molecular
   * weights assuming neutral gas */
  starform->subgrid_density_threshold =
      starform->subgrid_density_threshold_HpCM3 *
      phys_const->const_proton_mass * number_density_from_cgs *
      hydro_props->mu_neutral;

  /* Get the subgrid temperature criterion */
  starform->temperature_threshold = parser_get_param_double(
      parameter_file, "COLIBREStarFormation:temperature_threshold_K");

  /* Read in the different star formation models */

  char temp[32] = {0};
  parser_get_param_string(parameter_file, "COLIBREStarFormation:SF_model",
                          temp);

  if (strcmp(temp, "SchmidtLaw") == 0) {
    /* Schmidt model */
    starform->SF_law = colibre_star_formation_schmidt_law;

    /* Get the star formation efficiency */
    starform->schmidt_law.sfe = parser_get_param_double(
        parameter_file, "COLIBREStarFormation:star_formation_efficiency");

    /* Calculate the ff constant */
    const double ff_const = sqrt(3.0 * M_PI / (32.0 * G_newton));

    /* Calculate the constant */
    starform->schmidt_law.mdot_const = starform->schmidt_law.sfe / ff_const;
  } else if (strcmp(temp, "PressureLaw") == 0) {
    /* Pressure model */
    starform->SF_law = colibre_star_formation_pressure_law;

    /* Get the surface density unit Msun / pc^2 in internal units */
    const double Msun_per_pc2 =
        phys_const->const_solar_mass /
        (phys_const->const_parsec * phys_const->const_parsec);

    /* Get the SF surface density unit Msun / kpc^2 / yr in internal units */
    const double kpc = 1000. * phys_const->const_parsec;
    const double Msun_per_kpc2_per_year =
        phys_const->const_solar_mass / (kpc * kpc) / phys_const->const_year;

    /* Read the gas fraction from the file */
    starform->pressure_law.fgas = parser_get_opt_param_double(
        parameter_file, "COLIBREStarFormation:gas_fraction", 1.);

    /* Read the Kennicutt-Schmidt power law exponent */
    starform->pressure_law.KS_power_law = parser_get_param_double(
        parameter_file, "COLIBREStarFormation:KS_exponent");

    /* Calculate the power law of the corresponding star formation Schmidt law
     */
    starform->pressure_law.SF_power_law =
        (starform->pressure_law.KS_power_law - 1.) / 2.;

    /* Read the normalization of the KS law in KS law units */
    starform->pressure_law.KS_normalization_MSUNpYRpKPC2 =
        parser_get_param_double(
            parameter_file,
            "COLIBREStarFormation:KS_normalisation_Msun_p_yr_p_kpc2");

    /* Convert to internal units */
    starform->pressure_law.KS_normalization =
        starform->pressure_law.KS_normalization_MSUNpYRpKPC2 *
        Msun_per_kpc2_per_year;

    /* Calculate the starformation pre-factor (eq. 12 of Schaye & Dalla Vecchia
     * 2008) */
    starform->pressure_law.SF_normalization =
        starform->pressure_law.KS_normalization *
        pow(Msun_per_pc2, -starform->pressure_law.KS_power_law) *
        pow(hydro_gamma * starform->pressure_law.fgas / G_newton,
            starform->pressure_law.SF_power_law);
  }
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 * */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {

  /* Print the star formation properties */
  message("Star formation model is COLIBRE");
  message(
      "The star formation criteria are: minimum over density = %e maximal "
      "density = %e subgrid density criterion = %e alpha_virial = %e "
      "temperature threshold = %e",
      starform->min_over_den, starform->maximal_density_HpCM3,
      starform->subgrid_density_threshold_HpCM3, starform->alpha_virial,
      starform->temperature_threshold);
  switch (starform->SF_law) {
    case colibre_star_formation_schmidt_law:
      message(
          "Star formation law is a Schmidt law: Star formation efficiency = %e",
          starform->schmidt_law.sfe);
      break;
    case colibre_star_formation_pressure_law:
      message(
          "The star formation law is a pressure law (Schaye & Dalla Vecchia "
          "2008): "
          "Kennicutt-Schmidt law normalization = %e Msun/kpc^2/yr, slope of "
          "the Kennicutt-Schmidt law = %e and gas fraction = %e",
          starform->pressure_law.KS_normalization_MSUNpYRpKPC2,
          starform->pressure_law.KS_power_law, starform->pressure_law.fgas);
      break;
    default:
      error("Invalid star formation model!!!");
  }
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
    struct part* p, struct xpart* xp, const struct star_formation* cd,
    const struct cosmology* cosmo) {

  /* Calculate some things before hand */
  const float rho_inv = 1.f / p->rho;
  const float h_inv = 1.f / p->h;
  const float h_inv2 = h_inv * h_inv;

  p->sf_data.sigma_v2 *= rho_inv * h_inv * h_inv2 * cosmo->a2_inv;

  p->sf_data.stars_rho *= h_inv * h_inv2;
  p->sf_data.stars_sigma_v2 *= (1.f/p->sf_data.stars_rho) * 
      h_inv * h_inv2 * cosmo->a2_inv;
}

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
                                      const struct cosmology* cosmo) {

  /* If gas particles do not have neighbours they should not be able to form
   * stars if we set the velocity dispersion to infinity or FLT_MAX this
   * prevents this gas to form stars in any way. */
  p->sf_data.sigma_v2 = FLT_MAX;
}

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
star_formation_first_init_part(const struct phys_const* phys_const,
                               const struct unit_system* us,
                               const struct cosmology* cosmo,
                               const struct star_formation* data,
                               const struct part* p, struct xpart* xp) {}

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
    struct part* p, const struct star_formation* data) {

  /* Initialize the velocity dispersion as 0 before we do the SPH calculation of
   * the velocity dispersion */
  p->sf_data.sigma_v2 = 0.f;
  p->sf_data.scount = 0;
  p->sf_data.stars_rho = 0.f;
  p->sf_data.stars_sigma_v2 = 0.f;
}

/**
 * @brief Split the star formation content of a particle into n pieces
 *
 * We only need to split the SFR if it is positive, i.e. it is not
 * storing the redshift/time of last SF event.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void star_formation_split_part(
    struct part* p, struct xpart* xp, const double n) {

  if (xp->sf_data.SFR > 0.) xp->sf_data.SFR /= n;
}

/**
 * @brief Deal with the case where no spart are available for star formation.
 *
 * @param e The #engine.
 * @param p The #part.
 * @param xp The #xpart.
 */
__attribute__((always_inline)) INLINE static void
star_formation_no_spart_available(const struct engine* e, const struct part* p,
                                  const struct xpart* xp) {
  /* Nothing to do, we just skip it and deal with it next step */
}

/**
 * @brief Compute some information for the star formation model based
 * on all the particles that were read in.
 *
 * This is called once on start-up of the code.
 *
 * Nothing to do here for COLIBRE.
 *
 * @param star_form The #star_formation structure.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_stats(struct star_formation* star_form,
                                const struct engine* e) {}

#endif /* SWIFT_COLIBRE_STAR_FORMATION_H */
