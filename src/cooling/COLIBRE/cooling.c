/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
/**
 * @file src/cooling/COLIBRE/cooling.c
 * @brief COLIBRE cooling functions
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling.h"
#include "cooling_rates.h"
#include "cooling_struct.h"
#include "cooling_tables.h"
#include "entropy_floor.h"
#include "error.h"
#include "exp10.h"
#include "hydro.h"
#include "interpolate.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/* Maximum number of iterations for
 * bisection integration schemes */
static const int bisection_max_iterations = 150;

/* Tolerances for termination criteria. */
static const float explicit_tolerance = 0.05;
static const float bisection_tolerance = 1.0e-6;
static const float rounding_tolerance = 1.0e-4;
static const double bracket_factor = 1.5; /* sqrt(1.1) */

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * Also calls the additional H reionisation energy injection if need be.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The space data, including a pointer to array of particles
 */
void cooling_update(const struct cosmology *cosmo,
                    struct cooling_function_data *cooling, struct space *s) {

  /* Extra energy for reionization? */
  if (!cooling->H_reion_done) {

    /* Does this timestep straddle Hydrogen reionization? If so, we need to
     * input extra heat */
    if (cosmo->z <= cooling->H_reion_z && cosmo->z_old > cooling->H_reion_z) {

      if (s == NULL) error("Trying to do H reionization on an empty space!");

      /* Inject energy to all particles */
      cooling_Hydrogen_reionization(cooling, cosmo, s);

      /* Flag that reionization happened */
      cooling->H_reion_done = 1;
    }
  }
}

/**
 * @brief Bisection integration scheme
 *
 * @param u_ini_cgs Internal energy at beginning of hydro step in CGS.
 * @param n_H_cgs Hydrogen number density in CGS.
 * @param redshift Current redshift.
 * @param n_H_index Particle hydrogen number density index.
 * @param d_n_H Particle hydrogen number density offset.
 * @param met_index Particle metallicity index.
 * @param d_met Particle metallicity offset.
 * @param red_index Redshift index.
 * @param d_red Redshift offset.
 * @param Lambda_He_reion_cgs Cooling rate coming from He reionization.
 * @param ratefact_cgs Multiplication factor to get a cooling rate.
 * @param cooling #cooling_function_data structure.
 * @param abundance_ratio Array of ratios of metal abundance to solar.
 * @param dt_cgs timestep in CGS.
 * @param ID ID of the particle (for debugging).
 */
static INLINE double bisection_iter(
    const double u_ini_cgs, const double n_H_cgs, const double redshift,
    int n_H_index, float d_n_H, int met_index, float d_met, int red_index,
    float d_red, double Lambda_He_reion_cgs, double ratefact_cgs,
    const struct cooling_function_data *restrict cooling,
    const float abundance_ratio[chemistry_element_count + 3], double dt_cgs,
    long long ID) {

  /* Bracketing */
  double u_lower_cgs = u_ini_cgs;
  double u_upper_cgs = u_ini_cgs;

  /*************************************/
  /* Let's get a first guess           */
  /*************************************/

  double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      colibre_cooling_rate(log10(u_ini_cgs), redshift, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, met_index, d_met, red_index, d_red,
                           cooling, 0, 0, 0, 0);

  /*************************************/
  /* Let's try to bracket the solution */
  /*************************************/

  if (LambdaNet_cgs < 0) {

    /* we're cooling! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, 0, 0, 0, 0);

    int i = 0;
    while (u_lower_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs >
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs /= bracket_factor;
      u_upper_cgs /= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          colibre_cooling_rate(log10(u_lower_cgs), redshift, n_H_cgs,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, 0, 0, 0, 0);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "cooling",
          ID);
    }
  } else {

    /* we are heating! */
    u_lower_cgs /= bracket_factor;
    u_upper_cgs *= bracket_factor;

    /* Compute a new rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, 0, 0, 0, 0);

    int i = 0;
    while (u_upper_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs <
               0 &&
           i < bisection_max_iterations) {

      u_lower_cgs *= bracket_factor;
      u_upper_cgs *= bracket_factor;

      /* Compute a new rate */
      LambdaNet_cgs =
          Lambda_He_reion_cgs +
          colibre_cooling_rate(log10(u_upper_cgs), redshift, n_H_cgs,
                               abundance_ratio, n_H_index, d_n_H, met_index,
                               d_met, red_index, d_red, cooling, 0, 0, 0, 0);
      i++;
    }

    if (i >= bisection_max_iterations) {
      error(
          "particle %llu exceeded max iterations searching for bounds when "
          "heating",
          ID);
    }
  }

  /********************************************/
  /* We now have an upper and lower bound.    */
  /* Let's iterate by reducing the bracketing */
  /********************************************/

  /* bisection iteration */
  int i = 0;
  double u_next_cgs;

  do {

    /* New guess */
    u_next_cgs = 0.5 * (u_lower_cgs + u_upper_cgs);

    /* New rate */
    LambdaNet_cgs =
        Lambda_He_reion_cgs +
        colibre_cooling_rate(log10(u_next_cgs), redshift, n_H_cgs,
                             abundance_ratio, n_H_index, d_n_H, met_index,
                             d_met, red_index, d_red, cooling, 0, 0, 0, 0);

    /* Where do we go next? */
    if (u_next_cgs - u_ini_cgs - LambdaNet_cgs * ratefact_cgs * dt_cgs > 0.0) {
      u_upper_cgs = u_next_cgs;
    } else {
      u_lower_cgs = u_next_cgs;
    }

    i++;
  } while (fabs(u_upper_cgs - u_lower_cgs) / u_next_cgs > bisection_tolerance &&
           i < bisection_max_iterations);

  if (i >= bisection_max_iterations)
    error("Particle id %llu failed to converge", ID);

  return u_upper_cgs;
}

/**
 * @brief Apply the cooling function to a particle.
 *
 * We want to compute u_new such that u_new = u_old + dt * du/dt(u_new, X),
 * where X stands for the metallicity, density and redshift. These are
 * kept constant.
 *
 * We first compute du/dt(u_old). If dt * du/dt(u_old) is small enough, we
 * use an explicit integration and use this as our solution.
 *
 * Otherwise, we try to find a solution to the implicit time-integration
 * problem. This leads to the root-finding problem:
 *
 * f(u_new) = u_new - u_old - dt * du/dt(u_new, X) = 0
 *
 * A bisection scheme is used.
 * This is done by first bracketing the solution and then iterating
 * towards the solution by reducing the window down to a certain tolerance.
 * Note there is always at least one solution since
 * f(+inf) is < 0 and f(-inf) is > 0.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_properties the hydro_props struct
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The cooling time-step of this particle.
 * @param dt_therm The hydro time-step of this particle.
 */
void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm) {

  /* No cooling happens over zero time */
  if (dt == 0.) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--cooling runtime flag?");
#endif

  /* Get internal energy at the last kick step */
  const float u_start = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Get the change in internal energy due to hydro forces */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Get internal energy at the end of the next kick step (assuming dt does not
   * increase) */
  double u_0 = (u_start + hydro_du_dt * dt_therm);

  /* Check for minimal energy */
  u_0 = max(u_0, hydro_properties->minimal_internal_energy);

  /* Convert to CGS units */
  const double u_start_cgs = u_start * cooling->internal_energy_to_cgs;
  const double u_0_cgs = u_0 * cooling->internal_energy_to_cgs;
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Change in redshift over the course of this time-step
     (See cosmology theory document for the derivation) */
  const double delta_redshift = -dt * cosmo->H * cosmo->a_inv;

  /* Get this particle's abundance ratios compared to solar
   * Note that we need to add S and Ca that are in the tables but not tracked
   * by the particles themselves.
   * The order is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe] */
  float abundance_ratio[chemistry_element_count + 3];
  float logZZsol = abundance_ratio_to_solar(p, cooling, abundance_ratio);

  /* Get the Hydrogen and Helium mass fractions */
  const float XH =
      p->chemistry_data.smoothed_metal_mass_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction into Hydrogen number density */
  const double n_H =
      hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  /* ratefact = n_H * n_H / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  const double ratefact_cgs = n_H_cgs * (XH * cooling->inv_proton_mass_cgs);

  /* compute hydrogen number density, metallicity and redshift indices and
   * offsets (These are fixed for any value of u, so no need to recompute them)
   */

  float d_red, d_met, d_n_H;
  int red_index, met_index, n_H_index;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
               &red_index, &d_red);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Start by computing the cooling (heating actually) rate from Helium
     re-ionization as this needs to be added on no matter what */

  /* Get helium and hydrogen reheating term */
  const double Helium_reion_heat_cgs =
      eagle_helium_reionization_extraheat(cosmo->z, delta_redshift, cooling);

  /* Convert this into a rate */
  const double Lambda_He_reion_cgs =
      Helium_reion_heat_cgs / (dt_cgs * ratefact_cgs);

  /* Let's compute the internal energy at the end of the step */
  double u_final_cgs;

  /* First try an explicit integration (note we ignore the derivative) */
  const double LambdaNet_cgs =
      Lambda_He_reion_cgs +
      colibre_cooling_rate(log10(u_0_cgs), cosmo->z, n_H_cgs, abundance_ratio,
                           n_H_index, d_n_H, met_index, d_met, red_index, d_red,
                           cooling, 0, 0, 0, 0);

  /* if cooling rate is small, take the explicit solution */
  if (fabs(ratefact_cgs * LambdaNet_cgs * dt_cgs) <
      explicit_tolerance * u_0_cgs) {

    u_final_cgs = u_0_cgs + ratefact_cgs * LambdaNet_cgs * dt_cgs;

  } else {

    u_final_cgs =
        bisection_iter(u_0_cgs, n_H_cgs, cosmo->z, n_H_index, d_n_H, met_index,
                       d_met, red_index, d_red, Lambda_He_reion_cgs,
                       ratefact_cgs, cooling, abundance_ratio, dt_cgs, p->id);
  }

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u_cgs = u_final_cgs - u_start_cgs;

  /* Convert back to internal units */
  double delta_u = delta_u_cgs * cooling->internal_energy_from_cgs;

  /* We now need to check that we are not going to go below any of the limits */

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_floor = gas_internal_energy_from_entropy(rho, A_floor);

  /* Absolute minimum */
  const double u_minimal = hydro_properties->minimal_internal_energy;

  /* Largest of both limits */
  const double u_limit = max(u_minimal, u_floor);

  /* First, check whether we may end up below the minimal energy after
   * this step 1/2 kick + another 1/2 kick that could potentially be for
   * a time-step twice as big. We hence check for 1.5 delta_u. */
  if (u_start + 1.5 * delta_u < u_limit) {
    delta_u = (u_limit - u_start) / 1.5;
  }

  /* Second, check whether the energy used in the prediction could get negative.
   * We need to check for the 1/2 dt kick followed by a full time-step drift
   * that could potentially be for a time-step twice as big. We hence check
   * for 2.5 delta_u but this time against 0 energy not the minimum.
   * To avoid numerical rounding bringing us below 0., we add a tiny tolerance.
   */
  if (u_start + 2.5 * delta_u < 0.) {
    delta_u = -u_start / (2.5 + rounding_tolerance);
  }

  /* Turn this into a rate of change (including cosmology term) */
  const float cooling_du_dt = delta_u / dt_therm;

  /* Update the internal energy time derivative */
  hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cooling_du_dt * dt;
}

/**
 * @brief Computes the cooling time-step.
 *
 * The time-step is not set by the properties of cooling.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const #phys_const data struct.
 * @param us The internal system of units.
 * @param cosmo #cosmology struct.
 * @param hydro_props the properties of the hydro scheme.
 * @param p #part data.
 * @param xp extended particle data.
 */
__attribute__((always_inline)) INLINE float cooling_timestep(
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props, const struct part *restrict p,
    const struct xpart *restrict xp) {

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const #phys_const data structure.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
__attribute__((always_inline)) INLINE void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp) {

  xp->cooling_data.radiated_energy = 0.f;
}

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 *
 * The temperature returned is consistent with the cooling rates.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (cooling->Redshifts == NULL)
    error(
        "Cooling function has not been initialised. Did you forget the "
        "--temperature runtime flag?");
#endif

  /* Get physical internal energy */
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const double u_cgs = u * cooling->internal_energy_to_cgs;

  /* Get the Hydrogen mass fraction */
  const float XH =
      p->chemistry_data.smoothed_metal_mass_fraction[chemistry_element_H];

  /* Convert Hydrogen mass fraction into Hydrogen number density */
  const float rho = hydro_get_physical_density(p, cosmo);
  const double n_H = rho * XH / phys_const->const_proton_mass;
  const double n_H_cgs = n_H * cooling->number_density_to_cgs;

  const float logZZsol = log10(
      p->chemistry_data.smoothed_metal_mass_fraction_total / cooling->Zsol[0]);

  /* compute hydrogen number density, metallicity and redshift indices and
   * offsets  */

  float d_red, d_met, d_n_H;
  int red_index, met_index, n_H_index;

  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, cosmo->z,
               &red_index, &d_red);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs),
               &n_H_index, &d_n_H);

  /* Compute the log10 of the temperature by interpolating the table */
  const double log_10_T =
      colibre_convert_u_to_temp(log10(u_cgs), cosmo->z, n_H_index, d_n_H,
                                met_index, d_met, red_index, d_red, cooling);

  /* Undo the log! */
  return exp10(log_10_T);
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp #xpart data struct
 */
__attribute__((always_inline)) INLINE float cooling_get_radiated_energy(
    const struct xpart *restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Inject a fixed amount of energy to each particle in the simulation
 * to mimic Hydrogen reionization.
 *
 * @param cooling The properties of the cooling model.
 * @param cosmo The cosmological model.
 * @param s The #space containing the particles.
 */
void cooling_Hydrogen_reionization(const struct cooling_function_data *cooling,
                                   const struct cosmology *cosmo,
                                   struct space *s) {

  struct part *parts = s->parts;
  struct xpart *xparts = s->xparts;

  /* Energy to inject in internal units */
  const float extra_heat =
      cooling->H_reion_heat_cgs * cooling->internal_energy_from_cgs;

  message("Applying extra energy for H reionization!");

  /* Loop through particles and set new heat */
  for (size_t i = 0; i < s->nr_parts; i++) {

    struct part *p = &parts[i];
    struct xpart *xp = &xparts[i];

    const float old_u = hydro_get_physical_internal_energy(p, xp, cosmo);
    const float new_u = old_u + extra_heat;

    hydro_set_physical_internal_energy(p, xp, cosmo, new_u);
    hydro_set_drifted_physical_internal_energy(p, cosmo, new_u);
  }
}

/**
 * @brief Initialises properties stored in the cooling_function_data struct
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const #phys_const data structure
 * @param cooling #cooling_function_data struct to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          struct cooling_function_data *cooling) {

  /* read some parameters */

  /* Despite the names, the values of H_reion_heat_cgs and He_reion_heat_cgs
   * that are read in are actually in units of electron volts per proton mass.
   * We later convert to units just below */

  parser_get_param_string(parameter_file, "COLIBRECooling:dir_name",
                          cooling->cooling_table_path);

  cooling->H_reion_done = 0;
  cooling->H_reion_z =
      parser_get_param_float(parameter_file, "COLIBRECooling:H_reion_z");
  cooling->H_reion_heat_cgs =
      parser_get_param_float(parameter_file, "COLIBRECooling:H_reion_eV_p_H");
  cooling->He_reion_z_centre = parser_get_param_float(
      parameter_file, "COLIBRECooling:He_reion_z_centre");
  cooling->He_reion_z_sigma =
      parser_get_param_float(parameter_file, "COLIBRECooling:He_reion_z_sigma");
  cooling->He_reion_heat_cgs =
      parser_get_param_float(parameter_file, "COLIBRECooling:He_reion_eV_p_H");

  /* Optional parameters to correct the abundances */
  cooling->Ca_over_Si_ratio_in_solar = parser_get_opt_param_float(
      parameter_file, "COLIBRECooling:Ca_over_Si_in_solar", 1.f);
  cooling->S_over_Si_ratio_in_solar = parser_get_opt_param_float(
      parameter_file, "COLIBRECooling:S_over_Si_in_solar", 1.f);

  /* Convert H_reion_heat_cgs and He_reion_heat_cgs to cgs
   * (units used internally by the cooling routines). This is done by
   * multiplying by 'eV/m_H' in internal units, then converting to cgs units.
   * Note that the dimensions of these quantities are energy/mass = velocity^2
   */

  cooling->H_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  cooling->He_reion_heat_cgs *=
      phys_const->const_electron_volt / phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  read_cooling_header(cooling);
  read_cooling_tables(cooling);

  /* Compute conversion factors */
  cooling->internal_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  cooling->internal_energy_from_cgs = 1. / cooling->internal_energy_to_cgs;
  cooling->number_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Store some constants in CGS units */
  const double proton_mass_cgs =
      phys_const->const_proton_mass *
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  cooling->inv_proton_mass_cgs = 1. / proton_mass_cgs;
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Compute the coefficient at the front of the Compton cooling expression */
  const double radiation_constant =
      4. * phys_const->const_stefan_boltzmann / phys_const->const_speed_light_c;
  const double compton_coefficient =
      4. * radiation_constant * phys_const->const_thomson_cross_section *
      phys_const->const_boltzmann_k /
      (phys_const->const_electron_mass * phys_const->const_speed_light_c);
  const float dimension_coefficient[5] = {1, 2, -3, 0, -5};

  /* This should be ~1.0178085e-37 g cm^2 s^-3 K^-5 */
  const double compton_coefficient_cgs =
      compton_coefficient *
      units_general_cgs_conversion_factor(us, dimension_coefficient);

#ifdef SWIFT_DEBUG_CHECKS
  const double expected_compton_coefficient_cgs = 1.0178085e-37;
  if (fabs(compton_coefficient_cgs - expected_compton_coefficient_cgs) /
          expected_compton_coefficient_cgs >
      0.01)
    error("compton coefficient incorrect.");
#endif

  /* And now the Compton rate */
  cooling->compton_rate_cgs = compton_coefficient_cgs * cooling->T_CMB_0 *
                              cooling->T_CMB_0 * cooling->T_CMB_0 *
                              cooling->T_CMB_0;
}

/**
 * @brief Restore cooling tables (if applicable) after
 * restart
 *
 * @param cooling the #cooling_function_data structure
 * @param cosmo #cosmology structure
 */
void cooling_restore_tables(struct cooling_function_data *cooling,
                            const struct cosmology *cosmo) {

  read_cooling_header(cooling);
  read_cooling_tables(cooling);

  cooling_update(cosmo, cooling, /*space=*/NULL);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling #cooling_function_data struct.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  message("Cooling function is 'COLIBRE'.");
}

/**
 * @brief Clean-up the memory allocated for the cooling routines
 *
 * We simply free all the arrays.
 *
 * @param cooling the cooling data structure.
 */
void cooling_clean(struct cooling_function_data *cooling) {

  /* Free the side arrays */
  free(cooling->Redshifts);
  free(cooling->nH);
  free(cooling->Temp);
  free(cooling->Metallicity);
  free(cooling->Therm);
  free(cooling->LogAbundances);
  free(cooling->Abundances);
  free(cooling->Abundances_inv);
  free(cooling->atomicmass);
  free(cooling->LogMassFractions);
  free(cooling->MassFractions);

  /* Free the tables */
  free(cooling->table.Tcooling);
  free(cooling->table.Ucooling);
  free(cooling->table.Theating);
  free(cooling->table.Uheating);
  free(cooling->table.Telectron_fraction);
  free(cooling->table.Uelectron_fraction);
  free(cooling->table.T_from_U);
  free(cooling->table.U_from_T);
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the cooling routines. Helps debugging. */
  struct cooling_function_data cooling_copy = *cooling;
  cooling_copy.Redshifts = NULL;
  cooling_copy.nH = NULL;
  cooling_copy.Temp = NULL;
  cooling_copy.Metallicity = NULL;
  cooling_copy.Therm = NULL;
  cooling_copy.LogAbundances = NULL;
  cooling_copy.Abundances = NULL;
  cooling_copy.Abundances_inv = NULL;
  cooling_copy.atomicmass = NULL;
  cooling_copy.LogMassFractions = NULL;
  cooling_copy.MassFractions = NULL;

  cooling_copy.table.Tcooling = NULL;
  cooling_copy.table.Theating = NULL;
  cooling_copy.table.Telectron_fraction = NULL;
  cooling_copy.table.Ucooling = NULL;
  cooling_copy.table.Uheating = NULL;
  cooling_copy.table.Uelectron_fraction = NULL;
  cooling_copy.table.T_from_U = NULL;
  cooling_copy.table.U_from_T = NULL;

  restart_write_blocks((void *)&cooling_copy,
                       sizeof(struct cooling_function_data), 1, stream,
                       "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the cooling tables by
 * re-reading them.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo) {
  restart_read_blocks((void *)cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  cooling_restore_tables(cooling, cosmo);
}
