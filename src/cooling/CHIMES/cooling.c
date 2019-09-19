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
 * @file src/cooling/CHIMES/cooling.c
 * @brief CHIMES cooling functions
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "adiabatic_index.h" 
#include "chemistry.h"
#include "cooling.h"
#include "cooling_struct.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/* Tolerances for termination criteria. */
static const float rounding_tolerance = 1.0e-4;

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

  char chimes_data_dir[500]; 
  char string_buffer[500]; 

  /* read the parameters */
  if (sizeof(ChimesFloat) != sizeof(double)) 
    error("CHIMES ERROR: When we parse the CHIMES parameters, we assume that the ChimesFloat type has been set to double. However, sizeof(ChimesFloat) = %lu bytes, while sizeof(double) = %lu bytes.\n", sizeof(ChimesFloat), sizeof(double)); 


  /* Set paths to CHIMES data files. */ 
  parser_get_param_string(parameter_file, "CHIMESCooling:data_path", chimes_data_dir); 
  sprintf(cooling->ChimesGlobalVars.MainDataTablePath, "%s/chimes_main_data.hdf5", chimes_data_dir); 

  parser_get_param_string(parameter_file, "CHIMESCooling:EqmAbundanceTable", string_buffer); 
  sprintf(cooling->ChimesGlobalVars.EqAbundanceTablePath, "%s/EqAbundancesTables/%s", chimes_data_dir, string_buffer); 

  /* Define radiation field and shielding options. */ 
  /* Currently supported options: 
   * UV_field_flag == 0: No UV radiation. 
   *               == 1: Single, constant user-defined spectrum.
   * 
   * Shielding_flag == 0: No shielding. 
   */ 
  cooling->UV_field_flag = parser_get_param_int(parameter_file, "CHIMESCooling:UV_field_flag"); 
  cooling->Shielding_flag = parser_get_param_int(parameter_file, "CHIMESCooling:Shielding_flag"); 

  if (cooling->UV_field_flag == 0) 
    cooling->ChimesGlobalVars.N_spectra = 0; 
  else if (cooling->UV_field_flag == 1) 
    {
      cooling->ChimesGlobalVars.N_spectra = 1; 
      cooling->isotropic_photon_density = (ChimesFloat *) malloc(sizeof(ChimesFloat)); 
      
      parser_get_param_string(parameter_file, "CHIMESCooling:PhotoIonTable", string_buffer); 
      sprintf(cooling->ChimesGlobalVars.PhotoIonTablePath[0], "%s/%s", chimes_data_dir, string_buffer); 

      cooling->isotropic_photon_density[0] = parser_get_param_double(parameter_file, "CHIMESCooling:isotropic_photon_density"); 
    }
  else 
    error("CHIMESCooling: UV_field_flag %d not recognised.", cooling->UV_field_flag); 

  if (cooling->Shielding_flag == 0) 
    cooling->ChimesGlobalVars.cellSelfShieldingOn = 0; 
  else 
    error("CHIMESCooling: Shielding_flag %d not recognised.", cooling->Shielding_flag); 

  /* Cosmic ray ionisation rate of HI. */ 
  cooling->cosmic_ray_rate = parser_get_param_double(parameter_file, "CHIMESCooling:cosmic_ray_rate"); 

  /* CHIMES tolerance parameters */ 
  cooling->ChimesGlobalVars.relativeTolerance = parser_get_param_double(parameter_file, "CHIMESCooling:relativeTolerance"); 
  cooling->ChimesGlobalVars.absoluteTolerance = parser_get_param_double(parameter_file, "CHIMESCooling:absoluteTolerance"); 
  cooling->ChimesGlobalVars.thermalAbsoluteTolerance = parser_get_param_double(parameter_file, "CHIMESCooling:thermalAbsoluteTolerance"); 
  cooling->ChimesGlobalVars.explicitTolerance = parser_get_param_double(parameter_file, "CHIMESCooling:explicitTolerance"); 
  cooling->ChimesGlobalVars.scale_metal_tolerances = parser_get_param_double(parameter_file, "CHIMESCooling:scale_metal_tolerances"); 

  /* Maximum temperature for the molecular network */ 
  cooling->ChimesGlobalVars.T_mol = parser_get_param_double(parameter_file, "CHIMESCooling:T_mol"); 

  /* Determine which metal elements to include 
   * in the CHIMES network. Note that H and He 
   * are always included. */ 
  cooling->ChimesGlobalVars.element_included[0] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeCarbon"); 
  cooling->ChimesGlobalVars.element_included[1] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeNitrogen"); 
  cooling->ChimesGlobalVars.element_included[2] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeOxygen"); 
  cooling->ChimesGlobalVars.element_included[3] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeNeon"); 
  cooling->ChimesGlobalVars.element_included[4] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeMagnesium"); 
  cooling->ChimesGlobalVars.element_included[5] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeSilicon"); 
  cooling->ChimesGlobalVars.element_included[6] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeSulphur"); 
  cooling->ChimesGlobalVars.element_included[7] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeCalcium"); 
  cooling->ChimesGlobalVars.element_included[8] = parser_get_param_int(parameter_file, "CHIMESCooling:IncludeIron"); 
  
  /* For now, the CMB temperature in CHIMES is just set 
   * to the present day value. For cosmological runs, this 
   * will need to be updated for the current redshift. */ 
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  cooling->ChimesGlobalVars.cmb_temperature = cooling->T_CMB_0; 
  
  /* The following CHIMES parameters do not need 
   * to be modified by the user. These are just 
   * hard-coded for now. */ 
  cooling->ChimesGlobalVars.InitIonState = 1; 
  cooling->ChimesGlobalVars.grain_temperature = 10.0; 
  
  /* Physical velocity divergence isn't easily 
   * accessible, so just run with static 
   * molecular cooling for now. */ 
  cooling->ChimesGlobalVars.StaticMolCooling = 1; 

  /* Equilibrium mode: 
   * 0 --> Evolve in non-equilibrium. 
   * 1 --> Evolve with equilibrium abundances. 
   */
  cooling->ChemistryEqmMode = parser_get_param_int(parameter_file, "CHIMESCooling:ChemistryEqmMode"); 

  /* Flag switch thermal evolution on/off: 
   * 0 --> Fixed temperature (but chemical 
   *       abundances can still be evolved). 
   * 1 --> Enable thermal evolution. 
   */ 
  cooling->ThermEvolOn = parser_get_param_int(parameter_file, "CHIMESCooling:ThermEvolOn"); 
    
  /* Optional parameters to define S and Ca 
   * relative to Si. */ 
  cooling->S_over_Si_ratio_in_solar = parser_get_opt_param_float(parameter_file, "CHIMESCooling:S_over_Si_in_solar", 1.f); 
  cooling->Ca_over_Si_ratio_in_solar = parser_get_opt_param_float(parameter_file, "CHIMESCooling:Ca_over_Si_in_solar", 1.f); 

  /* Solar mass fractions of Si, S and Ca are 
   * hard-coded here, because we might not have 
   * access to the COLIBRE cooling tables to 
   * get these. */ 
  cooling->Si_solar_mass_fraction = 6.64948e-4; 
  cooling->S_solar_mass_fraction = 3.09171e-4; 
  cooling->Ca_solar_mass_fraction = 6.41451e-5; 

  /* CHIMES uses a solar metallicity of 0.0129. */ 
  cooling->Zsol = 0.0129; 

  /* Initialise the CHIMES module. */ 
  message("Initialising CHIMES cooling module."); 
  init_chimes(&cooling->ChimesGlobalVars); 
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling #cooling_function_data struct.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  message("Cooling function is 'CHIMES'.");
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * For now, we initialise the CHIMES abundance array to be singly ionised. 
 * This initial state is controlled by the ChimesGlobalVars.InitIonState 
 * parameter, which we have hard-coded to 1 in cooling_init_backend()). 
 * In the future, we will need to consider options to either compute the 
 * initial equilibrium abundances, or read them in from the ICs/snapshot. 
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_first_init_part(const struct phys_const* restrict phys_const,
			     const struct unit_system* restrict us,
			     const struct cosmology* restrict cosmo,
			     const struct cooling_function_data* data, 
			     struct part* restrict p,
			     struct xpart* restrict xp) {
  struct globalVariables ChimesGlobalVars = data->ChimesGlobalVars; 
  struct gasVariables ChimesGasVars; 

  /* Allocate memory to arrays within ChimesGasVars. */
  allocate_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 

  /* Set element abundances from 
   * metal mass fractions. */ 
  chimes_update_element_abundances(data, p, xp, &ChimesGasVars, 0); 

  /* Set initial values for CHIMES 
   * abundance array. */ 
  initialise_gas_abundances(&ChimesGasVars, &ChimesGlobalVars); 

  /* Copy abundances over to xp. */ 
  int i; 
  for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
    xp->cooling_data.chimes_abundances[i] = ChimesGasVars.abundances[i]; 

  /* Zero particle's radiated energy. */ 
  xp->cooling_data.radiated_energy = 0.f; 

  /* Free CHIMES memory. */ 
  free_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 
}

/**
 * @brief Update ChimesGasVars structure. 
 * 
 * Updates the ChimesGasVars structure with the various 
 * thermodynamics quantities of the gas particle that 
 * will be needed for the CHIMES chemistry solver. 
 *  
 * @param u_cgs The internal energy, in cgs units. 
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param hydro_properties the hydro_props struct
 * @param floor_props Properties of the entropy floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param ChimesGasVars CHIMES gasVariables structure. 
 * @param dt The cooling time-step of this particle.
 */
void chimes_update_gas_vars(const double u_cgs,
			    const struct phys_const *phys_const,
			    const struct unit_system *us,
			    const struct cosmology *cosmo,
			    const struct hydro_props *hydro_properties,
			    const struct entropy_floor_properties *floor_props,
			    const struct cooling_function_data *cooling,
			    struct part *restrict p, struct xpart *restrict xp,
			    struct gasVariables *ChimesGasVars, 
			    const float dt_cgs) {
  /* Physical constants that we will 
   * need, in cgs units */ 
  float dimension_k[5] = {1, 2, -2, 0, -1}; 
  double boltzmann_k_cgs = phys_const->const_boltzmann_k * units_general_cgs_conversion_factor(us, dimension_k); 
  double proton_mass_cgs = phys_const->const_proton_mass * units_cgs_conversion_factor(us, UNIT_CONV_MASS); 

  double mu = chimes_mu(cooling, p, xp); 
  ChimesGasVars->temperature = (ChimesFloat) u_cgs * hydro_gamma_minus_one * proton_mass_cgs * mu / boltzmann_k_cgs; 

#if defined(CHEMISTRY_COLIBRE) || defined(CHEMISTRY_EAGLE) 
  float const *metal_fraction = chemistry_get_metal_mass_fraction_for_cooling(p); 
  ChimesFloat XH = (ChimesFloat) metal_fraction[chemistry_element_H]; 
#else 
  /* Without COLIBRE or EAGLE chemistry, 
   * the metal abundances are unavailable. 
   * Set to primordial abundances. */ 
  ChimesFloat XH = 0.75; 
#endif  // CHEMISTRY_COLIBRE || CHEMISTRY_EAGLE 

  ChimesFloat nH = (ChimesFloat) hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass; 
  ChimesGasVars->nH_tot = nH * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY); 

  ChimesGasVars->TempFloor = (ChimesFloat) hydro_properties->minimal_temperature; 
  ChimesGasVars->cr_rate = cooling->cosmic_ray_rate; 
  ChimesGasVars->metallicity = (ChimesFloat) chemistry_get_total_metal_mass_fraction_for_cooling(p) / cooling->Zsol; 
  ChimesGasVars->hydro_timestep = (ChimesFloat) dt_cgs; 
  
  ChimesGasVars->cell_size = 0; 
  ChimesGasVars->constant_heating_rate = 0.0; 
  ChimesGasVars->ForceEqOn = cooling->ChemistryEqmMode; 
  ChimesGasVars->ThermEvolOn = cooling->ThermEvolOn; 
  ChimesGasVars->divVel = 0.0; 

  if (cooling->UV_field_flag == 1) 
    {
      /* Single, constant radiation field. */ 

      /* Copy over normalisation from cooling data. */ 
      ChimesGasVars->isotropic_photon_density[0] = cooling->isotropic_photon_density[0]; 
      
      /* Copy over spectrum parameters from 
       * chimes tables. */ 
      ChimesGasVars->G0_parameter[0] = chimes_table_spectra.G0_parameter[0]; 
      ChimesGasVars->H2_dissocJ[0] = chimes_table_spectra.H2_dissocJ[0]; 
    }

  /* Doppler broadening parameter, for 
   * H2 self-shielding, is hard-coded 
   * to 7.1 km/s for now. This is a 
   * typical value for GMCs in the 
   * Milky Way. */ 
  ChimesGasVars->doppler_broad = 7.1; 
}

/** 
 * @brief Update CHIMES element abundances. 
 * 
 * Updates the element abundances based on the 
 * particle's current mass fractions. If the 
 * element abundances have changed since the 
 * last time this routine was called, for example 
 * due to enrichment or metal diffusion, then the 
 * ion and molecule abundances will now be 
 * inconsistent with their corresponding total 
 * element abundances. Therefore, if mode == 1, 
 * we call the check_constraint_equations() 
 * routine from the CHIMES module, which adds all 
 * species of each element and compares to the 
 * total abundance of that element. If they are 
 * inconsistent, the species abundances are 
 * adjusted, preserving their relative fractions. 
 * If a metal is now present that was not present 
 * before (e.g. if a particle was primordial but 
 * has recently been enriched), then that metal is 
 * introduced as fully neutral. 
 * We do not call check_constraint_equations() if 
 * mode == 0 (for example, when we initialise the 
 * abundance arrays for the first time). 
 * 
 * @param p Pointer to the particle data.
 * @param ChimesGasVars CHIMES gasVariables structure. 
 */ 
void chimes_update_element_abundances(const struct cooling_function_data *cooling,
				      struct part *restrict p, 
				      struct xpart *restrict xp,
				      struct gasVariables *ChimesGasVars, 
				      const int mode) {
  struct globalVariables ChimesGlobalVars = cooling->ChimesGlobalVars; 
  int i; 

  /* Get element mass fractions */ 
#if defined(CHEMISTRY_COLIBRE) || defined(CHEMISTRY_EAGLE) 
  float const *metal_fraction = chemistry_get_metal_mass_fraction_for_cooling(p); 
  float XH = metal_fraction[chemistry_element_H]; 
  ChimesGasVars->element_abundances[0] = (ChimesFloat) metal_fraction[chemistry_element_He] / (4.0 * XH); 
  ChimesGasVars->element_abundances[1] = (ChimesFloat) metal_fraction[chemistry_element_C] / (12.0 * XH); 
  ChimesGasVars->element_abundances[2] = (ChimesFloat) metal_fraction[chemistry_element_N] / (14.0 * XH); 
  ChimesGasVars->element_abundances[3] = (ChimesFloat) metal_fraction[chemistry_element_O] / (16.0 * XH); 
  ChimesGasVars->element_abundances[4] = (ChimesFloat) metal_fraction[chemistry_element_Ne] / (20.0 * XH); 
  ChimesGasVars->element_abundances[5] = (ChimesFloat) metal_fraction[chemistry_element_Mg] / (24.0 * XH); 
  ChimesGasVars->element_abundances[6] = (ChimesFloat) metal_fraction[chemistry_element_Si] / (28.0 * XH); 
  ChimesGasVars->element_abundances[9] = (ChimesFloat) metal_fraction[chemistry_element_Fe] / (56.0 * XH); 
  
  ChimesGasVars->element_abundances[7] = (ChimesFloat) metal_fraction[chemistry_element_Si] * cooling->S_over_Si_ratio_in_solar * (cooling->S_solar_mass_fraction / cooling->Si_solar_mass_fraction) / (32.0 * XH); 
  ChimesGasVars->element_abundances[8] = (ChimesFloat) metal_fraction[chemistry_element_Si] * cooling->Ca_over_Si_ratio_in_solar * (cooling->Ca_solar_mass_fraction / cooling->Si_solar_mass_fraction) / (40.0 * XH); 

  /* Zero the abundances of any elements 
   * that are not included in the network. */ 
  for (i = 0; i < 9; i++) 
    ChimesGasVars->element_abundances[i + 1] *= (ChimesFloat) ChimesGlobalVars.element_included[i]; 
#else 
  /* Without COLIBRE or EAGLE chemistry, 
   * the metal abundances are unavailable. 
   * Set to primordial abundances. */ 
  ChimesGasVars->element_abundances[0] = 0.0833;  // He 

  for (i = 1; i < 10; i++) 
    ChimesGasVars->element_abundances[i] = 0.0; 
#endif  // CHEMISTRY_COLIBRE || CHEMISTRY_EAGLE 

  if (mode == 1) 
    check_constraint_equations(ChimesGasVars, &ChimesGlobalVars); 

} 

/**
 * @brief Apply the cooling function to a particle.
 *
 * We use the CHIMES module to integrate the cooling rate 
 * and chemical abundances over the time-step. 
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

  /* CHIMES structures. */
  struct globalVariables ChimesGlobalVars = cooling->ChimesGlobalVars; 
  struct gasVariables ChimesGasVars; 

  /* Allocate memory to arrays within ChimesGasVars. */ 
  allocate_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 

  /* Copy abundances over from xp to ChimesGasVars. */
  int i; 
  for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
    ChimesGasVars.abundances[i] = xp->cooling_data.chimes_abundances[i]; 


  /* Update element abundances from metal mass 
   * fractions. We need to do this here, and not 
   * later on in chimes_update_gas_vars(), because 
   * the element abundances need to be set in 
   * ChimesGasVars before we can calculate the 
   * mean molecular weight. */ 
  chimes_update_element_abundances(cooling, p, xp, &ChimesGasVars, 1); 

  /* Get internal energy at the last kick step */
  const float u_start = hydro_get_physical_internal_energy(p, xp, cosmo);

  /* Get the change in internal energy due to hydro forces */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  /* Get internal energy at the end of the next kick step (assuming dt does not
   * increase) */
  double u_0 = (u_start + hydro_du_dt * dt_therm);

  /* Check for minimal energy. Note that, to 
   * maintain consistency with the temperature 
   * floor that is imposed within CHIMES, we compute 
   * the minimal energy from the minimal temperature 
   * based on the particle's actual mean molecular 
   * weight mu, rather than an assumed constant mu. 
 */
  double mu = chimes_mu(cooling, p, xp); 
  double minimal_internal_energy; 
  minimal_internal_energy = hydro_properties->minimal_temperature; 
  minimal_internal_energy *= hydro_one_over_gamma_minus_one; 
  minimal_internal_energy *= (phys_const->const_boltzmann_k / phys_const->const_proton_mass); 
  minimal_internal_energy /= mu; 

  u_0 = max(u_0, minimal_internal_energy);

  /* Convert to CGS units */
  const double u_start_cgs = u_start * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  const double u_0_cgs = u_0 * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Update the ChimesGasVars structure with the 
   * particle's thermodynamic variables. */ 
  chimes_update_gas_vars(u_0_cgs, phys_const, us, cosmo, hydro_properties, floor_props, cooling, p, xp, &ChimesGasVars, dt_cgs); 

  /* Call CHIMES to integrate the chemistry 
   * and cooling over the time-step. */ 
  chimes_network(&ChimesGasVars, &ChimesGlobalVars); 

  /* Physical constants that we will 
   * need, in cgs units */ 
  float dimension_k[5] = {1, 2, -2, 0, -1}; 
  double boltzmann_k_cgs = phys_const->const_boltzmann_k * units_general_cgs_conversion_factor(us, dimension_k); 
  double proton_mass_cgs = phys_const->const_proton_mass * units_cgs_conversion_factor(us, UNIT_CONV_MASS); 

  /* Compute the internal energy at the end of the 
   * step using the final temperature from CHIMES. */
  double u_final_cgs;

  mu = chimes_mu(cooling, p, xp); 

  u_final_cgs = (double) ChimesGasVars.temperature; 
  u_final_cgs *= hydro_one_over_gamma_minus_one; 
  u_final_cgs *= boltzmann_k_cgs / proton_mass_cgs; 
  u_final_cgs /= mu; 

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u_cgs = u_final_cgs - u_start_cgs;

  /* Convert back to internal units */
  double delta_u = delta_u_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS); 

  /* We now need to check that we are not going to go below any of the limits */

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_floor = gas_internal_energy_from_entropy(rho, A_floor);

  /* Recompute new minimal internal energy 
  * from the new mean molecular weight. */
  minimal_internal_energy = hydro_properties->minimal_temperature; 
  minimal_internal_energy *= hydro_one_over_gamma_minus_one; 
  minimal_internal_energy *= (phys_const->const_boltzmann_k / phys_const->const_proton_mass); 
  minimal_internal_energy /= mu; 

  /* Largest of both limits */
  const double u_limit = max(minimal_internal_energy, u_floor);

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

  /* Copy abundances from ChimesGasVars back to xp. */ 
  for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
    xp->cooling_data.chimes_abundances[i] = ChimesGasVars.abundances[i]; 

  /* Free CHIMES memory. */ 
  free_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 

}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
float cooling_get_radiated_energy(const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy; 
}

/**
 *
 * @brief Calculate mean molecular weight from CHIMES abundances. 
 * 
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */ 
double chimes_mu(const struct cooling_function_data *cooling,
		 const struct part *restrict p, 
		 const struct xpart *restrict xp) {
  struct globalVariables ChimesGlobalVars = cooling->ChimesGlobalVars; 
  double numerator = 1.0; 

    /* Get element mass fractions */ 
#if defined(CHEMISTRY_COLIBRE) || defined(CHEMISTRY_EAGLE) 
  float const *metal_fraction = chemistry_get_metal_mass_fraction_for_cooling(p); 
  float XH = metal_fraction[chemistry_element_H]; 

  numerator += (double) metal_fraction[chemistry_element_He] / XH; 

  if (ChimesGlobalVars.element_included[0] == 1) 
    numerator += (double) metal_fraction[chemistry_element_C] / XH;

  if (ChimesGlobalVars.element_included[1] == 1) 
    numerator += (double) metal_fraction[chemistry_element_N] / XH; 

  if (ChimesGlobalVars.element_included[2] == 1) 
    numerator += (double) metal_fraction[chemistry_element_O] / XH; 

  if (ChimesGlobalVars.element_included[3] == 1) 
    numerator += (double) metal_fraction[chemistry_element_Ne] / XH; 

  if (ChimesGlobalVars.element_included[4] == 1) 
    numerator += (double) metal_fraction[chemistry_element_Mg] / XH; 

  if (ChimesGlobalVars.element_included[5] == 1) 
    numerator += (double) metal_fraction[chemistry_element_Si] / XH; 

  if (ChimesGlobalVars.element_included[8] == 1) 
    numerator += (double) metal_fraction[chemistry_element_Fe] / XH; 

  if (ChimesGlobalVars.element_included[6] == 1) 
    numerator += (double) metal_fraction[chemistry_element_Si] * cooling->S_over_Si_ratio_in_solar * (cooling->S_solar_mass_fraction / cooling->Si_solar_mass_fraction) / XH;

  if (ChimesGlobalVars.element_included[7] == 1) 
    numerator += (double) metal_fraction[chemistry_element_Si] * cooling->Ca_over_Si_ratio_in_solar * (cooling->Ca_solar_mass_fraction / cooling->Si_solar_mass_fraction) / XH;  
#else 
  /* Without COLIBRE or EAGLE chemistry, 
   * the metal abundances are unavailable. 
   * Set to primordial abundances. */ 
  numerator += 0.0833 * 4.0; 
#endif  // CHEMISTRY_COLIBRE || CHEMISTRY_EAGLE 

  int i; 
  double denominator = 0.0; 
  
  for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
    denominator += (double) xp->cooling_data.chimes_abundances[i]; 
  
  return numerator / denominator; 
} 

/**
 * @brief Compute the temperature of a #part based on the cooling function.
 * 
 * This uses the mean molecular weight computed from the 
 * CHIMES abundance array. 
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
float cooling_get_temperature(const struct phys_const* restrict phys_const,
			      const struct hydro_props* restrict hydro_props,
			      const struct unit_system* restrict us,
			      const struct cosmology* restrict cosmo,
			      const struct cooling_function_data* restrict cooling,
			      const struct part* restrict p, 
			      const struct xpart* restrict xp) {
  /* Physical constants, in cgs units */ 
  float dimension_k[5] = {1, 2, -2, 0, -1}; 
  double boltzmann_k_cgs = phys_const->const_boltzmann_k * units_general_cgs_conversion_factor(us, dimension_k); 
  double proton_mass_cgs = phys_const->const_proton_mass * units_cgs_conversion_factor(us, UNIT_CONV_MASS); 

  /* Mean molecular weight */
  const double mu = chimes_mu(cooling, p, xp); 

  /* Internal energy, in code units. */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo); 

  /* Convert to cgs. */ 
  const double u_cgs = u * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS); 

  /* Return particle temperature */
  return (float) hydro_gamma_minus_one * u_cgs * mu * (proton_mass_cgs / boltzmann_k_cgs); 
}
