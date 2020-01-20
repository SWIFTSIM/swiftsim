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
#include "colibre_tables_restrict.h" 
#include "colibre_subgrid.h" 
#include "entropy_floor.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief Initialises properties stored in the cooling_function_data struct
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const #phys_const data structure
 * @param hydro_props The properties of the hydro scheme. 
 * @param cooling #cooling_function_data struct to initialize
 */
void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling) {

  char chimes_data_dir[256]; 
  char string_buffer[196]; 

  /* read the parameters */

  /* Set paths to CHIMES data files. */ 
  parser_get_param_string(parameter_file, "CHIMESCooling:data_path", chimes_data_dir); 
  sprintf(cooling->ChimesGlobalVars.MainDataTablePath, "%s/chimes_main_data.hdf5", chimes_data_dir); 

  parser_get_param_string(parameter_file, "CHIMESCooling:EqmAbundanceTable", string_buffer); 
  sprintf(cooling->ChimesGlobalVars.EqAbundanceTablePath, "%s/EqAbundancesTables/%s", chimes_data_dir, string_buffer); 

  /* Define radiation field and shielding options. */ 
  /* Currently supported options: 
   * UV_field_flag == 0: No UV radiation. 
   *               == 1: Single user-defined spectrum. 
   *               == 2: COLIBRE UVB+ISRF (and scale cr_rate and f_dust).
   * 
   * use_redshift_dependent_UVB == 0: Constant UVB. 
   *                            == 1: Redshift-dependent cross sections. 
   *                            == 2: Redshift-dependent cross sections and eqm tables. 
   * 
   * Shielding_flag == 0: No shielding.
   *                == 1: Jeans shielding length. 
   *                == 2: COLIBRE shielding length. 
   */ 
  cooling->UV_field_flag = parser_get_param_int(parameter_file, "CHIMESCooling:UV_field_flag"); 
  cooling->Shielding_flag = parser_get_param_int(parameter_file, "CHIMESCooling:Shielding_flag"); 
  cooling->use_redshift_dependent_UVB = parser_get_param_int(parameter_file, "CHIMESCooling:use_redshift_dependent_UVB"); 

  cooling->radiation_field_normalisation_factor = 0.0; 

  if (cooling->use_redshift_dependent_UVB == 0) 
    {
      cooling->ChimesGlobalVars.redshift_dependent_UVB_index = -1; 
      cooling->ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0; 
    }
  else 
    {
      if (cooling->UV_field_flag == 0) 
	error("CHIMES ERROR: redshift-dependent UVB has been switched on, but UV_field_flag == 0. These options are incompatible.\n"); 
     
      cooling->ChimesGlobalVars.redshift_dependent_UVB_index = 0;  
      
      if (cooling->use_redshift_dependent_UVB == 1)
	cooling->ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0; 
      else if (cooling->use_redshift_dependent_UVB == 2)
	cooling->ChimesGlobalVars.use_redshift_dependent_eqm_tables = 1; 
      else 
	error("CHIMES ERROR: use_redshift_depdent_UVB = %d not recognised. \n", cooling->use_redshift_dependent_UVB); 
    }
  
  if (cooling->UV_field_flag == 0) 
    cooling->ChimesGlobalVars.N_spectra = 0; 
  else if ((cooling->UV_field_flag == 1) || (cooling->UV_field_flag == 2))
    {
      cooling->radiation_field_normalisation_factor = (ChimesFloat) parser_get_opt_param_double(parameter_file, "CHIMESCooling:radiation_field_normalisation_factor", 1.0); 

      if (cooling->UV_field_flag == 1) 
	{
	  cooling->ChimesGlobalVars.N_spectra = 1; 
 
	  parser_get_param_string(parameter_file, "CHIMESCooling:PhotoIonTable", string_buffer); 
	  sprintf(cooling->ChimesGlobalVars.PhotoIonTablePath[0], "%s/%s", chimes_data_dir, string_buffer); 
	}
      else if (cooling->UV_field_flag == 2) 
	{
	  cooling->ChimesGlobalVars.N_spectra = 2; 
      
	  parser_get_param_string(parameter_file, "CHIMESCooling:PhotoIonTable_UVB", string_buffer); 
	  sprintf(cooling->ChimesGlobalVars.PhotoIonTablePath[0], "%s/%s", chimes_data_dir, string_buffer); 
      
	  parser_get_param_string(parameter_file, "CHIMESCooling:PhotoIonTable_ISRF", string_buffer); 
	  sprintf(cooling->ChimesGlobalVars.PhotoIonTablePath[1], "%s/%s", chimes_data_dir, string_buffer); 
	}
    }
  else 
    error("CHIMESCooling: UV_field_flag %d not recognised.", cooling->UV_field_flag);

  if (cooling->Shielding_flag == 0) 
    cooling->ChimesGlobalVars.cellSelfShieldingOn = 0; 
  else if ((cooling->Shielding_flag == 1) || (cooling->Shielding_flag == 2)) 
    cooling->ChimesGlobalVars.cellSelfShieldingOn = 1; 
  else 
    error("CHIMESCooling: Shielding_flag %d not recognised.", cooling->Shielding_flag); 

  if ((cooling->Shielding_flag > 0) || (cooling->UV_field_flag == 2)) 
    {
      /* These parameters are needed if shielding 
       * is switched on, and also by the COLIBRE 
       * radiation field. */ 

      /* Factor to re-scale shielding length. */ 
      cooling->shielding_length_factor = parser_get_opt_param_double(parameter_file, "CHIMESCooling:shielding_length_factor", 1.0); 

      /* Maximum shielding length (in code units). 
       * If negative, do not impose a maximum. */ 
      cooling->max_shielding_length = parser_get_opt_param_double(parameter_file, "CHIMESCooling:max_shielding_length", -1.0); 
    }

  /* The redshift of hydrogren reionisation 
   * will be needed if we use a redshift-dependent 
   * UVB, and by the COLIBRE cooling tables if 
   * we use hybrid cooling. */ 
  cooling->ChimesGlobalVars.reionisation_redshift = (ChimesFloat) parser_get_param_float(parameter_file, "CHIMESCooling:H_reion_z"); 

  /* Parameters used for the COLIBRE ISRF and 
   * shielding length. These have just been 
   * hard-coded for now - the values are taken
   * from Ploeckinger et al. (in prep). */ 
  cooling->N_H0 = 3.65e20; 

  /* Flag to determine how to set initial 
   * CHIMES abundances: 
   * 0 -- Set each element to one ionisation 
   *      state, determined by the InitIonState 
   *      parameter. 
   * 1 -- Read abundances from eqm abundance tables. 
   * 2 -- Compute initial equilibrium abundances. 
   */ 
  cooling->init_abundance_mode = parser_get_param_int(parameter_file, "CHIMESCooling:init_abundance_mode"); 
  if (!((cooling->init_abundance_mode == 0) || (cooling->init_abundance_mode == 1) || (cooling->init_abundance_mode == 2)))
    error("CHIMESCooling: init_abundance_mode %d not recognised.", cooling->init_abundance_mode); 

  if (cooling->init_abundance_mode == 0) 
    cooling->InitIonState = parser_get_param_int(parameter_file, "CHIMESCooling:InitIonState"); 
  else 
    cooling->InitIonState = 0; 

  /* Cosmic ray ionisation rate of HI. */ 
  cooling->cosmic_ray_rate = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:cosmic_ray_rate"); 

  /* CHIMES tolerance parameters */ 
  cooling->ChimesGlobalVars.relativeTolerance = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:relativeTolerance"); 
  cooling->ChimesGlobalVars.absoluteTolerance = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:absoluteTolerance"); 
  cooling->ChimesGlobalVars.thermalAbsoluteTolerance = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:thermalAbsoluteTolerance"); 
  cooling->ChimesGlobalVars.explicitTolerance = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:explicitTolerance"); 
  cooling->ChimesGlobalVars.scale_metal_tolerances = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:scale_metal_tolerances"); 

  /* Maximum temperature for the molecular network */ 
  cooling->ChimesGlobalVars.T_mol = (ChimesFloat) parser_get_param_double(parameter_file, "CHIMESCooling:T_mol"); 
  
  /* For now, the CMB temperature in CHIMES is just set 
   * to the present day value. For cosmological runs, this 
   * will need to be updated for the current redshift. */ 

  /* CMB temperature at redshift zero. */ 
  cooling->T_CMB_0 = phys_const->const_T_CMB_0 *
                     units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Initially set the CMB temperature in 
   * CHIMES to the redshift zero value. 
   * This will be updated to the current 
   * redshift in cooling_update(). */ 
  cooling->ChimesGlobalVars.cmb_temperature = (ChimesFloat) cooling->T_CMB_0; 

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

  
  /* Flag to switch on extra debug print 
   * statements in CHIMES. 
   * 0 --> No extra debug prints. 
   * 1 --> If CVode returns a non-zero 
   *       flag (i.e. it produces a 
   *       CVODE error/warning), then 
   *       we print out everything in 
   *       the ChimesGasVars struct. 
   */ 
  cooling->ChimesGlobalVars.chimes_debug = parser_get_param_int(parameter_file, "CHIMESCooling:chimes_debug"); 
  
  /* The following CHIMES parameters do not need 
   * to be modified by the user. These are just 
   * hard-coded for now. */ 
  cooling->ChimesGlobalVars.grain_temperature = 10.0; 
  
  /* Physical velocity divergence isn't easily 
   * accessible, so just run with static 
   * molecular cooling for now. */ 
  cooling->ChimesGlobalVars.StaticMolCooling = 1; 
    
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

  /* Dust depletion factors in the solar neighbourhood. 
   * Taken from Ploeckinger et al. (in prep). */ 
  cooling->f_dust0_C = 0.34385; 
  cooling->f_dust0_O = 0.31766; 
  cooling->f_dust0_Mg = 0.94338; 
  cooling->f_dust0_Si = 0.94492; 
  cooling->f_dust0_Ca = 0.9999; 
  cooling->f_dust0_Fe = 0.99363; 

  /* Parameter that controls whether to reduce 
   * gas-phase metal abundances due to 
   * dust depletion. */ 
  cooling->colibre_metal_depletion = parser_get_param_int(parameter_file, "CHIMESCooling:colibre_metal_depletion"); 

  /* delta log U above the EOS below which 
   * we evolve the chemistry in eqm. */ 
  cooling->delta_logUEOS_apply_eqm = parser_get_param_float(parameter_file, "CHIMESCooling:delta_logUEOS_apply_eqm"); 

  /* Threshold in dt / t_cool above which we 
   * are in the rapid cooling regime. If negative, 
   * we never use this scheme (i.e. always drift 
   * the internal energies). */
  cooling->rapid_cooling_threshold = parser_get_param_double(parameter_file, "CHIMESCooling:rapid_cooling_threshold");

  /* Properties of the HII region model */ 
  cooling->HIIregion_temp = parser_get_param_float(parameter_file, "CHIMESCooling:HIIregion_temperature");
  cooling->HIIregion_ion_state = 1; 

  /* Switch for Hybrid cooling */ 
#ifdef COOLING_CHIMES_HYBRID 
  cooling->ChimesGlobalVars.hybrid_cooling_mode = 1; 
#else 
  cooling->ChimesGlobalVars.hybrid_cooling_mode = 0; 
#endif 

  if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 0) 
    {
      /* Use only the CHIMES network 
       * for cooling. */ 
      
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
    }
  else if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 1) 
    {
      /* Use CHIMES only for H and He. 
       * Metal cooling will be read in 
       * from the Cloudy cooling tables. */ 
      
      /* Since metals are not included in 
       * CHIMES here, we hard-code their 
       * include flags to zero. */ 
      int i; 
      for (i = 0; i < 9; i++) 
	cooling->ChimesGlobalVars.element_included[i] = 0; 

      /* Path to colibre cooling table */ 
      parser_get_param_string(parameter_file, "CHIMESCooling:colibre_table_path", cooling->colibre_table.cooling_table_path); 

      /* Redshift of H-reionisation is 
       * needed by the Colibre table 
       * to determine when to use 
       * the high-redshift bin. */ 
      cooling->colibre_table.H_reion_z = (float) cooling->ChimesGlobalVars.reionisation_redshift; 

      /* Set the S/Si and Ca/Si ratios to 
       * the values already read in for CHIMES. */ 
      cooling->colibre_table.S_over_Si_ratio_in_solar = cooling->S_over_Si_ratio_in_solar; 
      cooling->colibre_table.Ca_over_Si_ratio_in_solar = cooling->Ca_over_Si_ratio_in_solar; 

      /* Store some constants in CGS units */
      const float units_kB[5] = {1, 2, -2, 0, -1};
      const double kB_cgs = phys_const->const_boltzmann_k *
	units_general_cgs_conversion_factor(us, units_kB);
      cooling->colibre_table.log10_kB_cgs = log10(kB_cgs);

      /* Read the Colibre table. */ 
      message("Reading Colibre cooling table."); 
      read_cooling_header(&(cooling->colibre_table));
      read_cooling_tables(&(cooling->colibre_table));
    }
  else 
    error("CHIMES ERROR: hybrid_cooling mode %d not recognised. Allowed values are 0 (full CHIMES network) or 1 (Only H+He in CHIMES; metals from COLIBRE tables).", cooling->ChimesGlobalVars.hybrid_cooling_mode);

  /* Set redshift to a very high value, just 
   * while we initialise the CHIMES module. 
   * It will be set to the correct redshift in 
   * the cooling_update() routine. */ 
  cooling->ChimesGlobalVars.redshift = 1000.0; 

  /* Initialise the CHIMES module. */ 
  message("Initialising CHIMES cooling module."); 
  init_chimes(&cooling->ChimesGlobalVars); 

  if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 1) 
    {
      /* Create data structure for hybrid cooling, 
       * and store pointer to the Colibre table. */ 
      cooling->ChimesGlobalVars.hybrid_data = (void *) malloc(sizeof(struct global_hybrid_data_struct)); 
      struct global_hybrid_data_struct *myData; 
      myData = (struct global_hybrid_data_struct *) cooling->ChimesGlobalVars.hybrid_data; 
      myData->table = &(cooling->colibre_table); 

      /* Set the hybrid cooling function pointers. */
      cooling->ChimesGlobalVars.hybrid_cooling_fn = &colibre_metal_cooling_rate_temperature; 
      cooling->ChimesGlobalVars.allocate_gas_hybrid_data_fn = &chimes_allocate_gas_hybrid_data; 
      cooling->ChimesGlobalVars.free_gas_hybrid_data_fn = &chimes_free_gas_hybrid_data; 
    }
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling #cooling_function_data struct.
 */
void cooling_print_backend(const struct cooling_function_data *cooling) {

  if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 0) 
    message("Cooling function is 'CHIMES'.");
  else 
    message("Cooling function is 'CHIMES-HYBRID'.");
}

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The #space containing all the particles.
 */
void cooling_update(const struct cosmology* cosmo,
		    struct cooling_function_data* cooling,
		    struct space* s) {
  /* Update redshift stored in ChimesGlobalVars. */ 
  cooling->ChimesGlobalVars.redshift = cosmo->z; 

  /* Update T_CMB in CHIMES to current redshift. */ 
  cooling->ChimesGlobalVars.cmb_temperature = (ChimesFloat) cooling->T_CMB_0 * (1.0 + cosmo->z); 

  /* Update redshift-dependent UVB. */ 
  if (cooling->use_redshift_dependent_UVB) 
    interpolate_redshift_dependent_UVB(&(cooling->ChimesGlobalVars)); 
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
 * @param dt_cgs The cooling time-step of this particle.
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
  float dimension_G[5] = {-1, 3, -2, 0, 0}; 
  double newton_G_cgs = phys_const->const_newton_G * units_general_cgs_conversion_factor(us, dimension_G); 
  double proton_mass_cgs = phys_const->const_proton_mass * units_cgs_conversion_factor(us, UNIT_CONV_MASS); 

  double mu = chimes_mu(cooling, p, xp); 
  
  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_floor = gas_internal_energy_from_entropy(rho, A_floor);

  double u_actual, T_floor; 
  
  if (u_cgs < u_floor) 
    {
      /* Particle is below the entropy floor. 
       * Set internal energy to the floor. 
       * Chemistry will be evolved in equilibrium. */  
      u_actual = u_floor; 
      ChimesGasVars->ForceEqOn = 1; 
      T_floor = u_floor * hydro_gamma_minus_one * proton_mass_cgs * mu / boltzmann_k_cgs; 
    }
  else if (u_cgs < pow(10.0, cooling->delta_logUEOS_apply_eqm) * u_floor) 
    {
      /* Particle is above the entropy floor, but 
       * close enough that we will need to evolve 
       * the chemistry in equilibrium. */ 
      u_actual = u_cgs; 
      ChimesGasVars->ForceEqOn = 1; 
      T_floor = u_floor * hydro_gamma_minus_one * proton_mass_cgs * mu / boltzmann_k_cgs; 
    }
  else 
    {
      /* Particle is well above the entropy floor. 
       * Evolve chemistry as usual, according to the 
       * user-provided parameter. */ 
      u_actual = u_cgs; 
      ChimesGasVars->ForceEqOn = cooling->ChemistryEqmMode; 

      /* Set T_floor to minimal_temperature, not the 
       * entropy floor. When evolving chemistry in 
       * non-eq, a high T_floor can slow down the 
       * integration. Safer to evolve without and 
       * then re-impose entropy floor afterwards. */
      T_floor = hydro_properties->minimal_temperature; 
    }

  ChimesGasVars->temperature = (ChimesFloat) u_actual * hydro_gamma_minus_one * proton_mass_cgs * mu / boltzmann_k_cgs; 

#if defined(CHEMISTRY_COLIBRE) || defined(CHEMISTRY_EAGLE) 
  float const *metal_fraction = chemistry_get_metal_mass_fraction_for_cooling(p); 
  ChimesFloat XH = (ChimesFloat) metal_fraction[chemistry_element_H]; 

  ChimesGasVars->metallicity = 0.0; 
  float totmass = 0.0, metalmass = 0.0; 
  for (enum colibre_cooling_element elem = element_H; elem < element_OA; elem++) 
    {
      if ((elem != element_S) && (elem != element_Ca)) 
	{
	  totmass += metal_fraction[element_from_table_to_code(elem)]; 
	  if ((elem != element_H) && (elem != element_He)) 
	    metalmass += metal_fraction[element_from_table_to_code(elem)]; 
	}
      else if (elem == element_S) 
	{
	  totmass += metal_fraction[element_from_table_to_code(element_Si)] * cooling->S_over_Si_ratio_in_solar * cooling->S_solar_mass_fraction / cooling->Si_solar_mass_fraction; 
	  metalmass += metal_fraction[element_from_table_to_code(element_Si)] * cooling->S_over_Si_ratio_in_solar * cooling->S_solar_mass_fraction / cooling->Si_solar_mass_fraction; 
	}
      else if (elem == element_Ca) 
	{
	  totmass += metal_fraction[element_from_table_to_code(element_Si)] * cooling->Ca_over_Si_ratio_in_solar * cooling->Ca_solar_mass_fraction / cooling->Si_solar_mass_fraction; 
	  metalmass += metal_fraction[element_from_table_to_code(element_Si)] * cooling->Ca_over_Si_ratio_in_solar * cooling->Ca_solar_mass_fraction / cooling->Si_solar_mass_fraction; 
	}
    }
  ChimesGasVars->metallicity = (ChimesFloat) (metalmass / totmass) / cooling->Zsol; 
#else 
  /* Without COLIBRE or EAGLE chemistry, 
   * the metal abundances are unavailable. 
   * Set to primordial abundances. */ 
  ChimesFloat XH = 0.75; 
  ChimesGasVars->metallicity = 0.0; 
#endif  // CHEMISTRY_COLIBRE || CHEMISTRY_EAGLE 

  ChimesGasVars->dust_ratio = ChimesGasVars->metallicity; 

  ChimesFloat nH = (ChimesFloat) hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass; 
  ChimesGasVars->nH_tot = nH * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY); 

  ChimesGasVars->TempFloor = (ChimesFloat) T_floor;
  ChimesGasVars->cr_rate = cooling->cosmic_ray_rate; 
  ChimesGasVars->hydro_timestep = (ChimesFloat) dt_cgs; 

  ChimesGasVars->constant_heating_rate = 0.0; 
  ChimesGasVars->ThermEvolOn = cooling->ThermEvolOn; 
  ChimesGasVars->divVel = 0.0; 

  /* N_ref is used by both the COLIBRE ISRF
   * and the COLIBRE shielding length, and to 
   * scale metal depletion. Taken from 
   * Ploeckinger et al. (in prep). */ 
  double N_ref = 0.0; 
  if ((cooling->UV_field_flag == 2) || (cooling->Shielding_flag == 2) || (cooling->colibre_metal_depletion == 1)) 
    N_ref = calculate_colibre_N_ref(phys_const, us, cosmo, cooling, p, xp, mu); 

  if (cooling->colibre_metal_depletion == 1) 
    {
      /* Scale dust_ratio by N_ref, but 
       * only if N_ref < N_H0 */ 
      if (N_ref < cooling->N_H0) 
	ChimesGasVars->dust_ratio *= pow(N_ref / cooling->N_H0, 1.4); 
    }
      
  if (cooling->UV_field_flag == 1) 
    {
      /* Single, constant radiation field. */ 

      /* Copy over spectrum parameters from 
       * chimes tables. */ 
      ChimesGasVars->G0_parameter[0] = chimes_table_spectra.G0_parameter[0]; 
      ChimesGasVars->H2_dissocJ[0] = chimes_table_spectra.H2_dissocJ[0]; 
      ChimesGasVars->isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0]; 
      ChimesGasVars->isotropic_photon_density[0] *= cooling->radiation_field_normalisation_factor; 
    }
  else if (cooling->UV_field_flag == 2) 
    {
      /* COLIBRE radiation field */ 
      ChimesFloat J_over_J0; 

      /* Extra-galactic UVB */ 
      ChimesGasVars->G0_parameter[0] = chimes_table_spectra.G0_parameter[0]; 
      ChimesGasVars->H2_dissocJ[0] = chimes_table_spectra.H2_dissocJ[0]; 
      ChimesGasVars->isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0]; 

      /* ISRF */ 
      ChimesGasVars->G0_parameter[1] = chimes_table_spectra.G0_parameter[1]; 
      ChimesGasVars->H2_dissocJ[1] = chimes_table_spectra.H2_dissocJ[1]; 
      J_over_J0 = cooling->radiation_field_normalisation_factor * pow(N_ref / cooling->N_H0, 1.4); 
      ChimesGasVars->isotropic_photon_density[1] = chimes_table_spectra.isotropic_photon_density[1]; 

      /* low-density cut-off before reionisation */ 
      if ((cooling->ChimesGlobalVars.redshift > cooling->ChimesGlobalVars.reionisation_redshift) && (J_over_J0 > 0.0)) 
	ChimesGasVars->isotropic_photon_density[1] *= pow(10.0, -20.0 - ((-20.0 - log10(J_over_J0)) / (1.0 + exp(-2.0 * (log10(ChimesGasVars->nH_tot) + 4.0))))); 
      else 
	ChimesGasVars->isotropic_photon_density[1] *= J_over_J0; 
      
      /* Scale cr_rate by N_ref */ 
      ChimesGasVars->cr_rate *= cooling->radiation_field_normalisation_factor * pow(N_ref / cooling->N_H0, 1.4); 
    }

  if (cooling->Shielding_flag == 0) 
    ChimesGasVars->cell_size = 0; 
  else if (cooling->Shielding_flag == 1) 
    {
      /* Jeans length */ 
      ChimesGasVars->cell_size = sqrt(M_PI * hydro_gamma * boltzmann_k_cgs * ChimesGasVars->temperature / (mu * newton_G_cgs * (proton_mass_cgs * ChimesGasVars->nH_tot / XH) * proton_mass_cgs)); 
      ChimesGasVars->cell_size *= cooling->shielding_length_factor; 

      if (cooling->max_shielding_length > 0.0) 
	{
	  double max_shielding_length_cgs = cooling->max_shielding_length * units_cgs_conversion_factor(us, UNIT_CONV_LENGTH); 
	  if (ChimesGasVars->cell_size > max_shielding_length_cgs) 
	    ChimesGasVars->cell_size = max_shielding_length_cgs; 
	}
    }
  else if (cooling->Shielding_flag == 2) 
    ChimesGasVars->cell_size = cooling->shielding_length_factor * N_ref / ChimesGasVars->nH_tot; 

  /* Doppler broadening parameter, for 
   * H2 self-shielding, is hard-coded 
   * to 7.1 km/s for now. This is a 
   * typical value for GMCs in the 
   * Milky Way. */ 
  ChimesGasVars->doppler_broad = 7.1; 

  ChimesGasVars->InitIonState = cooling->InitIonState; 

  /* If using hybrid cooling, we need to 
   * set the abundance_ratio array using 
   * the corresponding routine from COLIBRE. */
  if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 1) 
    {
      struct gas_hybrid_data_struct *myData; 
      myData = (struct gas_hybrid_data_struct *) ChimesGasVars->hybrid_data; 
      abundance_ratio_to_solar(p, &(cooling->colibre_table), myData->abundance_ratio); 
    }
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
 * @param phys_const #phys_const data structure.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @param ChimesGasVars CHIMES gasVariables structure. 
 * @param mode Set to zero if particle not fully initialised.
 */ 
void chimes_update_element_abundances(const struct phys_const *phys_const,
				      const struct unit_system *us,
				      const struct cosmology *cosmo,
				      const struct cooling_function_data *cooling,
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

  double mu; 
  if (mode == 0) 
    {
      /* When chimes_update_element_abundances() 
       * is called on a particle that is being 
       * initialised for the first time, the 
       * chimes_abundance array has not yet been
       * initialised. We therefore set mu = 1. */ 
      mu = 1.0; 
    }
  else 
    mu = chimes_mu(cooling, p, xp); 

  double N_ref = calculate_colibre_N_ref(phys_const, us, cosmo, cooling, p, xp, mu); 
  double factor = min(pow(N_ref / cooling->N_H0, 1.4), 1.0); 
  if (cooling->colibre_metal_depletion == 1) 
    {
      /* Reduce gas-phase metal abundances 
       * due to dust depletion. */ 
      ChimesGasVars->element_abundances[1] *= (1.0 - (cooling->f_dust0_C * factor)); 
      ChimesGasVars->element_abundances[3] *= (1.0 - (cooling->f_dust0_O * factor)); 
      ChimesGasVars->element_abundances[5] *= (1.0 - (cooling->f_dust0_Mg * factor)); 
      ChimesGasVars->element_abundances[6] *= (1.0 - (cooling->f_dust0_Si * factor)); 
      ChimesGasVars->element_abundances[8] *= (1.0 - (cooling->f_dust0_Ca * factor)); 
      ChimesGasVars->element_abundances[9] *= (1.0 - (cooling->f_dust0_Fe * factor)); 
    }      
  
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
 * @param time Time since Big Bang 
 */
void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm, const double time) {

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
    ChimesGasVars.abundances[i] = (ChimesFloat) xp->cooling_data.chimes_abundances[i]; 

  /* Update element abundances from metal mass 
   * fractions. We need to do this here, and not 
   * later on in chimes_update_gas_vars(), because 
   * the element abundances need to be set in 
   * ChimesGasVars before we can calculate the 
   * mean molecular weight. */ 
  chimes_update_element_abundances(phys_const, us, cosmo, cooling, p, xp, &ChimesGasVars, 1); 

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
  const double u_0_cgs = u_0 * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  const double dt_cgs = dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Update the ChimesGasVars structure with the 
   * particle's thermodynamic variables. */ 
  chimes_update_gas_vars(u_0_cgs, phys_const, us, cosmo, hydro_properties, floor_props, cooling, p, xp, &ChimesGasVars, dt_cgs); 

  /* check if the particle is in an HII region. If it is, we 
   * immediately heat it up to 1e4 K (if required), and it 
   * will subsequently be evolved with equilibrium cooling. */ 
  int HII_flag = 0; 
  if ((time <= xp->tracers_data.HIIregion_timer_gas) &&
      (xp->tracers_data.HIIregion_timer_gas > 0.)) {
    ChimesGasVars.temperature = chimes_max(ChimesGasVars.temperature, (ChimesFloat) cooling->HIIregion_temp); 
    ChimesGasVars.TempFloor = chimes_max(ChimesGasVars.TempFloor, (ChimesFloat) cooling->HIIregion_temp); 
    ChimesGasVars.ForceEqOn = 1; 
    HII_flag = 1; 
  } else if ((time > xp->tracers_data.HIIregion_timer_gas) &&
             (xp->tracers_data.HIIregion_timer_gas > 0.)) {
    xp->tracers_data.HIIregion_timer_gas = -1.;
    xp->tracers_data.HIIregion_starid = -1;
  }
  
  /* Call CHIMES to integrate the chemistry 
   * and cooling over the time-step. */ 
  chimes_network(&ChimesGasVars, &ChimesGlobalVars); 

  if (HII_flag == 1) 
    {
      /* If particle is in an HII region, manually 
       * set it to be ionised, according to the 
       * HIIregion_ion_state parameter. */ 
      ChimesGasVars.InitIonState = cooling->HIIregion_ion_state; 
      initialise_gas_abundances(&ChimesGasVars, &ChimesGlobalVars); 
    }

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

  /* Convert back to internal units */ 
  double u_final = u_final_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS); 

  /* We now need to check that we are not going to go below any of the limits */

  /* Recompute new minimal internal energy 
  * from the new mean molecular weight. */
  minimal_internal_energy = hydro_properties->minimal_temperature; 
  minimal_internal_energy *= hydro_one_over_gamma_minus_one; 
  minimal_internal_energy *= (phys_const->const_boltzmann_k / phys_const->const_proton_mass); 
  minimal_internal_energy /= mu; 

  u_final = max(u_final, minimal_internal_energy); 

  /* Limit imposed by the entropy floor */
  const double A_floor = entropy_floor(p, cosmo, floor_props);
  const double rho = hydro_get_physical_density(p, cosmo);
  const double u_floor = gas_internal_energy_from_entropy(rho, A_floor);

  u_final = max(u_final, u_floor); 

  if ((u_final == u_floor) && (ChimesGasVars.ForceEqOn == 0)) 
    {
      /* The particle has reached the entropy 
       * floor, but it was evolved with non-
       * eqm chemistry (meaning that it started 
       * from a long way above the EOS). We 
       * now need to re-set its abundance array 
       * to be in chemical equilibrium. */ 
      const double u_floor_cgs = u_floor * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS); 
      chimes_update_gas_vars(u_floor_cgs, phys_const, us, cosmo, hydro_properties, floor_props, cooling, p, xp, &ChimesGasVars, dt_cgs); 
      ChimesGasVars.ForceEqOn = 1; 
      ChimesGasVars.ThermEvolOn = 0; 
      chimes_network(&ChimesGasVars, &ChimesGlobalVars);
    }

  /* Expected change in energy over the next kick step
     (assuming no change in dt) */
  const double delta_u = u_final - max(u_start, u_floor);

  /* Determine if we are in the slow- or rapid-cooling regime,
   * by comparing dt / t_cool to the rapid_cooling_threshold.
   *
   * Note that dt / t_cool = fabs(delta_u) / u_start. */
  const double dt_over_t_cool = fabs(delta_u) / max(u_start, u_floor);

   /* If rapid_cooling_threshold < 0, always use the slow-cooling
   * regime. */
  if ((cooling->rapid_cooling_threshold >= 0.0) &&
      (dt_over_t_cool >= cooling->rapid_cooling_threshold)) {

    /* Rapid-cooling regime. */ 

    /* Update the particle's u and du/dt */
    hydro_set_physical_internal_energy(p, xp, cosmo, u_final);
    hydro_set_drifted_physical_internal_energy(p, cosmo, u_final);
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.);

  } else {
    /* Slow-cooling regime. */ 

    /* Update du/dt so that we can subsequently drift internal energy. */
    const float cooling_du_dt = delta_u / dt_therm;

    /* Update the internal energy time derivative */
    hydro_set_physical_internal_energy_dt(p, cosmo, cooling_du_dt);
  }

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * (u_final - u_0); 
  
  /* Copy abundances from ChimesGasVars back to xp. */ 
  for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
    xp->cooling_data.chimes_abundances[i] = (double) ChimesGasVars.abundances[i]; 

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
 * @brief Split the cooling content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
void cooling_split_part(struct part *p, struct xpart *xp, double n) {

  xp->cooling_data.radiated_energy /= n;
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
    denominator += xp->cooling_data.chimes_abundances[i]; 
  
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

/**
 * @brief Calculate the N_ref column density. 
 * 
 * This routine returns the column density N_ref 
 * as defined in Ploeckinger et al. (in prep), 
 * which is used to scale the ISRF, cosmic rays, 
 * dust depletion and shielding column density 
 * in COLIBRE. 
 *
 * @param phys_const #phys_const data structure.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 * @param mu Mean molecular weight. 
 */
double calculate_colibre_N_ref(const struct phys_const *phys_const,
			       const struct unit_system *us,
			       const struct cosmology *cosmo,
			       const struct cooling_function_data *cooling,
			       struct part *restrict p, struct xpart* restrict xp, 
			       const double mu) {
  /* Physical constants that we will 
   * need, in cgs units */ 
  float dimension_k[5] = {1, 2, -2, 0, -1}; 
  double boltzmann_k_cgs = phys_const->const_boltzmann_k * units_general_cgs_conversion_factor(us, dimension_k); 
  float dimension_G[5] = {-1, 3, -2, 0, 0}; 
  double newton_G_cgs = phys_const->const_newton_G * units_general_cgs_conversion_factor(us, dimension_G); 
  double proton_mass_cgs = phys_const->const_proton_mass * units_cgs_conversion_factor(us, UNIT_CONV_MASS);

  /* Parameters that define N_ref. 
   * Taken from Ploeckinger et al. (in prep). */ 
  const double N_max = 1.0e24;         // cgs 
  const double N_min = 3.08567758e15;  // cgs 
  const double l_max_cgs = cooling->max_shielding_length * units_cgs_conversion_factor(us, UNIT_CONV_LENGTH); 
  const double log_T_min = 3.0; 
  const double log_T_max = 5.0; 

#if defined(CHEMISTRY_COLIBRE) || defined(CHEMISTRY_EAGLE) 
  float const *metal_fraction = chemistry_get_metal_mass_fraction_for_cooling(p); 
  ChimesFloat XH = (ChimesFloat) metal_fraction[chemistry_element_H]; 
#else 
  /* Without COLIBRE or EAGLE chemistry, 
   * the metal abundances are unavailable. 
   * Set to primordial abundances. */ 
  ChimesFloat XH = 0.75; 
#endif  // CHEMISTRY_COLIBRE || CHEMISTRY_EAGLE 

  /* Density*/ 
  const double nH = hydro_get_physical_density(p, cosmo) * XH / phys_const->const_proton_mass; 
  const double nH_cgs = nH * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY); 

  /* Internal energy */
  const double u = hydro_get_physical_internal_energy(p, xp, cosmo); 
  const double u_cgs = u * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS); 

  /* Temperature */ 
  const double temperature = u_cgs * hydro_gamma_minus_one * proton_mass_cgs * mu / boltzmann_k_cgs; 

  /* Jeans column density */ 
  double N_J = nH_cgs * sqrt(hydro_gamma * boltzmann_k_cgs * temperature / (mu *newton_G_cgs * (proton_mass_cgs * nH_cgs / XH) * proton_mass_cgs)); 

  double N_ref_prime = chimes_min(N_J, N_max); 
      
  if (l_max_cgs > 0.0) 
    N_ref_prime = chimes_min(N_ref_prime, l_max_cgs * nH_cgs); 

  double N_ref = pow(10.0, log10(N_ref_prime) - ((log10(N_ref_prime) - log10(N_min)) / (1.0 + exp(-5.0 * (log10(temperature) - ((log_T_min + log_T_max) / 2.0))))));
  
  return N_ref; 
}

/**
 * @brief Write a cooling struct to the given FILE as a stream of bytes.
 *
 * @param cooling the struct
 * @param stream the file stream
 */
void cooling_struct_dump(const struct cooling_function_data* cooling, 
			 FILE* stream) {

  /* Zero all pointers in the colibre_table within 
   * the cooling_function_data struct. */ 
  struct cooling_function_data cooling_copy = *cooling; 
  cooling_copy.colibre_table.Tcooling = NULL; 
  cooling_copy.colibre_table.Theating = NULL; 
  cooling_copy.colibre_table.Telectron_fraction = NULL; 
  cooling_copy.colibre_table.Redshifts = NULL; 
  cooling_copy.colibre_table.nH = NULL; 
  cooling_copy.colibre_table.Temp = NULL; 
  cooling_copy.colibre_table.Metallicity = NULL; 
  cooling_copy.colibre_table.LogAbundances = NULL; 
  cooling_copy.colibre_table.Abundances = NULL; 
  cooling_copy.colibre_table.Abundances_inv = NULL; 
  cooling_copy.colibre_table.atomicmass = NULL; 
  cooling_copy.colibre_table.atomicmass_inv = NULL; 
  cooling_copy.colibre_table.LogMassFractions = NULL; 
  cooling_copy.colibre_table.MassFractions = NULL; 
  cooling_copy.colibre_table.Zsol = NULL; 
  cooling_copy.colibre_table.Zsol_inv = NULL; 
  cooling_copy.ChimesGlobalVars.hybrid_data = NULL; 

  restart_write_blocks((void *) &cooling_copy,
                       sizeof(struct cooling_function_data), 1, stream,
                       "cooling", "cooling function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param cooling the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void cooling_struct_restore(struct cooling_function_data* cooling,
			    FILE* stream, const struct cosmology* cosmo) {

  restart_read_blocks((void *) cooling, sizeof(struct cooling_function_data), 1,
                      stream, NULL, "cooling function");

  if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 1) 
    {
      /* Read the Colibre table. */ 
      message("Reading Colibre cooling table."); 
      read_cooling_header(&(cooling->colibre_table));
      read_cooling_tables(&(cooling->colibre_table));
    }

  /* Initialise the CHIMES module. */ 
  message("Initialising CHIMES cooling module."); 
  init_chimes(&cooling->ChimesGlobalVars); 

  if (cooling->ChimesGlobalVars.hybrid_cooling_mode == 1) 
    {
      /* Create data structure for hybrid cooling, 
       * and store pointer to the Colibre table. */ 
      cooling->ChimesGlobalVars.hybrid_data = (void *) malloc(sizeof(struct global_hybrid_data_struct)); 
      struct global_hybrid_data_struct *myData; 
      myData = (struct global_hybrid_data_struct *) cooling->ChimesGlobalVars.hybrid_data; 
      myData->table = &(cooling->colibre_table); 

      /* Set the hybrid cooling function pointers. */
      cooling->ChimesGlobalVars.hybrid_cooling_fn = &colibre_metal_cooling_rate_temperature; 
      cooling->ChimesGlobalVars.allocate_gas_hybrid_data_fn = &chimes_allocate_gas_hybrid_data; 
      cooling->ChimesGlobalVars.free_gas_hybrid_data_fn = &chimes_free_gas_hybrid_data; 
    }
}

/**
 * @brief Converts cooling quantities of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 *
 * For this cooling module, this routine is used to set the cooling 
 * properties of the (x-)particles to a valid start state, in particular 
 * the CHIMES abundance array. 
 * 
 * This is controlled by the cooling->init_abundance_mode as follows: 
 * 0 -- Set each element to one ionisation state, determined by the 
 *      ChimesGlobalVars.InitIonState parameter. 
 * 1 -- Read abundances from eqm abundance tables. 
 * 2 -- Compute initial equilibrium abundances. 
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param phys_const #phys_const data structure.
 * @param us Internal system of units data structure.
 * @param floor_props Properties of the entropy floor.
 * @param cooling #cooling_function_data data structure. 
 */
void cooling_convert_quantities(struct part *restrict p, 
				struct xpart *restrict xp,
				const struct cosmology *cosmo, 
				const struct hydro_props *hydro_props, 
				const struct phys_const *phys_const, 
				const struct unit_system* us, 
				const struct entropy_floor_properties *floor_props, 
				const struct cooling_function_data* cooling) {
  struct globalVariables ChimesGlobalVars = cooling->ChimesGlobalVars; 
  struct gasVariables ChimesGasVars; 
  int i; 

  /* Allocate memory to arrays within ChimesGasVars. */
  allocate_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 

  /* Set element abundances from 
   * metal mass fractions. */ 
  chimes_update_element_abundances(phys_const, us, cosmo, cooling, p, xp, &ChimesGasVars, 0); 

  /* Zero particle's radiated energy. */ 
  xp->cooling_data.radiated_energy = 0.f; 

  /* Set initial values for CHIMES 
   * abundance array. */ 
  ChimesGasVars.InitIonState = cooling->InitIonState; 
  initialise_gas_abundances(&ChimesGasVars, &ChimesGlobalVars); 

  if (cooling->init_abundance_mode == 0) 
    {
      // Copy abundances over to xp. 
      for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
	xp->cooling_data.chimes_abundances[i] = (double) ChimesGasVars.abundances[i]; 
    }
  else if ((cooling->init_abundance_mode == 1) || (cooling->init_abundance_mode == 2)) 
    {
      // Copy initial abundances over to xp. 
      for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
	xp->cooling_data.chimes_abundances[i] = (double) ChimesGasVars.abundances[i]; 

      /* Get the particle's internal energy */ 
      double u_0 = hydro_get_physical_internal_energy(p, xp, cosmo); 
      double u_0_cgs = u_0 * units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

      /* If computing the eqm (init_abundance_mode == 2),
       * we will integrate the chemistry ten times for 
       * 1 Gyr per iteration. Multiple iterations 
       * are required so that the shielding column 
       * densities can be updated between each 
       * iteration. If reading from tables, we only 
       * need 1 iteration. */ 
      double dt_cgs = 3.15576e16; 
      int n_iterations; 
      if (cooling->init_abundance_mode == 1) 
	n_iterations= 1; 
      else 
	n_iterations= 10; 
  
      for (i = 0; i < n_iterations; i++) 
	{
	  /* Update element abundances. This 
	   * accounts for the dust depletion 
	   * factors. */ 
	  chimes_update_element_abundances(phys_const, us, cosmo, cooling, p, xp, &ChimesGasVars, 1); 

	  /* Update ChimesGasVars with the particle's 
	   * thermodynamic variables. */ 
	  chimes_update_gas_vars(u_0_cgs, phys_const, us, cosmo, hydro_props, floor_props, cooling, p, xp, &ChimesGasVars, dt_cgs); 
  
	  /* Set temperature evolution off, so that we
	   * compute equilibrium at fixed temperature. */ 
	  ChimesGasVars.ThermEvolOn = 0; 

	  /* Determine whether to use eqm tables 
	   * or compute the eqm state. */ 
	  if (cooling->init_abundance_mode == 1) 
	    ChimesGasVars.ForceEqOn = 1; 
	  else 
	    ChimesGasVars.ForceEqOn = 0; 

	  /* Integrate to chemical equilibrium. */
	  chimes_network(&ChimesGasVars, &ChimesGlobalVars); 
	}

      // Copy final abundances over to xp. 
      for (i = 0; i < ChimesGlobalVars.totalNumberOfSpecies; i++) 
	xp->cooling_data.chimes_abundances[i] = (double) ChimesGasVars.abundances[i]; 
    }
  else 
    error("CHIMESCooling: init_abundance_mode %d not recognised.", cooling->init_abundance_mode); 

  /* Free CHIMES memory. */ 
  free_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 
}
