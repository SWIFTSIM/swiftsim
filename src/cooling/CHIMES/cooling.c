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
  char chimes_eqm_abundance_table[500]; 

  /* read the parameters */
  if (sizeof(ChimesFloat) != sizeof(double)) 
    error("CHIMES ERROR: When we parse the CHIMES parameters, we assume that the ChimesFloat type has been set to double. However, sizeof(ChimesFloat) = %lu bytes, while sizeof(double) = %lu bytes.\n", sizeof(ChimesFloat), sizeof(double)); 


  /* Set paths to CHIMES data files. */ 
  parser_get_param_string(parameter_file, "CHIMESCooling:data_path", chimes_data_dir); 
  sprintf(cooling->ChimesGlobalVars.MainDataTablePath, "%s/chimes_main_data.hdf5", chimes_data_dir); 

  parser_get_param_string(parameter_file, "CHIMESCooling:EqmAbundanceTable", chimes_eqm_abundance_table); 
  sprintf(cooling->ChimesGlobalVars.EqAbundanceTablePath, "%s/EqAbundancesTables/%s", chimes_data_dir, chimes_eqm_abundance_table); 

  /* Define radiation field and shielding options. */ 
  /* Currently supported options: 
   * UV_field_flag == 0: No UV radiation. 
   * 
   * Shielding_flag == 0: No shielding. 
   */ 
  cooling->UV_field_flag = parser_get_param_int(parameter_file, "CHIMESCooling:UV_field_flag"); 
  cooling->Shielding_flag = parser_get_param_int(parameter_file, "CHIMESCooling:Shielding_flag"); 

  if (cooling->UV_field_flag == 0) 
    cooling->ChimesGlobalVars.N_spectra = 0; 
  else 
    error("CHIMESCooling: UV_field_flag %d not recognised.", cooling->UV_field_flag); 

  if (cooling->Shielding_flag == 0) 
    cooling->ChimesGlobalVars.cellSelfShieldingOn = 0; 
  else 
    error("CHIMESCooling: Shielding_flag %d not recognised.", cooling->Shielding_flag); 

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
  cooling->ChimesGlobalVars.StaticMolCooling = 0; 
  cooling->ChimesGlobalVars.InitIonState = 1; 
  cooling->ChimesGlobalVars.grain_temperature = 10.0; 

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
			     const struct part* restrict p,
			     struct xpart* restrict xp) {
  struct globalVariables ChimesGlobalVars = data->ChimesGlobalVars; 
  struct gasVariables ChimesGasVars; 
  ChimesGasVars.abundances = xp->cooling_data.chimes_abundances; 

  /* Get element mass fractions */ 
#if defined(CHEMISTRY_COLIBRE) || defined(CHEMISTRY_EAGLE) 
  float const *metal_fraction = chemistry_get_metal_mass_fraction_for_cooling(p); 
  float XH = metal_fraction[chemistry_element_H]; 
  ChimesGasVars.element_abundances[0] = metal_fraction[chemistry_element_He] / (4.0 * XH); 
  ChimesGasVars.element_abundances[1] = metal_fraction[chemistry_element_C] / (12.0 * XH); 
  ChimesGasVars.element_abundances[2] = metal_fraction[chemistry_element_N] / (14.0 * XH); 
  ChimesGasVars.element_abundances[3] = metal_fraction[chemistry_element_O] / (16.0 * XH); 
  ChimesGasVars.element_abundances[4] = metal_fraction[chemistry_element_Ne] / (20.0 * XH); 
  ChimesGasVars.element_abundances[5] = metal_fraction[chemistry_element_Mg] / (24.0 * XH); 
  ChimesGasVars.element_abundances[6] = metal_fraction[chemistry_element_Si] / (28.0 * XH); 
  ChimesGasVars.element_abundances[9] = metal_fraction[chemistry_element_Fe] / (56.0 * XH); 
  
  ChimesGasVars.element_abundances[7] = metal_fraction[chemistry_element_Si] * data->S_over_Si_ratio_in_solar * (data->S_solar_mass_fraction / data->Si_solar_mass_fraction) / (32.0 * XH); 
  ChimesGasVars.element_abundances[8] = metal_fraction[chemistry_element_Si] * data->Ca_over_Si_ratio_in_solar * (data->Ca_solar_mass_fraction / data->Si_solar_mass_fraction) / (40.0 * XH); 
#else 
  /* Without COLIBRE or EAGLE chemistry, 
   * the metal abundances are unavailable. 
   * Set to primordial abundances. */ 
  ChimesGasVars.element_abundances[0] = 0.0833;  // He 

  int i; 
  for (i = 1; i < 10; i++) 
    ChimesGasVars.element_abundances[i] = 0.0; 
#endif  // CHEMISTRY_COLIBRE || CHEMISTRY_EAGLE 

  initialise_gas_abundances(&ChimesGasVars, &ChimesGlobalVars); 
}
