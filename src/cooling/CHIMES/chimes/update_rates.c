/* This file contains functions to 
 * set up the reaction and 
 * photoionisation rates at the 
 * beginning of the run, based
 * on the given parameters.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chimes_vars.h"
#include "chimes_proto.h"

/* This routine is called at the beginning of the 
 * integration to set all rate coefficients, 
 * including those that will be held fixed during 
 * the integration. */ 
void set_initial_rate_coefficients(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data) 
{ 
  int i, j, NHI_index, NH_eff_index, NHeI_index, NHe_eff_index, NCO_index, NH2_index; 
  ChimesFloat dNHI, dNH_eff, dNHeI, dNHe_eff, dNCO, dNH2;
  ChimesFloat flux, G0, S1, S2, S3, log_NCO, log_NH2; 

  // NOTE: the individual E factors can be <1e-38, so they need to 
  // be in double precision. However, when we take the ratio 
  // (E1 + E2 + E3) / (E4 + E5 + E6) we can cast back to 
  // ChimesFloat 
  double E1, E2, E3, E4, E5, E6; 
  update_rate_coefficients(myGasVars, myGlobalVars, data, 1); 

  if (myGlobalVars->N_spectra > 0) 
    {
      /* The following rate coefficients 
       * do not vary throughout the course 
       * of the integration. They only 
       * need to be set once. */ 

	  // Compute shield factor
      if (myGlobalVars->cellSelfShieldingOn > 0)
	{
	  // photoion_fuv 
	  for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photoion_fuv_shield_factor[i] = (ChimesFloat) exp(-(chimes_table_photoion_fuv.gamma[i] * data.extinction)); 

	  // photoion_euv 
	  chimes_get_table_index(chimes_table_bins.Column_densities, chimes_table_bins.N_Column_densities, log10(chimes_max(data.HI_column, 1.0e-100)), &NHI_index, &dNHI);
	  chimes_get_table_index(chimes_table_bins.Column_densities, chimes_table_bins.N_Column_densities, log10(chimes_max(data.HI_column + (3.0 * data.H2_column), 1.0e-100)), &NH_eff_index, &dNH_eff);
	  chimes_get_table_index(chimes_table_bins.Column_densities, chimes_table_bins.N_Column_densities, log10(chimes_max(data.HeI_column, 1.0e-100)), &NHeI_index, &dNHeI);
	  chimes_get_table_index(chimes_table_bins.Column_densities, chimes_table_bins.N_Column_densities, log10(chimes_max(data.HeI_column + (0.75 * data.HeII_column), 1.0e-100)), &NHe_eff_index, &dNHe_eff);

	  for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	    {
	      for (j = 0; j < myGlobalVars->N_spectra; j++) 
		{
		  if (chimes_table_photoion_euv.E_thresh[i] < 15.4) 
		    S1 = pow(10.0, chimes_interpol_1d(chimes_table_photoion_euv.shieldFactor_1D[i][j][0], NHI_index, dNHI)); 
		  else 
		    S1 = 0.0; 

		  if (chimes_table_photoion_euv.E_thresh[i] < 54.42) 
		    S2 = pow(10.0, chimes_interpol_2d(chimes_table_photoion_euv.shieldFactor_2D[i][j][0], NH_eff_index, NHeI_index, dNH_eff, dNHeI)); 
		  else 
		    S2 = 0.0; 

		  S3 = pow(10.0, chimes_interpol_2d(chimes_table_photoion_euv.shieldFactor_2D[i][j][1], NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff)); 
		  
		  data.chimes_current_rates->photoion_euv_shield_factor[i][j] = S1 + S2 + S3; 
		}
	    }

	  if (myGasVars->ThermEvolOn) 
	    {
	      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
		{
		  for (j = 0; j < myGlobalVars->N_spectra; j++) 
		    {
		      if (chimes_table_photoion_euv.E_thresh[i] < 15.4) 
			{
			  E1 = pow(10.0, (double) chimes_interpol_1d(chimes_table_photoion_euv.shieldFactor_1D[i][j][1], NHI_index, dNHI)); 
			  E4 = pow(10.0, (double) chimes_interpol_1d(chimes_table_photoion_euv.shieldFactor_1D[i][j][2], NHI_index, dNHI)); 
			}
		      else 
			{
			  E1 = 0.0; 
			  E4 = 0.0; 
			}

		      if (chimes_table_photoion_euv.E_thresh[i] < 54.42) 
			{
			  E2 = pow(10.0, (double) chimes_interpol_2d(chimes_table_photoion_euv.shieldFactor_2D[i][j][2], NH_eff_index, NHeI_index, dNH_eff, dNHeI)); 
			  E5 = pow(10.0, (double) chimes_interpol_2d(chimes_table_photoion_euv.shieldFactor_2D[i][j][4], NH_eff_index, NHeI_index, dNH_eff, dNHeI)); 
			}
		      else 
			{
			  E2 = 0.0; 
			  E5 = 0.0; 
			}

		      E3 = pow(10.0, (double) chimes_interpol_2d(chimes_table_photoion_euv.shieldFactor_2D[i][j][3], NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff)); 
		      E6 = pow(10.0, (double) chimes_interpol_2d(chimes_table_photoion_euv.shieldFactor_2D[i][j][5], NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff)); 
		  
		      data.chimes_current_rates->photoion_euv_epsilon[i][j] = (ChimesFloat) ((E1 + E2 + E3) / (E4 + E5 + E6)); 
		    }
		}
	    }

	  // photodissoc_group1 
	  for (i = 0; i < chimes_table_photodissoc_group1.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photodissoc_group1_shield_factor[i] = (ChimesFloat) exp(-(chimes_table_photodissoc_group1.gamma[i] * data.extinction)); 

	  if (data.mol_flag_index == 1) 
	    {
	      // photodissoc_group2 
	      // All reactions use the same shield_factor 
	      if (data.extinction > 15) 
		data.chimes_current_rates->photodissoc_group2_shield_factor = (ChimesFloat) exp(-(chimes_table_photodissoc_group2.gamma_coeff[0] * data.extinction)); 
	      else 
		data.chimes_current_rates->photodissoc_group2_shield_factor = (ChimesFloat) exp(-(chimes_table_photodissoc_group2.gamma_coeff[1] * data.extinction) + (chimes_table_photodissoc_group2.gamma_coeff[2] * pow(data.extinction, 2.0))); 

	      // CO_photodissoc 
	      log_NCO = (ChimesFloat) log10(chimes_max(data.CO_column, 1.0e-100)); 
	      log_NH2 = (ChimesFloat) log10(chimes_max(data.H2_column, 1.0e-100)); 
	      chimes_get_table_index(chimes_table_bins.COself_column_densities, chimes_table_bins.N_COself_column_densities, log_NCO, &NCO_index, &dNCO); 
	      chimes_get_table_index(chimes_table_bins.H2CO_column_densities, chimes_table_bins.N_H2CO_column_densities, log_NH2, &NH2_index, &dNH2); 
	  
	      for (i = 0; i < chimes_table_CO_photodissoc.N_reactions[data.mol_flag_index]; i++) 
		{
		  data.chimes_current_rates->CO_photodissoc_shield_factor[i] = pow(10.0, chimes_interpol_2d(chimes_table_CO_photodissoc.self_shielding[i], NCO_index, NH2_index, dNCO, dNH2)); 
		  data.chimes_current_rates->CO_photodissoc_shield_factor[i] *= exp(- chimes_table_CO_photodissoc.gamma[i] * data.extinction);
		}
	    }
	}
      else 
	{
	  for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photoion_fuv_shield_factor[i] = 1.0; 

	  for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	    {
	      for (j = 0; j < myGlobalVars->N_spectra; j++) 
		data.chimes_current_rates->photoion_euv_shield_factor[i][j] = 1.0; 
	    }

	  if (myGasVars->ThermEvolOn) 
	    {
	      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
		{
		  for (j = 0; j < myGlobalVars->N_spectra; j++) 
		    {

		      /* self-shielding is switched off, so 
		       * take epsilon from the zeroth entry 
		       * of the shieldFactor tables. */
		      if (chimes_table_photoion_euv.E_thresh[i] < 15.4) 
			{
			  E1 = pow(10.0, (double) chimes_table_photoion_euv.shieldFactor_1D[i][j][1][0]); 
			  E4 = pow(10.0, (double) chimes_table_photoion_euv.shieldFactor_1D[i][j][2][0]); 
			}
		      else 
			{
			  E1 = 0.0; 
			  E4 = 0.0; 
			}

		      if (chimes_table_photoion_euv.E_thresh[i] < 54.42) 
			{
			  E2 = pow(10.0, (double) chimes_table_photoion_euv.shieldFactor_2D[i][j][2][0][0]); 
			  E5 = pow(10.0, (double) chimes_table_photoion_euv.shieldFactor_2D[i][j][4][0][0]); 
			}
		      else 
			{
			  E2 = 0.0; 
			  E5 = 0.0; 
			}

		      E3 = pow(10.0, (double) chimes_table_photoion_euv.shieldFactor_2D[i][j][3][0][0]); 
		      E6 = pow(10.0, (double) chimes_table_photoion_euv.shieldFactor_2D[i][j][5][0][0]); 
		  
		      data.chimes_current_rates->photoion_euv_epsilon[i][j] = (ChimesFloat) ((E1 + E2 + E3) / (E4 + E5 + E6)); 
		    }
		}
	    }
	  
	  for (i = 0; i < chimes_table_photodissoc_group1.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photodissoc_group1_shield_factor[i] = 1.0; 

	  if (data.mol_flag_index == 1) 
	    {
	      data.chimes_current_rates->photodissoc_group2_shield_factor = 1.0; 
	      
	      for (i = 0; i < chimes_table_CO_photodissoc.N_reactions[data.mol_flag_index]; i++) 
		data.chimes_current_rates->CO_photodissoc_shield_factor[i] = 1.0; 
	    }
	}

      // Zero the rate coefficients 
      for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->photoion_fuv_rate_coefficient[i] = 0.0; 

      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->photoion_euv_rate_coefficient[i] = 0.0; 

      if (myGasVars->ThermEvolOn) 
	{
	  for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photoion_fuv_heat_rate[i] = 0.0; 

	  for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	    data.chimes_current_rates->photoion_euv_heat_rate[i] = 0.0; 
	}

      for (i = 0; i < chimes_table_photoion_auger_fuv.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->photoion_auger_fuv_rate_coefficient[i] = 0.0; 

      for (i = 0; i < chimes_table_photoion_auger_euv.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->photoion_auger_euv_rate_coefficient[i] = 0.0; 

      for (i = 0; i < chimes_table_photodissoc_group1.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->photodissoc_group1_rate_coefficient[i] = 0.0; 

      for (i = 0; i < chimes_table_photodissoc_group2.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->photodissoc_group2_rate_coefficient[i] = 0.0; 
      
      for (i = 0; i < chimes_table_CO_photodissoc.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->CO_photodissoc_rate_coefficient[i] = 0.0; 

      // Sum over all spectra 
      for (j = 0; j < myGlobalVars->N_spectra; j++) 
	{
	  flux = myGasVars->isotropic_photon_density[j] * LIGHTSPEED; 
	  G0 = flux * myGasVars->G0_parameter[j]; 

	  for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++)
	      data.chimes_current_rates->photoion_fuv_rate_coefficient[i] += flux * data.chimes_current_rates->photoion_fuv_shield_factor[i] * chimes_table_photoion_fuv.sigmaPhot[i][j]; 
	    
	  for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	    data.chimes_current_rates->photoion_euv_rate_coefficient[i] += flux * data.chimes_current_rates->photoion_euv_shield_factor[i][j] * chimes_table_photoion_euv.sigmaPhot[i][j]; 

	  if (myGasVars->ThermEvolOn) 
	    {
	      for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++)
		data.chimes_current_rates->photoion_fuv_heat_rate[i] += flux * data.chimes_current_rates->photoion_fuv_shield_factor[i] * chimes_table_photoion_fuv.sigmaPhot[i][j] * chimes_table_photoion_fuv.epsilonPhot[i][j]; 

	      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++)
		data.chimes_current_rates->photoion_euv_heat_rate[i] += flux * data.chimes_current_rates->photoion_euv_shield_factor[i][j] * chimes_table_photoion_euv.sigmaPhot[i][j] * data.chimes_current_rates->photoion_euv_epsilon[i][j]; 
	    }

	  for (i = 0; i < chimes_table_photoion_auger_fuv.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photoion_auger_fuv_rate_coefficient[i] += flux * data.chimes_current_rates->photoion_fuv_shield_factor[chimes_table_photoion_auger_fuv.base_reaction[i]] * chimes_table_photoion_auger_fuv.sigmaPhot[i][j]; 

	  for (i = 0; i < chimes_table_photoion_auger_euv.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photoion_auger_euv_rate_coefficient[i] += flux * data.chimes_current_rates->photoion_euv_shield_factor[chimes_table_photoion_auger_euv.base_reaction[i]][j] * chimes_table_photoion_auger_euv.sigmaPhot[i][j]; 

	  for (i = 0; i < chimes_table_photodissoc_group1.N_reactions[data.mol_flag_index]; i++)
	    data.chimes_current_rates->photodissoc_group1_rate_coefficient[i] += G0 * data.chimes_current_rates->photodissoc_group1_shield_factor[i] * chimes_table_photodissoc_group1.rates[i]; 

	  if (data.mol_flag_index == 1) 
	    {
	      for (i = 0; i < chimes_table_photodissoc_group2.N_reactions[data.mol_flag_index]; i++)
		data.chimes_current_rates->photodissoc_group2_rate_coefficient[i] += G0 * data.chimes_current_rates->photodissoc_group2_shield_factor * chimes_table_photodissoc_group2.rates[i]; 

	      for (i = 0; i < chimes_table_CO_photodissoc.N_reactions[data.mol_flag_index]; i++)
		data.chimes_current_rates->CO_photodissoc_rate_coefficient[i] += G0 * data.chimes_current_rates->CO_photodissoc_shield_factor[i] * chimes_table_CO_photodissoc.rates[i]; 
	    }
	}
    }

  return; 
}


/* This routine updates the rate coefficients that 
 * can vary during the course of the integration. */ 
void update_rate_coefficients(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data, int mode) 
{
  int i, j, T_index, Psi_index, T_dust_index, NH2_index, b_index; 
  ChimesFloat log_T, log_Psi, log_T_dust, dT, dPsi, dT_dust, log_NH2, dNH2, log_b, db; 
  ChimesFloat flux, G0, n_over_cr; 

  // Determine table indices for interpolation 
  log_T = (ChimesFloat) log10(myGasVars->temperature); 
  chimes_get_table_index(chimes_table_bins.Temperatures, chimes_table_bins.N_Temperatures, log_T, &T_index, &dT);

  if (myGlobalVars->N_spectra > 0) 
    {
      G0 = 0.0; 
      for (i = 0; i < myGlobalVars->N_spectra; i++) 
	G0 += myGasVars->isotropic_photon_density[i] * LIGHTSPEED * myGasVars->G0_parameter[i]; 

      // In the following, we protect against division by zero, and 
      // against taking log(0). 
      log_Psi = (ChimesFloat) log10(chimes_max(G0 * exp(-data.extinction * G0_GAMMA) * pow(myGasVars->temperature, 0.5) / chimes_max(myGasVars->nH_tot * myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]], 1.0e-100), 1.0e-100));
    }
  else 
    log_Psi = chimes_table_bins.Psi[0]; 

  chimes_get_table_index(chimes_table_bins.Psi, chimes_table_bins.N_Psi, log_Psi, &Psi_index, &dPsi);
  
  if (mode) 
    {
      /* These reactions are only updated if mode == 1. 
       * This is the case if this routine has been been 
       * called from set_initial_rate_coefficients(), or 
       * if ThermEvolOn == 1 in the RhsFn. */
  
      // T_dependent reactions 
      for (i = 0; i < chimes_table_T_dependent.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->T_dependent_rate_coefficient[i] = pow(10.0, chimes_interpol_1d(chimes_table_T_dependent.rates[i], T_index, dT)); 

      // recombination_AB reactions 
      for (i = 0; i < chimes_table_recombination_AB.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->recombination_AB_rate_coefficient[i] = pow(10.0, chimes_interpol_1d(chimes_table_recombination_AB.rates[i][data.case_AB_index[i]], T_index, dT)); 

      // grain_recombination reactions 
      for (i = 0; i < chimes_table_grain_recombination.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->grain_recombination_rate_coefficient[i] = pow(10.0, chimes_interpol_2d(chimes_table_grain_recombination.rates[i], T_index, Psi_index, dT, dPsi)); 
  
      // H2_dust_formation 
      if (data.mol_flag_index == 1) 
	{
	  log_T_dust = (ChimesFloat) log10(myGlobalVars->grain_temperature); 
	  chimes_get_table_index(chimes_table_bins.Dust_Temperatures, chimes_table_bins.N_Dust_Temperatures, log_T_dust, &T_dust_index, &dT_dust); 

	  data.chimes_current_rates->H2_dust_formation_rate_coefficient = pow(10.0, chimes_interpol_2d(chimes_table_H2_dust_formation.rates, T_index, T_dust_index, dT, dT_dust)); 
	}

      if ((data.mol_flag_index == 1) && (myGlobalVars->N_spectra > 0)) 
	{
	  /* H2_photodissoc 
	   * The shielding factors depend on temperature, so 
	   * need to be updated when the temperature evolves. 
	   * However, the column densities are not updated. */
	  
	  if (myGlobalVars->cellSelfShieldingOn > 0) 
	    {
	      log_NH2 = (ChimesFloat) log10(chimes_max(data.H2_column, 1.0e-100)); 
	      log_b = (ChimesFloat) log10(chimes_max(myGasVars->doppler_broad, 1.0e-100)); 
	      chimes_get_table_index(chimes_table_bins.H2self_column_densities, chimes_table_bins.N_H2self_column_densities, log_NH2, &NH2_index, &dNH2); 
	      chimes_get_table_index(chimes_table_bins.b_turbulence, chimes_table_bins.N_b_turbulence, log_b, &b_index, &db); 
	      
	      for (i = 0; i < chimes_table_H2_photodissoc.N_reactions[data.mol_flag_index]; i++) 
		{
		  data.chimes_current_rates->H2_photodissoc_shield_factor[i] = pow(10.0, chimes_interpol_3d(chimes_table_H2_photodissoc.self_shielding[i], T_index, NH2_index, b_index, dT, dNH2, db)); 
		  data.chimes_current_rates->H2_photodissoc_shield_factor[i] *= exp(- chimes_table_H2_photodissoc.gamma[i] * data.extinction); 
		}
	    }
	  else 
	    {
	      for (i = 0; i < chimes_table_H2_photodissoc.N_reactions[data.mol_flag_index]; i++) 
		data.chimes_current_rates->H2_photodissoc_shield_factor[i] = 1.0; 
	    }
	  
	  // Zero the rate coefficients 
	  for (i = 0; i < chimes_table_H2_photodissoc.N_reactions[data.mol_flag_index]; i++) 
	    data.chimes_current_rates->H2_photodissoc_rate_coefficient[i] = 0.0; 
 
	  // Sum over all spectra 
	  for (j = 0; j < myGlobalVars->N_spectra; j++) 
	    {
	      flux = myGasVars->isotropic_photon_density[j] * LIGHTSPEED; 
	      for (i = 0; i < chimes_table_H2_photodissoc.N_reactions[1]; i++) 
		data.chimes_current_rates->H2_photodissoc_rate_coefficient[i] += flux * myGasVars->H2_dissocJ[j] * chimes_table_H2_photodissoc.rates[i] * data.chimes_current_rates->H2_photodissoc_shield_factor[i]; 
	    }
	}
    }

  if (data.mol_flag_index) 
    {
      // These groups only contain molecular reactions. 

      // H2 collis_dissoc. 
      if (mode) 
	{
	  /* k0, kLTE and the critical densities 
	   * are just functions of T, so they 
	   * only need to be interpolated when 
	   * mode == 1, i.e. when the initial rate 
	   * coefficients are set or if ThermEvolOn == 1. */ 
	  data.chimes_current_rates->H2_collis_dissoc_crit_H = pow(10.0, chimes_interpol_1d(chimes_table_H2_collis_dissoc.critical_density_H, T_index, dT)); 
	  data.chimes_current_rates->H2_collis_dissoc_crit_H2 = pow(10.0, chimes_interpol_1d(chimes_table_H2_collis_dissoc.critical_density_H2, T_index, dT)); 
	  data.chimes_current_rates->H2_collis_dissoc_crit_He = pow(10.0, chimes_interpol_1d(chimes_table_H2_collis_dissoc.critical_density_He, T_index, dT)); 

	  for (i = 0; i < chimes_table_H2_collis_dissoc.N_reactions[data.mol_flag_index]; i++) 
	    {
	      data.chimes_current_rates->H2_collis_dissoc_log_k0[i] = chimes_interpol_1d(chimes_table_H2_collis_dissoc.k0[i], T_index, dT); 
	      data.chimes_current_rates->H2_collis_dissoc_log_kLTE[i] = chimes_interpol_1d(chimes_table_H2_collis_dissoc.kLTE[i], T_index, dT); 
	    }
	}
      
      /* The rate coefficient itself also depends 
       * on the abundances of HI, H2 and HeI */ 
      n_over_cr = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] / data.chimes_current_rates->H2_collis_dissoc_crit_H; 
      n_over_cr += 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] / data.chimes_current_rates->H2_collis_dissoc_crit_H2; 
      n_over_cr += myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI]] / data.chimes_current_rates->H2_collis_dissoc_crit_He; 
      n_over_cr *= myGasVars->nH_tot; 

      for (i = 0; i < chimes_table_H2_collis_dissoc.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->H2_collis_dissoc_rate_coefficient[i] = pow(10.0, (((n_over_cr / (1.0 + n_over_cr)) * data.chimes_current_rates->H2_collis_dissoc_log_kLTE[i]) + ((1.0 / (1.0 + n_over_cr)) * data.chimes_current_rates->H2_collis_dissoc_log_k0[i]))); 


      // CO_cosmic_ray 
      for (i = 0; i < chimes_table_CO_cosmic_ray.N_reactions[data.mol_flag_index]; i++)
	data.chimes_current_rates->CO_cosmic_ray_rate_coefficient[i] = pow(10.0, chimes_interpol_1d(chimes_table_CO_cosmic_ray.rates[i], T_index, dT)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] * myGasVars->cr_rate;
    }
  
  return; 
} 

/* This routine updates the current rate, i.e. dx_i / dt, 
 * for each reaction. This typically takes the rate coefficient 
 * and multiplies it by the densities of each reactant, although 
 * the exact format depends on the reaction group. */ 
void update_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data) 
{ 
  int i, j, xHII_index; 
  ChimesFloat log_xHII, d_xHII; 

  // T_dependent reactions 
  for (i = 0; i < chimes_table_T_dependent.N_reactions[data.mol_flag_index]; i++) 
    {
      data.chimes_current_rates->T_dependent_rate[i] = data.chimes_current_rates->T_dependent_rate_coefficient[i]; 
      for (j = 0; j < 3; j++) 
	{
	  if (chimes_table_T_dependent.reactants[i][j] < 0) 
	    break; 
	  else 
	    data.chimes_current_rates->T_dependent_rate[i] *= myGasVars->nH_tot * myGasVars->abundances[chimes_table_T_dependent.reactants[i][j]]; 
	}
      
      /* This group contains 2-body and 3-body reactions, with rate coefficients 
       * in cm^3/s and cm^6/s, respectively. To get the current rate in units 
       * of s^-1 (i.e. dx_i/dt), we therefore need one and two factors of nH_tot, 
       * respectively. We therefore multiply by nH_tot for each reactant, and then 
       * divide out one factor of nH_tot at the end. */ 
      data.chimes_current_rates->T_dependent_rate[i] /= myGasVars->nH_tot; 
    }


  // constant reactions 
  /* All reactions in this group are 2-body, with rate coefficients in cm^3/s, 
   * so current_rate = rate_coefficient * x_0 * x_1 * nH_tot (in s^-1). */
  for (i = 0; i < chimes_table_constant.N_reactions[data.mol_flag_index]; i++) 
    data.chimes_current_rates->constant_rate[i] = chimes_table_constant.rates[i] * myGasVars->abundances[chimes_table_constant.reactants[i][0]] * myGasVars->abundances[chimes_table_constant.reactants[i][1]] * myGasVars->nH_tot; 


  // recombination_AB reactions 
  /* All reactions in this group are 2-body, with rate coefficients in cm^3/s, 
   * so current_rate = rate_coefficient * x_0 * x_1 * nH_tot (in s^-1). */
  for (i = 0; i < chimes_table_recombination_AB.N_reactions[data.mol_flag_index]; i++) 
    data.chimes_current_rates->recombination_AB_rate[i] = data.chimes_current_rates->recombination_AB_rate_coefficient[i] * myGasVars->abundances[chimes_table_recombination_AB.reactants[i][0]] * myGasVars->abundances[chimes_table_recombination_AB.reactants[i][1]] * myGasVars->nH_tot; 

  // grain_recombination reactions 
  // NOTE: Rate scales with dust_ratio 
  /* All reactions in this group are 2-body, with rate coefficients in cm^3/s, 
   * so current_rate = rate_coefficient * x_0 * x_1 * nH_tot (in s^-1). */
  for (i = 0; i < chimes_table_grain_recombination.N_reactions[data.mol_flag_index]; i++) 
    data.chimes_current_rates->grain_recombination_rate[i] = data.chimes_current_rates->grain_recombination_rate_coefficient[i] * myGasVars->abundances[chimes_table_grain_recombination.reactants[i][0]] * myGasVars->abundances[chimes_table_grain_recombination.reactants[i][1]] * myGasVars->nH_tot * myGasVars->dust_ratio; 


  // cosmic_ray reactions 
  if (myGasVars->cr_rate > 0.0) 
    {
      for (i = 0; i < chimes_table_cosmic_ray.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->cosmic_ray_rate[i] = myGasVars->cr_rate * chimes_table_cosmic_ray.rates[i] * myGasVars->abundances[chimes_table_cosmic_ray.reactants[i]]; 
      
      // secondary cosmic ray ionisation 
      log_xHII = (ChimesFloat) log10(chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]], 1.0e-100)); 
      chimes_get_table_index(chimes_table_bins.secondary_cosmic_ray_xHII, chimes_table_bins.N_secondary_cosmic_ray_xHII, log_xHII, &xHII_index, &d_xHII); 
      for (i = 0; i < 2; i++) 
	data.chimes_current_rates->cosmic_ray_rate[chimes_table_cosmic_ray.secondary_base_reaction[i]] *= 1.0 + pow(10.0, chimes_interpol_1d(chimes_table_cosmic_ray.secondary_ratio[i], xHII_index, d_xHII)); 

    }


  if (data.mol_flag_index) 
    {
      // H2_dust_formation 
      // NOTE: Rate scales with dust_ratio 
      data.chimes_current_rates->H2_dust_formation_rate = data.chimes_current_rates->H2_dust_formation_rate_coefficient * myGasVars->abundances[chimes_table_H2_dust_formation.reactants[0]] * myGasVars->dust_ratio * myGasVars->nH_tot; 

      // H2_collis_dissoc 
      for (i = 0; i < chimes_table_H2_collis_dissoc.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->H2_collis_dissoc_rate[i] = data.chimes_current_rates->H2_collis_dissoc_rate_coefficient[i] * myGasVars->abundances[chimes_table_H2_collis_dissoc.reactants[i][0]] * myGasVars->abundances[chimes_table_H2_collis_dissoc.reactants[i][1]] * myGasVars->nH_tot; 

      // CO_cosmic_ray
      // NOTE: the rate_coefficient contains a factor 1/sqrt(xCO). I have taken this
      // factor out of the rate_coefficient and into here, which gets multiplied by
      // the factor xCO from the reactant, hence the sqrt() when we multiply by the
      // reactant here. 
      for (i = 0; i < chimes_table_CO_cosmic_ray.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->CO_cosmic_ray_rate[i] = data.chimes_current_rates->CO_cosmic_ray_rate_coefficient[i] * sqrt(chimes_max(myGasVars->abundances[chimes_table_CO_cosmic_ray.reactants[i]], 0.0));
    }

  if (myGlobalVars->N_spectra > 0) 
    {
      // photoion_fuv 
      for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->photoion_fuv_rate[i] = data.chimes_current_rates->photoion_fuv_rate_coefficient[i] * myGasVars->abundances[chimes_table_photoion_fuv.reactants[i]]; 

      // photoion_euv 
      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->photoion_euv_rate[i] = data.chimes_current_rates->photoion_euv_rate_coefficient[i] * myGasVars->abundances[chimes_table_photoion_euv.reactants[i]]; 

      // photoion_auger_fuv 
      for (i = 0; i < chimes_table_photoion_auger_fuv.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->photoion_auger_fuv_rate[i] = data.chimes_current_rates->photoion_auger_fuv_rate_coefficient[i] * myGasVars->abundances[chimes_table_photoion_auger_fuv.reactants[i]]; 

      // photoion_auger_euv 
      for (i = 0; i < chimes_table_photoion_auger_euv.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->photoion_auger_euv_rate[i] = data.chimes_current_rates->photoion_auger_euv_rate_coefficient[i] * myGasVars->abundances[chimes_table_photoion_auger_euv.reactants[i]]; 


      // photodissoc_group1 
      for (i = 0; i < chimes_table_photodissoc_group1.N_reactions[data.mol_flag_index]; i++) 
	data.chimes_current_rates->photodissoc_group1_rate[i] = data.chimes_current_rates->photodissoc_group1_rate_coefficient[i] * myGasVars->abundances[chimes_table_photodissoc_group1.reactants[i]]; 

      if (data.mol_flag_index == 1) 
	{
	  // photodissoc_group2 
	  for (i = 0; i < chimes_table_photodissoc_group2.N_reactions[data.mol_flag_index]; i++) 
	    data.chimes_current_rates->photodissoc_group2_rate[i] = data.chimes_current_rates->photodissoc_group2_rate_coefficient[i] * myGasVars->abundances[chimes_table_photodissoc_group2.reactants[i]]; 

	  // H2_photodissoc 
	  for (i = 0; i < chimes_table_H2_photodissoc.N_reactions[data.mol_flag_index]; i++) 
	    data.chimes_current_rates->H2_photodissoc_rate[i] = data.chimes_current_rates->H2_photodissoc_rate_coefficient[i] * myGasVars->abundances[chimes_table_H2_photodissoc.reactants[i]]; 

	  // CO photodissoc 
	  for (i = 0; i < chimes_table_CO_photodissoc.N_reactions[data.mol_flag_index]; i++) 
	    data.chimes_current_rates->CO_photodissoc_rate[i] = data.chimes_current_rates->CO_photodissoc_rate_coefficient[i] * myGasVars->abundances[chimes_table_CO_photodissoc.reactants[i]]; 
	}
    }

  return; 
}

/* This routine loops through all reactions 
 * and updates the creation and destruction 
 * rates stored in the species structures. */
void update_rate_vector(struct Species_Structure *mySpecies, struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data) 
{
  int i, j; 

  /* NOTE: some reactions may involve species where the 
   * element abundance is zero, but the element_incl flag 
   * is one. These reactions are still computed below, 
   * and these species are included in the abundance array
   * but not the rate vector used within CVODE. However, 
   * such reactions will have a zero rate, because at 
   * least one reactant will have zero abundance, so it 
   * does not matter that we have not excised these reactions 
   * from the network, as they are adding zero to the 
   * total rates. */ 

  // T_dependent reactions 
  for (i = 0; i < chimes_table_T_dependent.N_reactions[data.mol_flag_index]; i++) 
    {
      for (j = 0; j < 3; j++) 
	{
	  if (chimes_table_T_dependent.reactants[i][j] < 0) 
	    break; 
	  else 
	    mySpecies[chimes_table_T_dependent.reactants[i][j]].destruction_rate += data.chimes_current_rates->T_dependent_rate[i]; 
	} 
      for (j = 0; j < 3; j++) 
	{
	  if (chimes_table_T_dependent.products[i][j] < 0) 
	    break; 
	  else 
	    mySpecies[chimes_table_T_dependent.products[i][j]].creation_rate += data.chimes_current_rates->T_dependent_rate[i]; 
	} 
    }

  // constant reactions 
  for (i = 0; i < chimes_table_constant.N_reactions[data.mol_flag_index]; i++) 
    {
      /* Only contains 2-body reactions, so 
       * there are always two reactants. */ 
      for (j = 0; j < 2; j++) 
	mySpecies[chimes_table_constant.reactants[i][j]].destruction_rate += data.chimes_current_rates->constant_rate[i]; 

      for (j = 0; j < 3; j++) 
	{
	  if (chimes_table_constant.products[i][j] < 0) 
	    break; 
	  else 
	    mySpecies[chimes_table_constant.products[i][j]].creation_rate += data.chimes_current_rates->constant_rate[i]; 
	} 
    }


  // recombination_AB reactions 
  for (i = 0; i < chimes_table_recombination_AB.N_reactions[data.mol_flag_index]; i++) 
    {
      /* All reactions in this group have two 
       * reactants and one product. */ 
      for (j = 0; j < 2; j++) 
	mySpecies[chimes_table_recombination_AB.reactants[i][j]].destruction_rate += data.chimes_current_rates->recombination_AB_rate[i]; 

      mySpecies[chimes_table_recombination_AB.products[i]].creation_rate += data.chimes_current_rates->recombination_AB_rate[i]; 
    }


  // grain_recombination reactions 
  for (i = 0; i < chimes_table_grain_recombination.N_reactions[data.mol_flag_index]; i++) 
    {
      /* All reactions in this group have two 
       * reactants and one product. */ 
      for (j = 0; j < 2; j++) 
	mySpecies[chimes_table_grain_recombination.reactants[i][j]].destruction_rate += data.chimes_current_rates->grain_recombination_rate[i]; 

      mySpecies[chimes_table_grain_recombination.products[i]].creation_rate += data.chimes_current_rates->grain_recombination_rate[i]; 
    }  


  // cosmic_ray reactions 
  if (myGasVars->cr_rate > 0.0) 
    {
      for (i = 0; i < chimes_table_cosmic_ray.N_reactions[data.mol_flag_index]; i++) 
	{
	  // All reactions have one reactant
	  mySpecies[chimes_table_cosmic_ray.reactants[i]].destruction_rate += data.chimes_current_rates->cosmic_ray_rate[i]; 

	  for (j = 0; j < 3; j++) 
	    {
	      if (chimes_table_cosmic_ray.products[i][j] < 0) 
		break; 
	      else 
		mySpecies[chimes_table_cosmic_ray.products[i][j]].creation_rate += data.chimes_current_rates->cosmic_ray_rate[i]; 
	    } 
	}
    }


  if (data.mol_flag_index) 
    {
      // H2_dust_formation 
      for (j = 0; j < 2; j++) 
	mySpecies[chimes_table_H2_dust_formation.reactants[j]].destruction_rate += data.chimes_current_rates->H2_dust_formation_rate; 
      
      mySpecies[chimes_table_H2_dust_formation.products[0]].creation_rate += data.chimes_current_rates->H2_dust_formation_rate; 

      // H2_collis_dissoc 
      for (i = 0; i < chimes_table_H2_collis_dissoc.N_reactions[data.mol_flag_index]; i++) 
	{
	  for (j = 0; j < 2; j++) 
	    mySpecies[chimes_table_H2_collis_dissoc.reactants[i][j]].destruction_rate += data.chimes_current_rates->H2_collis_dissoc_rate[i]; 

	  for (j = 0; j < 3; j++) 
	    mySpecies[chimes_table_H2_collis_dissoc.products[i][j]].creation_rate += data.chimes_current_rates->H2_collis_dissoc_rate[i]; 
	}

      // CO_cosmic_ray reactions 
      if (myGasVars->cr_rate > 0.0) 
	{
	  for (i = 0; i < chimes_table_CO_cosmic_ray.N_reactions[data.mol_flag_index]; i++) 
	    {
	      // All reactions have one reactant
	      mySpecies[chimes_table_CO_cosmic_ray.reactants[i]].destruction_rate += data.chimes_current_rates->CO_cosmic_ray_rate[i]; 
	      for (j = 0; j < 2; j++) 
		mySpecies[chimes_table_CO_cosmic_ray.products[i][j]].creation_rate += data.chimes_current_rates->CO_cosmic_ray_rate[i];
	    }
	}
    }

  if (myGlobalVars->N_spectra > 0) 
    {
      // photoion_fuv 
      for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++) 
	{
	  mySpecies[chimes_table_photoion_fuv.reactants[i]].destruction_rate += data.chimes_current_rates->photoion_fuv_rate[i]; 
	  mySpecies[chimes_table_photoion_fuv.products[i][0]].creation_rate += data.chimes_current_rates->photoion_fuv_rate[i]; 
	  mySpecies[chimes_table_photoion_fuv.products[i][1]].creation_rate += data.chimes_current_rates->photoion_fuv_rate[i]; 
	}

      // photoion_euv 
      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	{
	  mySpecies[chimes_table_photoion_euv.reactants[i]].destruction_rate += data.chimes_current_rates->photoion_euv_rate[i]; 
	  mySpecies[chimes_table_photoion_euv.products[i][0]].creation_rate += data.chimes_current_rates->photoion_euv_rate[i]; 
	  mySpecies[chimes_table_photoion_euv.products[i][1]].creation_rate += data.chimes_current_rates->photoion_euv_rate[i]; 
	}

      // photoion_auger_fuv 
      for (i = 0; i < chimes_table_photoion_auger_fuv.N_reactions[data.mol_flag_index]; i++) 
	{
	  mySpecies[chimes_table_photoion_auger_fuv.reactants[i]].destruction_rate += data.chimes_current_rates->photoion_auger_fuv_rate[i]; 
	  mySpecies[chimes_table_photoion_auger_fuv.products[i][0]].creation_rate += data.chimes_current_rates->photoion_auger_fuv_rate[i]; 
	  mySpecies[chimes_table_photoion_auger_fuv.products[i][1]].creation_rate += data.chimes_current_rates->photoion_auger_fuv_rate[i] * chimes_table_photoion_auger_fuv.number_of_electrons[i]; 
	}

      // photoion_auger_euv 
      for (i = 0; i < chimes_table_photoion_auger_euv.N_reactions[data.mol_flag_index]; i++) 
	{
	  mySpecies[chimes_table_photoion_auger_euv.reactants[i]].destruction_rate += data.chimes_current_rates->photoion_auger_euv_rate[i]; 
	  mySpecies[chimes_table_photoion_auger_euv.products[i][0]].creation_rate += data.chimes_current_rates->photoion_auger_euv_rate[i]; 
	  mySpecies[chimes_table_photoion_auger_euv.products[i][1]].creation_rate += data.chimes_current_rates->photoion_auger_euv_rate[i] * chimes_table_photoion_auger_euv.number_of_electrons[i]; 
	}

      // photodissoc_group1 
      for (i = 0; i < chimes_table_photodissoc_group1.N_reactions[data.mol_flag_index]; i++) 
	{
	  mySpecies[chimes_table_photodissoc_group1.reactants[i]].destruction_rate += data.chimes_current_rates->photodissoc_group1_rate[i]; 
	  mySpecies[chimes_table_photodissoc_group1.products[i][0]].creation_rate += data.chimes_current_rates->photodissoc_group1_rate[i]; 
	  mySpecies[chimes_table_photodissoc_group1.products[i][1]].creation_rate += data.chimes_current_rates->photodissoc_group1_rate[i]; 
	}
      
      
      if (data.mol_flag_index == 1) 
	{
	  // photodissoc_group2 
	  for (i = 0; i < chimes_table_photodissoc_group2.N_reactions[data.mol_flag_index]; i++) 
	    {
	      mySpecies[chimes_table_photodissoc_group2.reactants[i]].destruction_rate += data.chimes_current_rates->photodissoc_group2_rate[i]; 
	      mySpecies[chimes_table_photodissoc_group2.products[i][0]].creation_rate += data.chimes_current_rates->photodissoc_group2_rate[i]; 
	      mySpecies[chimes_table_photodissoc_group2.products[i][1]].creation_rate += data.chimes_current_rates->photodissoc_group2_rate[i]; 
	    }

	  // H2_photodissoc 
	  for (i = 0; i < chimes_table_H2_photodissoc.N_reactions[data.mol_flag_index]; i++) 
	    {
	      mySpecies[chimes_table_H2_photodissoc.reactants[i]].destruction_rate += data.chimes_current_rates->H2_photodissoc_rate[i]; 
	      mySpecies[chimes_table_H2_photodissoc.products[i][0]].creation_rate += data.chimes_current_rates->H2_photodissoc_rate[i]; 
	      mySpecies[chimes_table_H2_photodissoc.products[i][1]].creation_rate += data.chimes_current_rates->H2_photodissoc_rate[i]; 
	    }

	  // CO_photodissoc 
	  for (i = 0; i < chimes_table_CO_photodissoc.N_reactions[data.mol_flag_index]; i++) 
	    {
	      mySpecies[chimes_table_CO_photodissoc.reactants[i]].destruction_rate += data.chimes_current_rates->CO_photodissoc_rate[i]; 
	      mySpecies[chimes_table_CO_photodissoc.products[i][0]].creation_rate += data.chimes_current_rates->CO_photodissoc_rate[i]; 
	      mySpecies[chimes_table_CO_photodissoc.products[i][1]].creation_rate += data.chimes_current_rates->CO_photodissoc_rate[i]; 
	    }
	}
    }

  return; 
}


  /* This function determines which species are to be 
   * included in the CVODE integration. Here we exclude 
   * species that contain elements whose abundances are 
   * zero, even if the element is included in the overall
   * network. */
void set_species_structures(struct Species_Structure *mySpecies, struct gasVariables *myGasVars, int *total_network, int *nonmolecular_network, struct globalVariables *myGlobalVars)
{
  int i;
  int inclSpeciesFlags[10];
  int inclSpeciesFlag_CO;

  for (i = 0; i < 10; i++)
    {
      if (myGasVars->element_abundances[i] > METALS_MINIMUM_THRESHOLD)
	inclSpeciesFlags[i] = 1;
      else
	inclSpeciesFlags[i] = 0;
    }

  if (myGasVars->element_abundances[1] > METALS_MINIMUM_THRESHOLD && myGasVars->element_abundances[3] > METALS_MINIMUM_THRESHOLD)
    inclSpeciesFlag_CO = 1;
  else
    inclSpeciesFlag_CO = 0;
 	
  /* Determine which species are included */ 
  for (i = myGlobalVars->speciesIndices[sp_elec]; i <= myGlobalVars->speciesIndices[sp_Hm]; i++)
    {
      mySpecies[i].include_species = 1;
      mySpecies[i].element_abundance = 1.0;
    }

  for (i = myGlobalVars->speciesIndices[sp_HeI]; i <= myGlobalVars->speciesIndices[sp_HeIII]; i++)
    {
      mySpecies[i].include_species = inclSpeciesFlags[0];
      mySpecies[i].element_abundance = myGasVars->element_abundances[0];
    }

  if (myGlobalVars->element_included[0] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_CI]; i <= myGlobalVars->speciesIndices[sp_Cm]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[1];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[1];
	}
    }

  if (myGlobalVars->element_included[1] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_NI]; i <= myGlobalVars->speciesIndices[sp_NVIII]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[2];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[2];
	}
    }

  if (myGlobalVars->element_included[2] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_OI]; i <= myGlobalVars->speciesIndices[sp_Om]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[3];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[3];
	}
    }

  if (myGlobalVars->element_included[3] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_NeI]; i <= myGlobalVars->speciesIndices[sp_NeXI]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[4];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[4];
	}
    }

  if (myGlobalVars->element_included[4] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_MgI]; i <= myGlobalVars->speciesIndices[sp_MgXIII]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[5];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[5];
	}
    }

  if (myGlobalVars->element_included[5] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_SiI]; i <= myGlobalVars->speciesIndices[sp_SiXV]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[6];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[6];
	}
    }

  if (myGlobalVars->element_included[6] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_SI]; i <= myGlobalVars->speciesIndices[sp_SXVII]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[7];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[7];
	}
    }

  if (myGlobalVars->element_included[7] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_CaI]; i <= myGlobalVars->speciesIndices[sp_CaXXI]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[8];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[8];
	}
    }

  if (myGlobalVars->element_included[8] == 1)
    {
      for (i = myGlobalVars->speciesIndices[sp_FeI]; i <= myGlobalVars->speciesIndices[sp_FeXXVII]; i++)
	{
	  mySpecies[i].include_species = inclSpeciesFlags[9];
	  mySpecies[i].element_abundance = myGasVars->element_abundances[9];
	}
    }

  mySpecies[myGlobalVars->speciesIndices[sp_H2]].include_species = 1;
  mySpecies[myGlobalVars->speciesIndices[sp_H2p]].include_species = 1;
  mySpecies[myGlobalVars->speciesIndices[sp_H3p]].include_species = 1;

  mySpecies[myGlobalVars->speciesIndices[sp_H2]].element_abundance = 1.0;
  mySpecies[myGlobalVars->speciesIndices[sp_H2p]].element_abundance = 1.0;
  mySpecies[myGlobalVars->speciesIndices[sp_H3p]].element_abundance = 1.0;

  if (myGlobalVars->element_included[2] == 1)
    {
      mySpecies[myGlobalVars->speciesIndices[sp_OH]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[sp_H2O]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[sp_O2]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[sp_OHp]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[sp_H2Op]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[sp_H3Op]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[sp_O2p]].include_species = inclSpeciesFlags[3];

      mySpecies[myGlobalVars->speciesIndices[sp_OH]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[sp_H2O]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[sp_O2]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[sp_OHp]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[sp_H2Op]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[sp_H3Op]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[sp_O2p]].element_abundance = myGasVars->element_abundances[3];
    }

  if (myGlobalVars->element_included[0] == 1)
    {
      mySpecies[myGlobalVars->speciesIndices[sp_C2]].include_species = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH]].include_species = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH2]].include_species = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH3p]].include_species = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CHp]].include_species = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH2p]].include_species = inclSpeciesFlags[1];

      mySpecies[myGlobalVars->speciesIndices[sp_C2]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH2]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH3p]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CHp]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CH2p]].element_abundance = myGasVars->element_abundances[1];
    }

  if (myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    {
      mySpecies[myGlobalVars->speciesIndices[sp_HCOp]].include_species = inclSpeciesFlag_CO;
      mySpecies[myGlobalVars->speciesIndices[sp_CO]].include_species = inclSpeciesFlag_CO;
      mySpecies[myGlobalVars->speciesIndices[sp_COp]].include_species = inclSpeciesFlag_CO;
      mySpecies[myGlobalVars->speciesIndices[sp_HOCp]].include_species = inclSpeciesFlag_CO;

      mySpecies[myGlobalVars->speciesIndices[sp_HCOp]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_CO]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_COp]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[sp_HOCp]].element_abundance = myGasVars->element_abundances[1];
    }

  /* Now loop through this array and determine the total
   * number of species that are included in the network. */
  *total_network = 0;
  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    {
      if (mySpecies[i].include_species == 1)
	*total_network += 1;
    }

  /* Now subtract from this the number of 
   * molecules included in the network. */
  *nonmolecular_network = *total_network;
  for (i = sp_H2; i <= sp_O2p; i++)
    {
      if (myGlobalVars->speciesIndices[i] > -1)
	{
	  if (mySpecies[myGlobalVars->speciesIndices[i]].include_species == 1)
	    *nonmolecular_network -= 1;
	}
    }
}

void redshift_dependent_UVB_copy_lowz_to_hiz(struct globalVariables *myGlobalVars) 
{
  int i, j, k, l; 

  for (i = 0; i < chimes_table_photoion_fuv.N_reactions[1]; i++) 
    {
      chimes_table_redshift_dependent_UVB.photoion_fuv_sigmaPhot[i][1] = chimes_table_redshift_dependent_UVB.photoion_fuv_sigmaPhot[i][0]; 
      chimes_table_redshift_dependent_UVB.photoion_fuv_epsilonPhot[i][1] = chimes_table_redshift_dependent_UVB.photoion_fuv_epsilonPhot[i][0]; 
    }

  for (i = 0; i < chimes_table_photoion_euv.N_reactions[1]; i++) 
    {
      chimes_table_redshift_dependent_UVB.photoion_euv_sigmaPhot[i][1] = chimes_table_redshift_dependent_UVB.photoion_euv_sigmaPhot[i][0]; 

      for (j = 0; j < 3; j++) 
	{
	  for (k = 0; k < chimes_table_bins.N_Column_densities; k++) 
	    chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_1D[i][1][j][k] = chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_1D[i][0][j][k];
	}

      for (j = 0; j < 6; j++) 
	{
	  for (k = 0; k < chimes_table_bins.N_Column_densities; k++) 
	    {
	      for (l = 0; l < chimes_table_bins.N_Column_densities; l++) 
		chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_2D[i][1][j][k][l] = chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_2D[i][0][j][k][l];
	    }
	}
    }

  for (i = 0; i < chimes_table_photoion_auger_fuv.N_reactions[1]; i++) 
    chimes_table_redshift_dependent_UVB.photoion_auger_fuv_sigmaPhot[i][1] = chimes_table_redshift_dependent_UVB.photoion_auger_fuv_sigmaPhot[i][0]; 

  for (i = 0; i < chimes_table_photoion_auger_euv.N_reactions[1]; i++) 
    chimes_table_redshift_dependent_UVB.photoion_auger_euv_sigmaPhot[i][1] = chimes_table_redshift_dependent_UVB.photoion_auger_euv_sigmaPhot[i][0]; 

  chimes_table_redshift_dependent_UVB.isotropic_photon_density[1] = chimes_table_redshift_dependent_UVB.isotropic_photon_density[0]; 
  chimes_table_redshift_dependent_UVB.G0_parameter[1] = chimes_table_redshift_dependent_UVB.G0_parameter[0]; 
  chimes_table_redshift_dependent_UVB.H2_dissocJ[1] = chimes_table_redshift_dependent_UVB.H2_dissocJ[0]; 

  if (myGlobalVars->use_redshift_dependent_eqm_tables == 1) 
    {
      for (i = 0; i < chimes_table_eqm_abundances.N_Temperatures; i++) 
	{
	  for (j = 0; j < chimes_table_eqm_abundances.N_Densities; j++) 
	    {
	      for (k = 0; k < chimes_table_eqm_abundances.N_Metallicities; k++) 
		{
		  for (l = 0; l < myGlobalVars->totalNumberOfSpecies; l++) 
		    chimes_table_redshift_dependent_UVB.eqm_abundances[1].Abundances[l][i][j][k] = chimes_table_redshift_dependent_UVB.eqm_abundances[0].Abundances[l][i][j][k]; 
		}
	    }
	}
    }
}  

void interpolate_redshift_dependent_UVB(struct globalVariables *myGlobalVars) 
{
  ChimesFloat low_z, hi_z, dz, dz_m; 
  ChimesFloat redshift = myGlobalVars->redshift; 
  int N_reactions_all, i, j, k, l; 
  int z_index_low = chimes_table_redshift_dependent_UVB.z_index_low; 
  int z_index_hi = chimes_table_redshift_dependent_UVB.z_index_hi; 
  int Nz = chimes_table_redshift_dependent_UVB.N_redshifts; 
  int spectrum_index = myGlobalVars->redshift_dependent_UVB_index; 
  int first_UVB_load_flag = 0; 

  /* First, determine whether we need to load 
   * a new UVB table. */ 
  if ((z_index_low < 0) || (z_index_hi < 0)) 
    {
      // No UVB tables have yet been read in
      first_UVB_load_flag = 1; 

      if (redshift >= chimes_table_redshift_dependent_UVB.redshift_bins[Nz - 1]) 
	{
	  // Redshift is higher than the highest bin 
	  low_z = chimes_table_redshift_dependent_UVB.redshift_bins[Nz - 1]; 
	  hi_z = low_z; 
	  load_redshift_dependent_UVB(low_z, 0, myGlobalVars); 
	  redshift_dependent_UVB_copy_lowz_to_hiz(myGlobalVars); 

	  z_index_low = Nz - 1; 
	  z_index_hi = Nz - 1; 
	  chimes_table_redshift_dependent_UVB.z_index_low = z_index_low; 
	  chimes_table_redshift_dependent_UVB.z_index_hi = z_index_hi; 
	}
      else if (redshift <= chimes_table_redshift_dependent_UVB.redshift_bins[0]) 
	{
	  // Redshift is lower than the lowest bin 
	  low_z = chimes_table_redshift_dependent_UVB.redshift_bins[0]; 
	  hi_z = low_z; 
	  load_redshift_dependent_UVB(low_z, 0, myGlobalVars); 
	  redshift_dependent_UVB_copy_lowz_to_hiz(myGlobalVars); 

	  z_index_low = 0; 
	  z_index_hi = 0; 
	  chimes_table_redshift_dependent_UVB.z_index_low = z_index_low; 
	  chimes_table_redshift_dependent_UVB.z_index_hi = z_index_hi; 
	}
      else 
	{
	  z_index_hi = 0; 
	  while (chimes_table_redshift_dependent_UVB.redshift_bins[z_index_hi] <= redshift) 
	    z_index_hi += 1; 
	  
	  z_index_low = z_index_hi - 1; 

	  low_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_low]; 
	  hi_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_hi]; 
	  
	  load_redshift_dependent_UVB(low_z, 0, myGlobalVars); 
	  load_redshift_dependent_UVB(hi_z, 1, myGlobalVars); 

	  chimes_table_redshift_dependent_UVB.z_index_low = z_index_low; 
	  chimes_table_redshift_dependent_UVB.z_index_hi = z_index_hi; 
	}
    }
  else 
    {
      low_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_low]; 
      hi_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_hi]; 
      
      if (redshift < low_z)
	{
	  /* Current redshift has moved 
	   * below low_z. */ 
	  if ((z_index_low == 0) && (z_index_hi > z_index_low)) 
	    {
	      /* Current redshift has moved 
	       * below lowest redshift bin. */ 
	      z_index_hi = z_index_low; 
	      hi_z = low_z; 
	      chimes_table_redshift_dependent_UVB.z_index_hi = z_index_hi; 
	      
	      redshift_dependent_UVB_copy_lowz_to_hiz(myGlobalVars); 
	    }
	  else if (z_index_low > 0) 
	    {
	      // Determine new z_index_low
	      z_index_low = 0; 
	      while (chimes_table_redshift_dependent_UVB.redshift_bins[z_index_low] <= redshift) 
		z_index_low += 1; 
	      z_index_low -= 1; 

	      if (z_index_hi - z_index_low == 2) 
		{
		  /* We have only moved down 
		   * a single redshift bin. */ 
		  z_index_hi = z_index_low + 1; 
		  low_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_low]; 
		  hi_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_hi]; 

		  redshift_dependent_UVB_copy_lowz_to_hiz(myGlobalVars); 
		  load_redshift_dependent_UVB(low_z, 0, myGlobalVars); 
		}
	      else 
		{
		  /* We have moved down 
		   * multiple redshift bins. */ 
		  z_index_hi = z_index_low + 1; 
		  low_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_low]; 
		  hi_z = chimes_table_redshift_dependent_UVB.redshift_bins[z_index_hi]; 
		  
		  load_redshift_dependent_UVB(low_z, 0, myGlobalVars); 
		  load_redshift_dependent_UVB(hi_z, 1, myGlobalVars); 
		}
	    }
	}
    }

  // Interpolate tables to current redshift
  if (z_index_low == z_index_hi) 
    {
      dz = 0.0; 
      dz_m = 1.0; 
    }
  else 
    {
      dz = (redshift - low_z) / (hi_z - low_z); 
      dz_m = 1.0 - dz; 
    }

  if (!((redshift > myGlobalVars->reionisation_redshift) && (first_UVB_load_flag == 0))) 
    {
      // photoion_fuv 
      N_reactions_all = chimes_table_photoion_fuv.N_reactions[1]; 
      for (i = 0; i < N_reactions_all; i++) 
	{
	  chimes_table_photoion_fuv.sigmaPhot[i][spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.photoion_fuv_sigmaPhot[i][0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.photoion_fuv_sigmaPhot[i][1]) * dz)); 
	  chimes_table_photoion_fuv.epsilonPhot[i][spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.photoion_fuv_epsilonPhot[i][0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.photoion_fuv_epsilonPhot[i][1]) * dz)); 
	}

      // photoion_euv 
      N_reactions_all = chimes_table_photoion_euv.N_reactions[1]; 
      for (i = 0; i < N_reactions_all; i++) 
	{
	  chimes_table_photoion_euv.sigmaPhot[i][spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.photoion_euv_sigmaPhot[i][0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.photoion_euv_sigmaPhot[i][1]) * dz)); 

	  for (j = 0; j < 3; j++) 
	    {
	      for (k = 0; k < chimes_table_bins.N_Column_densities; k++) 
		chimes_table_photoion_euv.shieldFactor_1D[i][spectrum_index][j][k] = (chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_1D[i][0][j][k] * dz_m) + (chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_1D[i][1][j][k] * dz); 
	    }

	  for (j = 0; j < 6; j++) 
	    {
	      for (k = 0; k < chimes_table_bins.N_Column_densities; k++) 
		{
		  for (l = 0; l < chimes_table_bins.N_Column_densities; l++) 
		    chimes_table_photoion_euv.shieldFactor_2D[i][spectrum_index][j][k][l] = (chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_2D[i][0][j][k][l] * dz_m) + (chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_2D[i][1][j][k][l] * dz); 
		}
	    }
	}

      // photoion_auger_fuv 
      N_reactions_all = chimes_table_photoion_auger_fuv.N_reactions[1]; 
      for (i = 0; i < N_reactions_all; i++) 
	chimes_table_photoion_auger_fuv.sigmaPhot[i][spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.photoion_auger_fuv_sigmaPhot[i][0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.photoion_auger_fuv_sigmaPhot[i][1]) * dz)); 

      // photoion_auger_euv 
      N_reactions_all = chimes_table_photoion_auger_euv.N_reactions[1]; 
      for (i = 0; i < N_reactions_all; i++) 
	chimes_table_photoion_auger_euv.sigmaPhot[i][spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.photoion_auger_euv_sigmaPhot[i][0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.photoion_auger_euv_sigmaPhot[i][1]) * dz)); 

      if (redshift > myGlobalVars->reionisation_redshift)
	chimes_table_spectra.isotropic_photon_density[spectrum_index] = 0.0; 
      else 
	chimes_table_spectra.isotropic_photon_density[spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.isotropic_photon_density[0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.isotropic_photon_density[1]) * dz)); 

      chimes_table_spectra.G0_parameter[spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.G0_parameter[0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.G0_parameter[1]) * dz)); 
      chimes_table_spectra.H2_dissocJ[spectrum_index] = pow(10.0, (log10(chimes_table_redshift_dependent_UVB.H2_dissocJ[0]) * dz_m) + (log10(chimes_table_redshift_dependent_UVB.H2_dissocJ[1]) * dz)); 


      if (myGlobalVars->use_redshift_dependent_eqm_tables == 1) 
	{
	  for (i = 0; i < chimes_table_eqm_abundances.N_Temperatures; i++) 
	    {
	      for (j = 0; j < chimes_table_eqm_abundances.N_Densities; j++) 
		{
		  for (k = 0; k < chimes_table_eqm_abundances.N_Metallicities; k++) 
		    {
		      for (l = 0; l < myGlobalVars->totalNumberOfSpecies; l++) 
			chimes_table_eqm_abundances.Abundances[l][i][j][k] = (chimes_table_redshift_dependent_UVB.eqm_abundances[0].Abundances[l][i][j][k] * dz_m) + (chimes_table_redshift_dependent_UVB.eqm_abundances[1].Abundances[l][i][j][k] * dz); 
		    }
		}
	    }
	}
    }
}
