/****************************************************************************
 * This file is part of CHIMES.
 * Copyright (c) 2020 Alexander Richings (alexander.j.richings@durham.ac.uk)
 *
 * CHIMES is free software: you can redistribute it and/or modify
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
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "chimes_proto.h"
#include "chimes_vars.h"

/**
 * @brief Calculates the total number density.
 *
 * Sums the number densities of all ions and molecules
 * in the network to get the total number density.
 *
 * @param my_abundances Abundance array.
 * @param nH Total hydrogen number density.
 * @param myGlobalVars The #globalVariables struct.
 */
ChimesFloat calculate_total_number_density(
    ChimesFloat *my_abundances, ChimesFloat nH,
    struct globalVariables *myGlobalVars) {
  int i;
  ChimesFloat result = 0.0;

  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    result += my_abundances[i] * nH;

  return result;
}

/**
 * @brief Calculates the mean molecular weight.
 *
 * Calculates the mean molecular weight from the abundance
 * array. This is defined as:
 * mu = (1.0 + sum(element_abundances[i] * Z[i])) / sum(abundances[i]).
 * Where element_abundances[] are the abundances of He and the metals
 * relative to H, Z[] are the atomic masses of each element, and
 * abundances[] are the abundances of the individual ions and molecules
 * relative to H.
 *
 * @param myGasVars The gasVariables struct.
 * @param myGlobalVars The globalVariables struct.
 */
ChimesFloat calculate_mean_molecular_weight(
    struct gasVariables *myGasVars, struct globalVariables *myGlobalVars) {
  ChimesFloat denominator = 0.0;
  int i;

  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    denominator += myGasVars->abundances[i];

  return (1.0 + myGasVars->element_abundances[0] * 4.0 +
          myGasVars->element_abundances[1] * 12.0 +
          myGasVars->element_abundances[2] * 14.0 +
          myGasVars->element_abundances[3] * 16.0 +
          myGasVars->element_abundances[4] * 20.0 +
          myGasVars->element_abundances[5] * 24.0 +
          myGasVars->element_abundances[6] * 28.0 +
          myGasVars->element_abundances[7] * 32.0 +
          myGasVars->element_abundances[8] * 40.0 +
          myGasVars->element_abundances[9] * 56.0) /
         denominator;
}

/**
 * @brief Compton cooling rate.
 *
 * Calculates the Compton cooling rate from the
 * CMB radiation. It returns Lambda/nH^2, in
 * units of erg cm^3 s^-1
 *
 * @param T Gas temperature.
 * @param Tcmb CMB temperatures.
 * @param xe Electron abundance.
 * @param nH Total hydrogen number density.
 */
ChimesFloat compton_cooling(ChimesFloat T, ChimesFloat Tcmb, ChimesFloat xe,
                            ChimesFloat nH) {
  return 1.017e-37 * pow(Tcmb, 4) * (T - Tcmb) * xe / nH; /* Lambda/nH^2 */
}

/**
 * @brief OH rotational cooling.
 *
 * Calculates the cooling rate from rotational
 * transitions of the OH molecule. This is calculated
 * from Hollenbach and McKee 1979, ApJS, 41, 555 (see
 * their equations 6.21 and 6.22, and the parameters
 * in their Table 3).
 * This function returns Lambda/nH^2, in units of
 * erg cm^3 s^-1.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 * @param data The #UserData struct.
 */
ChimesFloat OH_rotational_cooling(struct gasVariables *myGasVars,
                                  struct globalVariables *myGlobalVars,
                                  struct UserData data) {
  ChimesFloat dv, N_tau, tau_T, c_tau, n_cr, ym;

  // Thermal velocity dispersion in cgs
  dv = sqrt(
      3.0 * BOLTZMANNCGS * myGasVars->temperature /
      (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars)));

  N_tau = 1.485e11 * dv;
  tau_T = 4.0 * data.OH_column /
          (10.0 * (myGasVars->temperature / 27.0) * 6.8e-4 * N_tau);

  c_tau = tau_T *
          pow(2 * PI * log(2.13 + pow(tau_T / 2.718281828459, 2.0)), 0.5) /
          ((exp(-data.extinction) / (1.0 + pow(data.extinction, 2.0))) +
           2.0 * data.extinction * pow(log(1 + (tau_T / 2.718281828459)), 0.5) *
               pow(log(tau_T / (data.extinction * 2.718281828459)), 0.5));

  if (isnan(c_tau) != 0) c_tau = 0.0;

  n_cr = 1.5e10 * pow(myGasVars->temperature / 1.0e3, 0.5);
  ym = log(1.0 + (c_tau / (1.0 + 10.0 * (n_cr / myGasVars->nH_tot))));

  // Lambda / nH^2 (erg cm^3 s^-1)
  return myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] *
         (2.0 * pow(BOLTZMANNCGS * myGasVars->temperature, 2.0) * 2.3e-2 /
          (myGasVars->nH_tot * 27.0 * BOLTZMANNCGS)) *
         ((2.0 + ym + 0.6 * pow(ym, 2.0)) /
          (1.0 + c_tau + (n_cr / myGasVars->nH_tot) +
           1.5 * pow(n_cr / myGasVars->nH_tot, 0.5)));
}

/**
 * @brief Update the cooling rates.
 *
 * Calculates the rates of the various cooling
 * and heating channels and stores them in the
 * #chimes_current_rates_struct struct within
 * the #UserData struct.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 * @param data The #UserData struct.
 */
void update_cooling_rates(struct gasVariables *myGasVars,
                          struct globalVariables *myGlobalVars,
                          struct UserData data) {
  int i, T_index, nHI_index, ne_index, nHII_index;
  ChimesFloat dT, dnHI, dne, dnHII;
  ChimesFloat log_T, log_nHI, log_ne, log_nHII;

  log_T = (ChimesFloat)log10(myGasVars->temperature);
  chimes_get_table_index(chimes_table_bins.Temperatures,
                         chimes_table_bins.N_Temperatures, log_T, &T_index,
                         &dT);

  for (i = 0; i < chimes_table_cooling.N_coolants; i++)
    data.chimes_current_rates->cooling_rate[i] = pow(
        10.0, chimes_interpol_1d(chimes_table_cooling.rates[i], T_index, dT));

  if (chimes_table_cooling.N_coolants_2d > 0) {
    if (log_T <
        chimes_table_bins
            .cool_2d_Temperatures[chimes_table_bins.N_cool_2d_Temperatures -
                                  1]) {
      chimes_get_table_index(chimes_table_bins.cool_2d_Temperatures,
                             chimes_table_bins.N_cool_2d_Temperatures, log_T,
                             &T_index, &dT);

      log_ne = (ChimesFloat)log10(chimes_max(
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] *
              myGasVars->nH_tot,
          1.0e-100));
      chimes_get_table_index(chimes_table_bins.cool_2d_ElectronDensities,
                             chimes_table_bins.N_cool_2d_ElectronDensities,
                             log_ne, &ne_index, &dne);

      for (i = 0; i < chimes_table_cooling.N_coolants_2d; i++)
        data.chimes_current_rates->cooling_rate_2d[i] =
            pow(10.0, chimes_interpol_2d(chimes_table_cooling.rates_2d[i],
                                         T_index, ne_index, dT, dne));
    } else {
      chimes_get_table_index(chimes_table_bins.cool_hiT_2d_Temperatures,
                             chimes_table_bins.N_cool_hiT_2d_Temperatures,
                             log_T, &T_index, &dT);

      for (i = 0; i < chimes_table_cooling.N_coolants_2d; i++)
        data.chimes_current_rates->cooling_rate_2d[i] =
            pow(10.0, chimes_interpol_1d(chimes_table_cooling.rates_hiT_2d[i],
                                         T_index, dT));
    }
  }

  if (chimes_table_cooling.N_coolants_4d > 0) {
    if (log_T <
        chimes_table_bins
            .cool_4d_Temperatures[chimes_table_bins.N_cool_4d_Temperatures -
                                  1]) {
      chimes_get_table_index(chimes_table_bins.cool_4d_Temperatures,
                             chimes_table_bins.N_cool_4d_Temperatures, log_T,
                             &T_index, &dT);

      log_nHI = (ChimesFloat)log10(chimes_max(
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] *
              myGasVars->nH_tot,
          1.0e-100));
      log_ne = (ChimesFloat)log10(chimes_max(
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] *
              myGasVars->nH_tot,
          1.0e-100));
      log_nHII = (ChimesFloat)log10(chimes_max(
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]] *
              myGasVars->nH_tot,
          1.0e-100));

      chimes_get_table_index(chimes_table_bins.cool_4d_HIDensities,
                             chimes_table_bins.N_cool_4d_HIDensities, log_nHI,
                             &nHI_index, &dnHI);
      chimes_get_table_index(chimes_table_bins.cool_4d_ElectronDensities,
                             chimes_table_bins.N_cool_4d_ElectronDensities,
                             log_ne, &ne_index, &dne);
      chimes_get_table_index(chimes_table_bins.cool_4d_HIIDensities,
                             chimes_table_bins.N_cool_4d_HIIDensities, log_nHII,
                             &nHII_index, &dnHII);

      for (i = 0; i < chimes_table_cooling.N_coolants_4d; i++)
        data.chimes_current_rates->cooling_rate_4d[i] =
            pow(10.0, chimes_interpol_4d(chimes_table_cooling.rates_4d[i],
                                         T_index, nHI_index, ne_index,
                                         nHII_index, dT, dnHI, dne, dnHII));
    } else {
      chimes_get_table_index(chimes_table_bins.cool_hiT_4d_Temperatures,
                             chimes_table_bins.N_cool_hiT_4d_Temperatures,
                             log_T, &T_index, &dT);

      for (i = 0; i < chimes_table_cooling.N_coolants_4d; i++)
        data.chimes_current_rates->cooling_rate_4d[i] =
            pow(10.0, chimes_interpol_1d(chimes_table_cooling.rates_hiT_4d[i],
                                         T_index, dT));
    }
  }

  return;
}

/**
 * @brief Calculates the total cooling rate.
 *
 * Sums all cooling and heating channels to calculate
 * the total cooling rate. This function can return
 * either the net cooling rate, cooling only, or heating
 * only, depending on the mode argument. The rate is
 * returned in units of erg cm^-3 s^-1.
 * If the net cooling rate is returned, a positive value
 * indicates net cooling, while a negative value indicates
 * net heating. If only the cooling or heating rate is
 * returned, the value will always be positive.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 * @param data The #UserData struct.
 * @param mode An integer flag. 0: return net cooling; 1: return cooling only;
 * 2: return heating only.
 */
ChimesFloat calculate_total_cooling_rate(struct gasVariables *myGasVars,
                                         struct globalVariables *myGlobalVars,
                                         struct UserData data, int mode) {
  int i, xHII_index, Psi_index, T_index, T_mol_index, N_CO_rot_index,
      N_CO_vib_index;
  int N_H2O_rot_index, N_H2O_vib_index, T_H2O_index;
  ChimesFloat x_elec, log_xHII, d_xHII, log_Psi, dPsi, log_T, dT, G0;
  ChimesFloat H2_lowDens, H2_LTE, dT_mol, xHI, xH2, H2_crit_density;
  ChimesFloat dN_CO_rot, dN_CO_vib, log_N_eff;
  ChimesFloat CO_rot_L0, CO_rot_Llte, CO_rot_nhalf, CO_rot_a, CO_rot_neff;
  ChimesFloat CO_vib_L0, CO_vib_Llte, CO_vib_neff, H2O_rot_neff, H2O_vib_neff;
  ChimesFloat dN_H2O_rot, dN_H2O_vib, dT_H2O;
  ChimesFloat H2O_rot_L0, H2O_rot_Llte, H2O_rot_nhalf, H2O_rot_a, H2O_vib_L0,
      H2O_vib_Llte;
  ChimesFloat cr_secondary;
  ChimesFloat total_cooling = 0.0;
  ChimesFloat total_heating = 0.0;

  update_cooling_rates(myGasVars, myGlobalVars, data);

  x_elec = myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]];

  for (i = 0; i < chimes_table_cooling.N_coolants; i++)
    total_cooling += data.chimes_current_rates->cooling_rate[i] *
                     myGasVars->abundances[chimes_table_cooling.coolants[i]] *
                     x_elec;

  for (i = 0; i < chimes_table_cooling.N_coolants_2d; i++)
    total_cooling +=
        data.chimes_current_rates->cooling_rate_2d[i] *
        myGasVars->abundances[chimes_table_cooling.coolants_2d[i]] * x_elec;

  for (i = 0; i < chimes_table_cooling.N_coolants_4d; i++)
    total_cooling +=
        data.chimes_current_rates->cooling_rate_4d[i] *
        myGasVars->abundances[chimes_table_cooling.coolants_4d[i]] /
        myGasVars->nH_tot;

  // Photoheating
  if (myGlobalVars->N_spectra > 0) {
    for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index];
         i++)
      total_heating +=
          data.chimes_current_rates->photoion_fuv_heat_rate[i] *
          myGasVars->abundances[chimes_table_photoion_fuv.reactants[i]] /
          myGasVars->nH_tot;

    for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index];
         i++)
      total_heating +=
          data.chimes_current_rates->photoion_euv_heat_rate[i] *
          myGasVars->abundances[chimes_table_photoion_euv.reactants[i]] /
          myGasVars->nH_tot;
  }

  // Cosmic ray heating
  for (i = 0; i < chimes_table_cosmic_ray.N_reactions[data.mol_flag_index]; i++)
    total_heating += 3.2e-11 * data.chimes_current_rates->cosmic_ray_rate[i] /
                     myGasVars->nH_tot;

  // Correct for secondary cosmic rays
  log_xHII = (ChimesFloat)log10(chimes_max(
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]], 1.0e-100));
  chimes_get_table_index(chimes_table_bins.secondary_cosmic_ray_xHII,
                         chimes_table_bins.N_secondary_cosmic_ray_xHII,
                         log_xHII, &xHII_index, &d_xHII);
  for (i = 0; i < 2; i++) {
    cr_secondary =
        pow(10.0, chimes_interpol_1d(chimes_table_cosmic_ray.secondary_ratio[i],
                                     xHII_index, d_xHII));
    total_heating -= 3.2e-11 *
                     data.chimes_current_rates->cosmic_ray_rate
                         [chimes_table_cosmic_ray.secondary_base_reaction[i]] *
                     cr_secondary / (myGasVars->nH_tot * (1.0 + cr_secondary));
  }

  // Compton cooling from the CMB
  total_cooling += compton_cooling(
      myGasVars->temperature, myGlobalVars->cmb_temperature,
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]],
      myGasVars->nH_tot);

  if (data.mol_flag_index == 1) {
    // H2 rovibrational cooling
    log_T = (ChimesFloat)log10(myGasVars->temperature);
    chimes_get_table_index(chimes_table_bins.mol_cool_Temperatures,
                           chimes_table_bins.N_mol_cool_Temperatures, log_T,
                           &T_mol_index, &dT_mol);
    H2_lowDens =
        pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_H2,
                                     T_mol_index, dT_mol)) *
        myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]];
    H2_lowDens +=
        pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_HI,
                                     T_mol_index, dT_mol)) *
        myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]];
    H2_lowDens +=
        pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_HII,
                                     T_mol_index, dT_mol)) *
        myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]];
    H2_lowDens +=
        pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_HeI,
                                     T_mol_index, dT_mol)) *
        myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI]];
    H2_lowDens +=
        pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_elec,
                                     T_mol_index, dT_mol)) *
        myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]];
    H2_lowDens *= myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]];

    if (H2_lowDens > 0.0) {
      H2_LTE = pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2_cool_LTE,
                                            T_mol_index, dT_mol)) *
               myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] /
               myGasVars->nH_tot;
      total_cooling += H2_LTE / (1.0 + (H2_LTE / H2_lowDens));
    }

    // H2 collis dissoc heating
    xHI = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]];
    xH2 = myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]];

    total_cooling +=
        7.2e-12 *
        data.chimes_current_rates->H2_collis_dissoc_rate_coefficient
            [chimes_table_H2_collis_dissoc.Heating_reaction_index] *
        xHI * xH2;
    total_cooling +=
        7.2e-12 *
        data.chimes_current_rates->T_dependent_rate_coefficient
            [chimes_table_T_dependent.H2_collis_dissoc_heating_reaction_index] *
        pow(xH2, 2.0);

    // Gas-phase H2 formation
    if (xHI + xH2 == 0)
      H2_crit_density = 0.0;
    else
      H2_crit_density =
          (xHI + xH2) /
          ((xHI / data.chimes_current_rates->H2_collis_dissoc_crit_H) +
           (xH2 / data.chimes_current_rates->H2_collis_dissoc_crit_H2));

    total_heating +=
        (((2.93e-12 *
           data.chimes_current_rates->T_dependent_rate
               [chimes_table_T_dependent.H2_form_heating_reaction_index]) +
          (5.65e-12 *
           data.chimes_current_rates->constant_rate
               [chimes_table_constant.H2_form_heating_reaction_index])) *
         (1.0 / (myGasVars->nH_tot + H2_crit_density)));

    // Dust-catalysed H2 formation
    total_heating +=
        7.16e-12 *
        (data.chimes_current_rates->H2_dust_formation_rate /
         myGasVars->nH_tot) *
        (myGasVars->nH_tot / (myGasVars->nH_tot + H2_crit_density));

    // CO cooling
    if ((myGlobalVars->element_included[0] == 1) &&
        (myGlobalVars->element_included[2] == 1)) {
      // N_eff units: cm^-2 per km s^-1
      if (myGlobalVars->StaticMolCooling == 1)
        log_N_eff = (ChimesFloat)log10(
            chimes_max(1.0e5 * data.CO_column /
                           (sqrt(3.0 * BOLTZMANNCGS * myGasVars->temperature /
                                 (PROTON_MASS * calculate_mean_molecular_weight(
                                                    myGasVars, myGlobalVars)))),
                       1.0e-100));
      else
        log_N_eff = (ChimesFloat)log10(chimes_max(
            1.0e5 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] *
                myGasVars->nH_tot /
                chimes_max(fabs(myGasVars->divVel), 1.0e-100),
            1.0e-100));

      chimes_get_table_index(chimes_table_bins.CO_cool_rot_ColumnDensities,
                             chimes_table_bins.N_CO_cool_rot_ColumnDensities,
                             log_N_eff, &N_CO_rot_index, &dN_CO_rot);
      chimes_get_table_index(chimes_table_bins.CO_cool_vib_ColumnDensities,
                             chimes_table_bins.N_CO_cool_vib_ColumnDensities,
                             log_N_eff, &N_CO_vib_index, &dN_CO_vib);

      CO_rot_L0 =
          pow(10.0, chimes_interpol_1d(chimes_table_cooling.CO_cool_rot_L0,
                                       T_mol_index, dT_mol));
      CO_rot_Llte =
          pow(10.0, chimes_interpol_2d(chimes_table_cooling.CO_cool_rot_Llte,
                                       T_mol_index, N_CO_rot_index, dT_mol,
                                       dN_CO_rot));
      CO_rot_nhalf =
          pow(10.0, chimes_interpol_2d(chimes_table_cooling.CO_cool_rot_nhalf,
                                       T_mol_index, N_CO_rot_index, dT_mol,
                                       dN_CO_rot));
      CO_rot_a = pow(10.0, chimes_interpol_2d(
                               chimes_table_cooling.CO_cool_rot_a, T_mol_index,
                               N_CO_rot_index, dT_mol, dN_CO_rot));

      CO_rot_neff =
          myGasVars->nH_tot *
          (xH2 + 9.857 * pow((myGasVars->temperature / 1.0e3), 0.25) * xHI +
           680.13 * pow(myGasVars->temperature, -0.25) * x_elec);

      total_cooling +=
          xH2 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] /
          ((1.0 / CO_rot_L0) + (CO_rot_neff / CO_rot_Llte) +
           (1.0 / CO_rot_L0) * pow((CO_rot_neff / CO_rot_nhalf), CO_rot_a) *
               (1.0 - (CO_rot_nhalf * CO_rot_L0 / CO_rot_Llte)));

      CO_vib_L0 =
          pow(10.0, chimes_interpol_1d(chimes_table_cooling.CO_cool_vib_L0,
                                       T_mol_index, dT_mol));
      CO_vib_Llte =
          pow(10.0, chimes_interpol_2d(chimes_table_cooling.CO_cool_vib_Llte,
                                       T_mol_index, N_CO_vib_index, dT_mol,
                                       dN_CO_vib));

      CO_vib_neff =
          myGasVars->nH_tot *
          (xH2 + 50.0 * xHI +
           9035.09 * exp(68.0 / pow(myGasVars->temperature, 1.0 / 3.0)) *
               pow(myGasVars->temperature / 300.0, 0.938) * x_elec);

      total_cooling +=
          xH2 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] /
          ((1.0 / CO_vib_L0) + (CO_vib_neff / CO_vib_Llte));
    }

    // H2O and OH cooling
    if (myGlobalVars->element_included[2] == 1) {
      if (myGlobalVars->StaticMolCooling == 1)
        log_N_eff = (ChimesFloat)log10(
            chimes_max(1.0e5 * data.H2O_column /
                           (sqrt(3.0 * BOLTZMANNCGS * myGasVars->temperature /
                                 (PROTON_MASS * calculate_mean_molecular_weight(
                                                    myGasVars, myGlobalVars)))),
                       1.0e-100));
      else
        log_N_eff = (ChimesFloat)log10(chimes_max(
            1.0e5 *
                myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] *
                myGasVars->nH_tot /
                chimes_max(fabs(myGasVars->divVel), 1.0e-100),
            1.0e-100));

      // rotational cooling
      chimes_get_table_index(chimes_table_bins.H2O_cool_rot_ColumnDensities,
                             chimes_table_bins.N_H2O_cool_rot_ColumnDensities,
                             log_N_eff, &N_H2O_rot_index, &dN_H2O_rot);

      H2O_rot_neff =
          myGasVars->nH_tot *
          (xH2 + 10.0 * xHI +
           pow(10.0,
               (-8.02 + (15.749 / pow(myGasVars->temperature, 1.0 / 6.0)) -
                (47.137 / pow(myGasVars->temperature, 1.0 / 3.0)) +
                (76.648 / pow(myGasVars->temperature, 0.5)) -
                (60.191 / pow(myGasVars->temperature, 2.0 / 3.0)))) *
               x_elec / (7.4e-12 * pow(myGasVars->temperature, 0.5)));
      if (log_T >= 2.0) {
        chimes_get_table_index(chimes_table_bins.H2O_cool_hiT_Temperatures,
                               chimes_table_bins.N_H2O_cool_hiT_Temperatures,
                               log_T, &T_H2O_index, &dT_H2O);

        H2O_rot_L0 = pow(
            10.0, chimes_interpol_1d(chimes_table_cooling.H2O_cool_rot_hiT_L0,
                                     T_H2O_index, dT_H2O));
        H2O_rot_Llte =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2O_cool_rot_hiT_Llte,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));
        H2O_rot_nhalf =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2O_cool_rot_hiT_nhalf,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));
        H2O_rot_a =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2O_cool_rot_hiT_a, T_H2O_index,
                          N_H2O_rot_index, dT_H2O, dN_H2O_rot));

        total_cooling +=
            myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * xH2 /
            ((1.0 / H2O_rot_L0) + (H2O_rot_neff / H2O_rot_Llte) +
             (1.0 / H2O_rot_L0) *
                 pow((H2O_rot_neff / H2O_rot_nhalf), H2O_rot_a) *
                 (1.0 - (H2O_rot_nhalf * H2O_rot_L0 / H2O_rot_Llte)));
      } else {
        chimes_get_table_index(chimes_table_bins.H2O_cool_lowT_Temperatures,
                               chimes_table_bins.N_H2O_cool_lowT_Temperatures,
                               log_T, &T_H2O_index, &dT_H2O);

        H2O_rot_L0 = pow(
            10.0,
            chimes_interpol_1d(chimes_table_cooling.H2Oortho_cool_rot_lowT_L0,
                               T_H2O_index, dT_H2O));
        H2O_rot_Llte =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2Oortho_cool_rot_lowT_Llte,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));
        H2O_rot_nhalf =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2Oortho_cool_rot_lowT_nhalf,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));
        H2O_rot_a =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2Oortho_cool_rot_lowT_a,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));

        total_cooling +=
            0.75 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] *
            xH2 /
            ((1.0 / H2O_rot_L0) + (H2O_rot_neff / H2O_rot_Llte) +
             (1.0 / H2O_rot_L0) *
                 pow((H2O_rot_neff / H2O_rot_nhalf), H2O_rot_a) *
                 (1.0 - (H2O_rot_nhalf * H2O_rot_L0 / H2O_rot_Llte)));

        H2O_rot_L0 = pow(
            10.0,
            chimes_interpol_1d(chimes_table_cooling.H2Opara_cool_rot_lowT_L0,
                               T_H2O_index, dT_H2O));
        H2O_rot_Llte =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2Opara_cool_rot_lowT_Llte,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));
        H2O_rot_nhalf =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2Opara_cool_rot_lowT_nhalf,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));
        H2O_rot_a =
            pow(10.0, chimes_interpol_2d(
                          chimes_table_cooling.H2Opara_cool_rot_lowT_a,
                          T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot));

        total_cooling +=
            0.25 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] *
            xH2 /
            ((1.0 / H2O_rot_L0) + (H2O_rot_neff / H2O_rot_Llte) +
             (1.0 / H2O_rot_L0) *
                 pow((H2O_rot_neff / H2O_rot_nhalf), H2O_rot_a) *
                 (1.0 - (H2O_rot_nhalf * H2O_rot_L0 / H2O_rot_Llte)));
      }

      // vibrational cooling
      chimes_get_table_index(chimes_table_bins.H2O_cool_vib_ColumnDensities,
                             chimes_table_bins.N_H2O_cool_vib_ColumnDensities,
                             log_N_eff, &N_H2O_vib_index, &dN_H2O_vib);

      H2O_vib_L0 =
          pow(10.0, chimes_interpol_1d(chimes_table_cooling.H2O_cool_vib_L0,
                                       T_H2O_index, dT_H2O));
      H2O_vib_Llte =
          pow(10.0, chimes_interpol_2d(chimes_table_cooling.H2O_cool_vib_Llte,
                                       T_H2O_index, N_H2O_vib_index, dT_H2O,
                                       dN_H2O_vib));

      H2O_vib_neff =
          myGasVars->nH_tot *
          (xH2 + 10.0 * xHI +
           4.0625e8 * exp(47.5 / pow(myGasVars->temperature, 1.0 / 3.0)) *
               pow(myGasVars->temperature, -0.5) * x_elec);

      total_cooling +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * xH2 /
          ((1.0 / H2O_vib_L0) + (H2O_vib_neff / H2O_vib_Llte));

      total_cooling += OH_rotational_cooling(myGasVars, myGlobalVars, data);
    }

    chimes_get_table_index(chimes_table_bins.Temperatures,
                           chimes_table_bins.N_Temperatures, log_T, &T_index,
                           &dT);

    if (myGlobalVars->N_spectra > 0) {
      // H2 photodissoc heating
      total_heating += 6.4e-13 *
                       data.chimes_current_rates->H2_photodissoc_rate[0] /
                       myGasVars->nH_tot;

      // H2 UV pumping
      total_heating += 2.7e-11 *
                       data.chimes_current_rates->H2_photodissoc_rate[0] *
                       (1.0 / (myGasVars->nH_tot + H2_crit_density));

      // Photoelectric dust heating & grain recombination cooling
      // Note that dust processes are only included when
      // the molecular network is switched on.
      G0 = 0.0;
      for (i = 0; i < myGlobalVars->N_spectra; i++)
        G0 += myGasVars->isotropic_photon_density[i] * LIGHTSPEED *
              myGasVars->G0_parameter[i];

      G0 *= exp(-data.extinction * G0_GAMMA);

      if ((G0 > 0.0) && (x_elec > 0.0)) {
        // Include a factor phi_pah = 0.5 in the denominator
        log_Psi = (ChimesFloat)log10(chimes_max(
            G0 * pow(myGasVars->temperature, 0.5) /
                chimes_max(myGasVars->nH_tot * x_elec * 0.5, 1.0e-100),
            1.0e-100));

        chimes_get_table_index(chimes_table_bins.Psi, chimes_table_bins.N_Psi,
                               log_Psi, &Psi_index, &dPsi);

        total_heating +=
            pow(10.0,
                chimes_interpol_2d(chimes_table_cooling.photoelectric_heating,
                                   T_index, Psi_index, dT, dPsi)) *
            G0 * myGasVars->dust_ratio / myGasVars->nH_tot;
        total_cooling +=
            pow(10.0,
                chimes_interpol_2d(chimes_table_cooling.grain_recombination,
                                   T_index, Psi_index, dT, dPsi)) *
            myGasVars->dust_ratio *
            myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]];
      }
    }

    // Gas-grain energy transfer
    total_cooling +=
        pow(10.0, chimes_interpol_1d(chimes_table_cooling.gas_grain_transfer,
                                     T_index, dT)) *
        myGasVars->dust_ratio *
        (myGasVars->temperature - myGlobalVars->grain_temperature);
  }

  if (myGlobalVars->hybrid_cooling_mode == 1) {
    if (myGlobalVars->hybrid_cooling_fn == NULL) {
      printf(
          "CHIMES error: hybrid_cooling_mode == %d, but the hybrid_cooling_fn "
          "has not been set. \n",
          myGlobalVars->hybrid_cooling_mode);
      exit(EXIT_FAILURE);
    }

    total_heating += (ChimesFloat)(*myGlobalVars->hybrid_cooling_fn)(
        myGasVars, myGlobalVars);
  }

  // Convert units to erg/cm3/s
  total_cooling *= pow(myGasVars->nH_tot, 2.0);
  total_heating *= pow(myGasVars->nH_tot, 2.0);
  total_heating += myGasVars->constant_heating_rate;

  if (mode == 0)
    return total_cooling - total_heating;
  else if (mode == 1)
    return total_cooling;
  else if (mode == 2)
    return total_heating;
  else {
    printf(
        "CHIMES error: calculate_total_cooling_rate() called with mode == %d. "
        "Allowed modes: 0 (net cooling), 1 (cooling only), 2 (heating only).\n",
        mode);
    exit(EXIT_FAILURE);
  }
}

/**
 * @brief Evolve the cooling for equilibrium abundances.
 *
 * If the abundances have been set to equilibrium, from the
 * pre-computed equilibrium tables, and thermal evolution
 * is switched on, then this routine is used to evolve
 * the temperature, using a bisection integration method.
 * As the temperature evolves, the equilibrium abundances
 * are re-computed from the tables, and are then used
 * to update the net cooling rate.
 *
 * @param data The #UserData struct.
 */
void do_equilibrium_cooling(struct UserData data) {
  int i, maxIter;
  ChimesFloat u, u_old, u_upper, u_lower, du, LambdaNet, dt;

  dt = data.myGasVars->hydro_timestep;

  data.myGasVars->temperature =
      chimes_max(data.myGasVars->temperature, data.myGasVars->TempFloor);

  set_equilibrium_abundances_from_tables(data);

  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                           data.myGasVars->ThermEvolOn);
  update_rates(data.myGasVars, data.myGlobalVars, data);

  LambdaNet =
      -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0);

  if (data.myGasVars->temperature <= data.myGasVars->TempFloor &&
      LambdaNet <= 0.0)
    return;

  u = data.myGasVars->temperature * 1.5 *
      calculate_total_number_density(data.myGasVars->abundances,
                                     data.myGasVars->nH_tot,
                                     data.myGlobalVars) *
      BOLTZMANNCGS;
  u_old = u;
  u_upper = u;
  u_lower = u;

  /* If the cooling rate is small, take explicit solution. */
  if (fabs(LambdaNet * dt) < 0.10 * u_old) {
    u = u_old + LambdaNet * dt;
    data.myGasVars->temperature =
        chimes_max(u / (1.5 *
                        calculate_total_number_density(
                            data.myGasVars->abundances, data.myGasVars->nH_tot,
                            data.myGlobalVars) *
                        BOLTZMANNCGS),
                   data.myGasVars->TempFloor);
    set_equilibrium_abundances_from_tables(data);

    // Check that explicit solution is valid
    update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                             data.myGasVars->ThermEvolOn);
    update_rates(data.myGasVars, data.myGlobalVars, data);

    LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars,
                                              data, 0);

    if (fabs(LambdaNet * dt) < 0.10 * u_old && fabs(LambdaNet * dt) < 0.10 * u)
      return;
    else {
      /* New cooling rate has increased, and explicit
       * solution is no longer valid. Reset and
       * continue with implicit solution. */
      u = u_old;
      data.myGasVars->temperature =
          chimes_max(u / (1.5 *
                          calculate_total_number_density(
                              data.myGasVars->abundances,
                              data.myGasVars->nH_tot, data.myGlobalVars) *
                          BOLTZMANNCGS),
                     data.myGasVars->TempFloor);
      set_equilibrium_abundances_from_tables(data);
    }
  }

  i = 0;
  maxIter = 150;

  /* Bracketing */
  if (u - u_old - LambdaNet * dt < 0.0) /* heating */
  {
    u_upper *= sqrt(1.2);
    u_lower /= sqrt(1.2);

    data.myGasVars->temperature =
        u_upper / (1.5 *
                   calculate_total_number_density(data.myGasVars->abundances,
                                                  data.myGasVars->nH_tot,
                                                  data.myGlobalVars) *
                   BOLTZMANNCGS);
    set_equilibrium_abundances_from_tables(data);

    update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                             data.myGasVars->ThermEvolOn);
    update_rates(data.myGasVars, data.myGlobalVars, data);

    LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars,
                                              data, 0);

    while (u_upper - u_old - LambdaNet * dt < 0.0 && i < maxIter) {
      u_upper *= 1.2;
      u_lower *= 1.2;

      data.myGasVars->temperature =
          u_upper / (1.5 *
                     calculate_total_number_density(data.myGasVars->abundances,
                                                    data.myGasVars->nH_tot,
                                                    data.myGlobalVars) *
                     BOLTZMANNCGS);
      set_equilibrium_abundances_from_tables(data);

      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                               data.myGasVars->ThermEvolOn);
      update_rates(data.myGasVars, data.myGlobalVars, data);

      LambdaNet = -calculate_total_cooling_rate(data.myGasVars,
                                                data.myGlobalVars, data, 0);

      i++;
    }

    if (i == maxIter)
      printf("WARNING: Problem with eqm cooling finding the upper bound.\n");
  } else /* cooling */
  {
    u_upper *= sqrt(1.2);
    u_lower /= sqrt(1.2);

    data.myGasVars->temperature =
        u_lower / (1.5 *
                   calculate_total_number_density(data.myGasVars->abundances,
                                                  data.myGasVars->nH_tot,
                                                  data.myGlobalVars) *
                   BOLTZMANNCGS);
    if (data.myGasVars->temperature <= data.myGasVars->TempFloor) {
      data.myGasVars->temperature = data.myGasVars->TempFloor;
      u_lower = data.myGasVars->TempFloor * 1.5 *
                calculate_total_number_density(data.myGasVars->abundances,
                                               data.myGasVars->nH_tot,
                                               data.myGlobalVars) *
                BOLTZMANNCGS;
      set_equilibrium_abundances_from_tables(data);

      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                               data.myGasVars->ThermEvolOn);
      update_rates(data.myGasVars, data.myGlobalVars, data);

      LambdaNet = -calculate_total_cooling_rate(data.myGasVars,
                                                data.myGlobalVars, data, 0);
    } else {
      set_equilibrium_abundances_from_tables(data);

      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                               data.myGasVars->ThermEvolOn);
      update_rates(data.myGasVars, data.myGlobalVars, data);

      LambdaNet = -calculate_total_cooling_rate(data.myGasVars,
                                                data.myGlobalVars, data, 0);

      while (u_lower - u_old - LambdaNet * dt > 0.0 && i < maxIter) {
        u_upper /= 1.2;
        u_lower /= 1.2;

        data.myGasVars->temperature =
            u_lower / (1.5 *
                       calculate_total_number_density(
                           data.myGasVars->abundances, data.myGasVars->nH_tot,
                           data.myGlobalVars) *
                       BOLTZMANNCGS);
        if (data.myGasVars->temperature <= data.myGasVars->TempFloor) {
          data.myGasVars->temperature = data.myGasVars->TempFloor;
          u_lower = data.myGasVars->TempFloor * 1.5 *
                    calculate_total_number_density(data.myGasVars->abundances,
                                                   data.myGasVars->nH_tot,
                                                   data.myGlobalVars) *
                    BOLTZMANNCGS;
          set_equilibrium_abundances_from_tables(data);

          update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                                   data.myGasVars->ThermEvolOn);
          update_rates(data.myGasVars, data.myGlobalVars, data);

          LambdaNet = -calculate_total_cooling_rate(data.myGasVars,
                                                    data.myGlobalVars, data, 0);
          break;
        }

        set_equilibrium_abundances_from_tables(data);

        update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                                 data.myGasVars->ThermEvolOn);
        update_rates(data.myGasVars, data.myGlobalVars, data);

        LambdaNet = -calculate_total_cooling_rate(data.myGasVars,
                                                  data.myGlobalVars, data, 0);

        i++;
      }
      if (i == maxIter)
        printf("WARNING: Problem with eqm cooling finding the lower bound.\n");
    }
    if (u_lower - u_old - LambdaNet * dt > 0.0 && i < maxIter)
      return; /* u_lower reached TempFloor, but is still above converged
                 solution. */
  }

  /* Iterate to convergence */
  i = 0;

  do {
    u = 0.5 * (u_lower + u_upper);

    data.myGasVars->temperature =
        u / (1.5 *
             calculate_total_number_density(data.myGasVars->abundances,
                                            data.myGasVars->nH_tot,
                                            data.myGlobalVars) *
             BOLTZMANNCGS);
    set_equilibrium_abundances_from_tables(data);

    update_rate_coefficients(data.myGasVars, data.myGlobalVars, data,
                             data.myGasVars->ThermEvolOn);
    update_rates(data.myGasVars, data.myGlobalVars, data);

    LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars,
                                              data, 0);

    if (u - u_old - LambdaNet * dt > 0.0)
      u_upper = u;
    else
      u_lower = u;

    du = u_upper - u_lower;
    i++;
  } while (fabs(du / u) > 1.0e-6 && i < maxIter);

  if (i >= maxIter) printf("WARNING: eqm cooling failed to converge.\n");

  return;
}
