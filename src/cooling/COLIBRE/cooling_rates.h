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

/* Config parameters. */
#ifndef SWIFT_COLIBRE_COOLING_RATES_H
#define SWIFT_COLIBRE_COOLING_RATES_H

#include "../config.h"

/* Local includes. */
#include "chemistry_struct.h"
#include "cooling_tables.h"
#include "error.h"
#include "exp10.h"
#include "interpolate.h"

__attribute__((always_inline)) INLINE int element_from_table_to_code(int i) {

#ifdef SWIFT_DEBUG_CHECKS
  if ((i >= colibre_cooling_N_elementtypes) || (i < 0))
    error("Outside range of elements in cooling tables");
#endif

  switch (i) {
    case element_H:
      return chemistry_element_H;
    case element_He:
      return chemistry_element_He;
    case element_C:
      return chemistry_element_C;
    case element_N:
      return chemistry_element_N;
    case element_O:
      return chemistry_element_O;
    case element_Ne:
      return chemistry_element_Ne;
    case element_Mg:
      return chemistry_element_Mg;
    case element_Si:
      return chemistry_element_Si;
      /* S and Ca are not tracked individually; their abundance is
       * assumed to be the same as Si (with respect to solar) */
    case element_S:
      return chemistry_element_Si;
    case element_Ca:
      return chemistry_element_Si;
    case element_Fe:
      return chemistry_element_Fe;
      /* other elements, if used, scale with metallicity */
    case element_OA:
      return -1;
  }

  return -1;
}

/**
 * @brief Compute ratio of mass fraction to solar mass fraction
 * for each element carried by a given particle.
 *
 * The solar abundances are taken from the tables themselves.
 *
 * The COLIBRE chemistry model does not track S and Ca. We assume
 * that their abundance with respect to solar is the same as
 * the ratio for Si.
 * We optionally apply a correction if the user asked for a different
 * ratio.
 *
 * We also re-order the elements such that they match the order of the
 * tables. This is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe].
 *
 * The solar abundances table (from the cooling struct) is arranged as
 * [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe].
 *
 * @param p Pointer to #part struct.
 * @param cooling #cooling_function_data struct.
 * @param ratio_solar (return) Array of ratios to solar abundances.
 */
__attribute__((always_inline)) INLINE float abundance_ratio_to_solar(
    const struct part *p, const struct cooling_function_data *cooling,
    float ratio_solar[colibre_cooling_N_elementtypes]) {

  float totmass = 0., metalmass = 0.;
  int met_index;
  float d_met, logZZsol, ZZsol, Mfrac;
  float log_nx_nH_sol, log_nx_nH_min, log_nx_nH_max, log_nx_nH;

  /* from mass fractions to abundances (nx/nH) */
  for (int i = 0; i < colibre_cooling_N_elementtypes; i++) {

    if ((i != element_S) && (i != element_Ca) && (i != element_OA)) {
      int indx1d = row_major_index_2d(cooling->indxZsol, i,
                                      colibre_cooling_N_metallicity,
                                      colibre_cooling_N_elementtypes);
      Mfrac = p->chemistry_data
                  .smoothed_metal_mass_fraction[element_from_table_to_code(i)];
      ratio_solar[i] =
          Mfrac /
          p->chemistry_data
              .smoothed_metal_mass_fraction[element_from_table_to_code(
                  element_H)] *
          cooling->atomicmass[element_H] * cooling->atomicmass_inv[i] *
          cooling->Abundances_inv[indx1d];
      totmass += Mfrac;
      if (i > element_He) metalmass += Mfrac;
    } else if (i == element_S) {
      /* S scales with Si */
      ratio_solar[element_S] =
          ratio_solar[element_Si] * cooling->S_over_Si_ratio_in_solar;
      int indx1d = row_major_index_2d(cooling->indxZsol, element_Si,
                                      colibre_cooling_N_metallicity,
                                      colibre_cooling_N_elementtypes);
      /* mass fraction S */
      int indxS = row_major_index_2d(cooling->indxZsol, element_S,
                                     colibre_cooling_N_metallicity,
                                     colibre_cooling_N_elementtypes);
      Mfrac = cooling->S_over_Si_ratio_in_solar *
              cooling->atomicmass[element_S] *
              cooling->atomicmass_inv[element_Si] * cooling->Abundances[indxS] *
              cooling->Abundances_inv[indx1d] *
              p->chemistry_data
                  .smoothed_metal_mass_fraction[element_from_table_to_code(
                      element_Si)];
      totmass += Mfrac;
      metalmass += Mfrac;
    } else if (i == element_Ca) {
      /* Ca scales with Si */
      ratio_solar[element_Ca] =
          ratio_solar[element_Si] * cooling->Ca_over_Si_ratio_in_solar;
      int indx1d = row_major_index_2d(cooling->indxZsol, element_Si,
                                      colibre_cooling_N_metallicity,
                                      colibre_cooling_N_elementtypes);
      /* mass fraction Ca*/
      int indxCa = row_major_index_2d(cooling->indxZsol, element_Ca,
                                      colibre_cooling_N_metallicity,
                                      colibre_cooling_N_elementtypes);
      Mfrac = cooling->Ca_over_Si_ratio_in_solar *
              cooling->atomicmass[element_Ca] *
              cooling->atomicmass_inv[element_Si] *
              cooling->Abundances[indxCa] * cooling->Abundances_inv[indx1d] *
              p->chemistry_data
                  .smoothed_metal_mass_fraction[element_from_table_to_code(
                      element_Si)];
      totmass += Mfrac;
      metalmass += Mfrac;
    }
  }

  ZZsol = metalmass / totmass * cooling->Zsol_inv[0];
  logZZsol = log10f(ZZsol);
  /* All other elements (element_OA): scale with metallicity */
  ratio_solar[element_OA] = ZZsol;

  /* at this point ratio_solar is (nx/nH) / (nx/nH)_sol */
  /* to multiply with the tables, we want the individual abundance ratio
   * relative */
  /* to what is used in the tables for each metallicity */

  /* for example: for a metallicity of 1 per cent solar, the metallicity bin */
  /* for logZZsol = -2 has already the reduced cooling rates for each element;
   */
  /* it should therefore NOT be multiplied by 0.01 again */
  /* BUT: if e.g. Carbon is twice as abundant as the solar abundance ratio, */
  /* i.e. nC / nH = 0.02 * (nC/nH)_sol for the overall metallicity of 0.01, */
  /* the Carbon cooling rate is multiplied by 2 */

  for (int i = 0; i < colibre_cooling_N_elementtypes; i++) {
    get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
                 &met_index, &d_met);
    log_nx_nH_min = cooling->LogAbundances[row_major_index_2d(
        met_index, i, colibre_cooling_N_metallicity,
        colibre_cooling_N_elementtypes)];
    log_nx_nH_max = cooling->LogAbundances[row_major_index_2d(
        met_index + 1, i, colibre_cooling_N_metallicity,
        colibre_cooling_N_elementtypes)];
    log_nx_nH_sol = cooling->LogAbundances[row_major_index_2d(
        cooling->indxZsol, i, colibre_cooling_N_metallicity,
        colibre_cooling_N_elementtypes)];
    log_nx_nH =
        (log_nx_nH_min * (1.f - d_met) + log_nx_nH_max * d_met) - log_nx_nH_sol;
    ratio_solar[i] *= exp10f(-log_nx_nH);
  }

  /* at this point ratio_solar is (nx/nH) / (nx/nH)_table [Z], */
  /* the metallicity dependent abundance ratio for solar abundances */

  return logZZsol;
}

/**
 * @brief Computes the extra heat from Helium reionisation at a given redshift.
 *
 * We follow the implementation of Wiersma et al. 2009, MNRAS, 399, 574-600,
 * section. 2. The calculation returns energy in CGS.
 *
 * Note that delta_z is negative.
 *
 * @param z The current redshift.
 * @param delta_z The change in redhsift over the course of this time-step.
 * @param cooling The #cooling_function_data used in the run.
 * @return Helium reionization energy in CGS units.
 */
__attribute__((always_inline)) INLINE double
eagle_helium_reionization_extraheat(
    double z, double delta_z, const struct cooling_function_data *cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Recover the values we need */
  const double z_centre = cooling->He_reion_z_centre;
  const double z_sigma = cooling->He_reion_z_sigma;
  const double heat_cgs = cooling->He_reion_heat_cgs;

  double extra_heat = 0.;

  /* Integral of the Gaussian between z and z - delta_z */
  extra_heat += erf((z - delta_z - z_centre) / (M_SQRT2 * z_sigma));
  extra_heat -= erf((z - z_centre) / (M_SQRT2 * z_sigma));

  /* Multiply by the normalisation factor */
  extra_heat *= heat_cgs * 0.5;

  return extra_heat;
}

/**
 * @brief Computes the log_10 of the temperature corresponding to a given
 * internal energy, hydrogen number density, metallicity and redshift
 *
 * @param log_10_u_cgs Log base 10 of internal energy in cgs.
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param met_index Index along the metallicity dimension.
 * @param d_met Offset between metallicity and table[met_index].
 * @param red_index Index along the redshift dimension.
 * @param d_red Offset between redshift and table[red_index].
 * @param cooling #cooling_function_data structure.
 *
 * @return log_10 of the temperature.
 *
 * TO DO: outside table ranges, it uses at the moment the minimum, maximu value
 */
__attribute__((always_inline)) INLINE double colibre_convert_u_to_temp(
    const double log_10_u_cgs, const float redshift, int n_H_index, float d_n_H,
    int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the internal energy axis */
  int u_index;
  float d_u;

  get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_10_u_cgs,
               &u_index, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  float log_10_T;

  /* Temperature from internal energy */
  log_10_T = interpolation_4d(
      cooling->table.T_from_U, red_index, u_index, met_index, n_H_index, d_red,
      d_u, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density);

  /* Special case for temperatures below the start of the table */
  if (u_index == 0 && d_u == 0.f) {

    /* The temperature is multiplied by u / 10^T[0]
     * where T[0] is the first entry in the table */
    log_10_T += log_10_u_cgs - cooling->Temp[0];
  }

  return log_10_T;
}

/**
 * @brief Computes the log_10 of the internal energy corresponding to a given
 * temperature, hydrogen number density, metallicity and redshift
 *
 * @param log_10_T Log base 10 of temperature in K
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param met_index Index along the metallicity dimension.
 * @param d_met Offset between metallicity and table[met_index].
 * @param red_index Index along the redshift dimension.
 * @param d_red Offset between redshift and table[red_index].
 * @param cooling #cooling_function_data structure.
 *
 * @return log_10 of the internal energy in cgs
 *
 * TO DO: outside table ranges, it uses at the moment the minimum, maximu value
 */

__attribute__((always_inline)) INLINE double colibre_convert_temp_to_u(
    const double log_10_T, const float redshift, int n_H_index, float d_n_H,
    int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the internal energy axis */
  int T_index;
  float d_T;

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_10_T, &T_index,
               &d_T);

  /* Interpolate internal energy table to return internal energy for current
   * temperature (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  float log_10_U;

  /* Internal energy from temperature*/
  log_10_U = interpolation_4d(
      cooling->table.U_from_T, red_index, T_index, met_index, n_H_index, d_red,
      d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_temperature, colibre_cooling_N_metallicity,
      colibre_cooling_N_density);

  return log_10_U;
}

/**
 * @brief Computes the mean particle mass for a given
 * metallicity, temperature, redshift, and density.
 *
 * @param log_T_cgs Log base 10 of temperature in K
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param ZZsol Metallicity relative to the solar value from the tables
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @retura linear electron density in cm-3 (NOT the electron fraction)
 */

INLINE double colibre_meanparticlemass_temperature(
    double log_T_cgs, double redshift, double n_H_cgs, float ZZsol,
    int n_H_index, float d_n_H, int met_index, float d_met, int red_index,
    float d_red, const struct cooling_function_data *restrict cooling) {

  /* Get index of T along the temperature axis */
  int T_index;
  float d_T;

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);

  double mu = interpolation_4d(
      cooling->table.Tmu, red_index, T_index, met_index, n_H_index, d_red, d_T,
      d_met, d_n_H, colibre_cooling_N_redshifts, colibre_cooling_N_temperature,
      colibre_cooling_N_metallicity, colibre_cooling_N_density);

  return mu;
}

/**
 * @brief Computes the electron density for a given element
 * abundance ratio, internal energy, redshift, and density.
 *
 * @param log_u_cgs Log base 10 of internal energy in cgs [erg g-1]
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param ZZsol Metallicity relative to the solar value from the tables
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @retura linear electron density in cm-3 (NOT the electron fraction)
 */

INLINE double colibre_electron_density(
    double log_u_cgs, double redshift, double n_H_cgs, float ZZsol,
    const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
    float d_n_H, int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the internal energy axis */
  int U_index;
  float d_U;

  get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_u_cgs,
               &U_index, &d_U);

  /* n_e / n_H */
  double electron_fraction = interpolation4d_plus_summation(
      cooling->table.Uelectron_fraction, abundance_ratio, element_H,
      colibre_cooling_N_electrontypes - 4, red_index, U_index, met_index,
      n_H_index, d_red, d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  double ne = electron_fraction * n_H_cgs;

  return ne;
}

/**
 * @brief Computes the electron density for a given element
 * abundance ratio, temperature, redshift, and density.
 *
 * @param log_T_cgs Log base 10 of temperature
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param ZZsol Metallicity relative to the solar value from the tables
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @retura linear electron density in cm-3 (NOT the electron fraction)
 */

INLINE double colibre_electron_density_temperature(
    double log_T_cgs, double redshift, double n_H_cgs, float ZZsol,
    const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
    float d_n_H, int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the internal energy axis */
  int T_index;
  float d_T;

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);

  /* n_e / n_H */
  double electron_fraction = interpolation4d_plus_summation(
      cooling->table.Telectron_fraction, abundance_ratio, element_H,
      colibre_cooling_N_electrontypes - 4, red_index, T_index, met_index,
      n_H_index, d_red, d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  double ne = electron_fraction * n_H_cgs;

  return ne;
}

/**
 * @brief Computes the net cooling rate (cooling - heating) for a given element
 * abundance ratio, internal energy, redshift, and density. The unit of the net
 * cooling rate is Lambda / nH**2 [erg cm^3 s-1] and all input values are in
 * cgs. The Compton cooling is not taken from the tables but calculated
 * analytically and added separately
 *
 * @param log_u_cgs Log base 10 of internal energy in cgs [erg g-1]
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @param onlyicool if true / 1 only plot cooling channel icool
 * @param onlyiheat if true / 1 only plot cooling channel iheat
 * @param icool cooling channel to be used
 * @param iheat heating channel to be used
 *
 * Throughout the code: onlyicool = onlyiheat = icool = iheat = 0
 * These are only used for testing: examples/CoolingRates/CoolingRatesCOLIBRE
 */

INLINE double colibre_cooling_rate(
    double log_u_cgs, double redshift, double n_H_cgs,
    const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
    float d_n_H, int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling, int onlyicool,
    int onlyiheat, int icool, int iheat) {

  double cooling_rate, heating_rate, Compton_cooling_rate, temp, logtemp;
  double net_cooling_rate, electron_fraction;

  float weights_cooling[colibre_cooling_N_cooltypes - 2];
  float weights_heating[colibre_cooling_N_heattypes - 2];

  const double zp1 = 1. + redshift;
  const double zp1p2 = zp1 * zp1;
  const double zp1p4 = zp1p2 * zp1p2;

  /* CMB temperature at this redshift */
  const double T_CMB = cooling->T_CMB_0 * zp1;

  /* set weights for cooling rates */
  if (onlyicool == 0) {
    /* standard case */
    for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      if (i < colibre_cooling_N_elementtypes) {
        weights_cooling[i] = abundance_ratio[i];
      } else if (i == cooltype_Compton) {
        weights_cooling[i] =
            0.f; /* added analytically later, do not use value from table*/
      } else {
        weights_cooling[i] = 1.f; /* use same abundances as in the tables */
      }
    }
  } else {
    /* only return cooling rate from one cooling channel: icool */
    for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      weights_cooling[i] = 0.f;
    }
    if (icool >= 0) {
      if (icool < colibre_cooling_N_elementtypes) {
        weights_cooling[icool] = abundance_ratio[icool];
      } else if (icool != cooltype_Compton) {
        weights_cooling[icool] = 1.f;
      }
    }
  }

  /* set weights for heating rates */
  if (onlyiheat == 0) {
    /* standard case */
    for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
      if (i < colibre_cooling_N_elementtypes) {
        weights_heating[i] = abundance_ratio[i];
      } else {
        weights_heating[i] = 1.f; /* use same abundances as in the tables */
      }
    }
  } else {
    /* only return heating rate from one heating channel: iheat */
    for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
      weights_heating[i] = 0.f;
    }
    if (iheat >= 0) {
      if (iheat < colibre_cooling_N_elementtypes) {
        weights_heating[iheat] = abundance_ratio[iheat];
      } else {
        weights_heating[iheat] = 1.f;
      }
    }
  }

  /* Get index of u along the internal energy axis */
  int U_index;
  float d_U;
  get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_u_cgs,
               &U_index, &d_U);

  /* n_e / n_H */
  electron_fraction = interpolation4d_plus_summation(
      cooling->table.Uelectron_fraction, abundance_ratio, element_H,
      colibre_cooling_N_electrontypes - 4, red_index, U_index, met_index,
      n_H_index, d_red, d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  /* Lambda / n_H**2 */
  cooling_rate = interpolation4d_plus_summation(
      cooling->table.Ucooling, weights_cooling, element_H,
      colibre_cooling_N_cooltypes - 3, red_index, U_index, met_index, n_H_index,
      d_red, d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_cooltypes);

  /* Gamma / n_H**2 */
  heating_rate = interpolation4d_plus_summation(
      cooling->table.Uheating, weights_heating, element_H,
      colibre_cooling_N_heattypes - 3, red_index, U_index, met_index, n_H_index,
      d_red, d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_heattypes);

  /* Temperature from internal energy */
  logtemp = interpolation_4d(
      cooling->table.T_from_U, red_index, U_index, met_index, n_H_index, d_red,
      d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density);

  temp = exp10(logtemp);

  Compton_cooling_rate = 0.;
  if (onlyicool == 0 || (onlyicool == 1 && icool == cooltype_Compton)) {
    /* Analytic Compton cooling rate: Lambda_Compton / n_H**2 */
    Compton_cooling_rate = cooling->compton_rate_cgs * (temp - T_CMB) * zp1p4 *
                           electron_fraction / n_H_cgs;
  }

  net_cooling_rate = heating_rate - cooling_rate - Compton_cooling_rate;

  return net_cooling_rate;
}

/**
 * @brief Computes the net cooling rate (cooling - heating) for a given element
 * abundance ratio, temperature, redshift, and density. The unit of the net
 * cooling rate is Lambda / nH**2 [erg cm^3 s-1] and all input values are in
 * cgs. The Compton cooling is not taken from the tables but calculated
 * analytically and added separately
 *
 * @param log_T_cgs Log base 10 of temperature in K
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @param onlyicool if true / 1 only plot cooling channel icool
 * @param onlyiheat if true / 1 only plot cooling channel iheat
 * @param icool cooling channel to be used
 * @param iheat heating channel to be used
 *
 * Throughout the code: onlyicool = onlyiheat = icool = iheat = 0
 * These are only used for testing: examples/CoolingRates/CoolingRatesCOLIBRE
 */

INLINE double colibre_cooling_rate_temperature(
    double log_T_cgs, double redshift, double n_H_cgs,
    const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
    float d_n_H, int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling, int onlyicool,
    int onlyiheat, int icool, int iheat) {

  /* Get index of T along the internal energy axis */
  int T_index;
  float d_T;

  double cooling_rate, heating_rate, Compton_cooling_rate, temp;
  double net_cooling_rate, electron_fraction;

  float weights_cooling[colibre_cooling_N_cooltypes - 2];
  float weights_heating[colibre_cooling_N_heattypes - 2];

  const double zp1 = 1. + redshift;
  const double zp1p2 = zp1 * zp1;
  const double zp1p4 = zp1p2 * zp1p2;

  /* CMB temperature at this redshift */
  const double T_CMB = cooling->T_CMB_0 * zp1;

  int i;

  if (onlyicool == 0) {

    /* set weights for cooling rates */
    for (i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      if (i < colibre_cooling_N_elementtypes)
        weights_cooling[i] = abundance_ratio[i];
      if (i == cooltype_H2)
        weights_cooling[i] = 1.; /* use same H2 abundance as in tables */
      if (i == cooltype_molecules) weights_cooling[i] = 1.;
      if (i == cooltype_HD)
        weights_cooling[i] = 1.; /* use same HD abundance as in tables */
      if (i == cooltype_NetFFH) weights_cooling[i] = 1.;
      if (i == cooltype_NetFFM) weights_cooling[i] = 1.;
      if (i == cooltype_eeBrems)
        weights_cooling[i] = 1.; /* use same electron abundance */
      if (i == cooltype_Compton)
        weights_cooling[i] = 0.; /* added analytically */
      if (i == cooltype_Dust) weights_cooling[i] = 1.;
    }
  } else {
    for (i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      weights_cooling[i] = 0.;
    }
    if (icool >= 0) {
      if (icool < colibre_cooling_N_elementtypes) {
        weights_cooling[icool] = abundance_ratio[icool];
      } else if (icool != cooltype_Compton) {
        weights_cooling[icool] = 1.;
      }
    }
  }

  if (onlyiheat == 0) {
    /* set weights for heating rates */
    for (i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
      if (i < colibre_cooling_N_elementtypes)
        weights_heating[i] = abundance_ratio[i];
      if (i == heattype_H2) weights_heating[i] = 1.;
      if (i == heattype_COdiss) weights_heating[i] = 1.;
      if (i == heattype_CosmicRay) weights_heating[i] = 1.;
      if (i == heattype_UTA) weights_heating[i] = 1.;
      if (i == heattype_line) weights_heating[i] = 1.;
      if (i == heattype_Hlin) weights_heating[i] = 1.;
      if (i == heattype_ChaT) weights_heating[i] = 1.;
      if (i == heattype_HFF) weights_heating[i] = 1.;
      if (i == heattype_Compton) weights_heating[i] = 1.;
      if (i == heattype_Dust) weights_heating[i] = 1.;
    }
  } else {
    for (i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      weights_heating[i] = 0.;
    }
    if (iheat >= 0) {
      if (iheat < colibre_cooling_N_elementtypes) {
        weights_heating[iheat] = abundance_ratio[iheat];
      } else {
        weights_heating[iheat] = 1.;
      }
    }
  }

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);

  /* n_e / n_H */
  electron_fraction = interpolation4d_plus_summation(
      cooling->table.Telectron_fraction, abundance_ratio, element_H,
      colibre_cooling_N_electrontypes - 4, red_index, T_index, met_index,
      n_H_index, d_red, d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  /* Lambda / n_H**2 */
  cooling_rate = interpolation4d_plus_summation(
      cooling->table.Tcooling, weights_cooling, element_H,
      colibre_cooling_N_cooltypes - 3, red_index, T_index, met_index, n_H_index,
      d_red, d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_cooltypes);

  /* Gamma / n_H**2 */
  heating_rate = interpolation4d_plus_summation(
      cooling->table.Theating, weights_heating, element_H,
      colibre_cooling_N_heattypes - 3, red_index, T_index, met_index, n_H_index,
      d_red, d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_heattypes);

  temp = pow(10., log_T_cgs);

  Compton_cooling_rate = 0.;
  if (onlyicool == 0 || (onlyicool == 1 && icool == cooltype_Compton)) {
    /* Analytic Compton cooling rate: Lambda_Compton / n_H**2 */
    Compton_cooling_rate = cooling->compton_rate_cgs * (temp - T_CMB) * zp1p4 *
                           electron_fraction / n_H_cgs;
  }

  net_cooling_rate = heating_rate - cooling_rate - Compton_cooling_rate;

  return net_cooling_rate;
}

#endif /* SWIFT_COLIBRE_COOLING_RATES_H */
