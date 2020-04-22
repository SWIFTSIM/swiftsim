#ifndef SWIFT_COLIBRE_TABLES_RESTRICT_H
#define SWIFT_COLIBRE_TABLES_RESTRICT_H

/**
 * @file src/cooling/CHIMES/colibre_tables_restrict.h
 * @brief COLIBRE cooling tables - restricted.
 */

/* Config parameters */
#include "../config.h"

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "exp10.h"

/**
 * @brief Compute ratio of mass fraction to solar mass fraction
 * for each element carried by a given particle.
 *
 * The solar abundances are taken from the tables themselves.
 *
 * The COLIBRE chemistry model does not track S and Ca. We assume
 * that their abundance with respect to solar is the same as
 * the ratio for Si.
 *
 * The other un-tracked elements are scaled with the total metallicity.
 *
 * We optionally apply a correction if the user asked for a different
 * ratio.
 *
 * We also re-order the elements such that they match the order of the
 * tables. This is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, OA].
 *
 * The solar abundances table (from the cooling struct) is arranged as
 * [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe].
 *
 * @param p Pointer to #part struct.
 * @param table Colibre cooling table structure.
 * @param ratio_solar (return) Array of ratios to solar abundances.
 */
__attribute__((always_inline)) INLINE static float abundance_ratio_to_solar(
    const struct part *p, const struct colibre_cooling_tables *table,
    float ratio_solar[colibre_cooling_N_elementtypes]) {

  /* Get the particle's metal mass fractions (M_x / M) */
  // const float *Z_mass_frac =
  // chemistry_get_metal_mass_fraction_for_cooling(p);  static float
  // *Z_mass_frac = chemistry_get_metal_mass_fraction_for_cooling(p);
  float const *Z_mass_frac = chemistry_get_metal_mass_fraction_for_cooling(p);

  /* Convert mass fractions to abundances (nx/nH) and compute metal mass */
  float totmass = 0., metalmass = 0.;
  for (enum colibre_cooling_element elem = element_H; elem < element_OA;
       elem++) {

    /* Normal elements: Get the abundance from the particle carried arrays */
    if ((elem != element_S) && (elem != element_Ca)) {

      const int indx1d = cooling_row_major_index_2d(
          table->indxZsol, elem, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes);

      const float Mfrac = Z_mass_frac[element_from_table_to_code(elem)];

      /* ratio_X = ((M_x/M) / (M_H/M)) * (m_H / m_X) * (1 / Z_sun_X) */
      ratio_solar[elem] =
          (Mfrac / Z_mass_frac[element_from_table_to_code(element_H)]) *
          table->atomicmass[element_H] * table->atomicmass_inv[elem] *
          table->Abundances_inv[indx1d];

      totmass += Mfrac;
      if (elem > element_He) metalmass += Mfrac;

      /* Special case: S scales with Si */
    } else if (elem == element_S) {

      ratio_solar[element_S] =
          ratio_solar[element_Si] * table->S_over_Si_ratio_in_solar;

      const int indx1d = cooling_row_major_index_2d(
          table->indxZsol, element_Si, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes);
      /* mass fraction S */
      const int indxS = cooling_row_major_index_2d(
          table->indxZsol, element_S, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes);

      /* ratio_S = ((M_S/M) / (M_H/M)) * (m_H / m_S) * (1 / Z_sun_S) */
      const float Mfrac =
          table->S_over_Si_ratio_in_solar * table->atomicmass[element_S] *
          table->atomicmass_inv[element_Si] * table->Abundances[indxS] *
          table->Abundances_inv[indx1d] *
          Z_mass_frac[element_from_table_to_code(element_Si)];

      totmass += Mfrac;
      metalmass += Mfrac;

      /* Special case: Ca scales with Si */
    } else if (elem == element_Ca) {

      ratio_solar[element_Ca] =
          ratio_solar[element_Si] * table->Ca_over_Si_ratio_in_solar;

      const int indx1d = cooling_row_major_index_2d(
          table->indxZsol, element_Si, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes);
      /* mass fraction Ca */
      const int indxCa = cooling_row_major_index_2d(
          table->indxZsol, element_Ca, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes);

      /* ratio_Ca = ((M_Ca/M) / (M_H/M)) * (m_H / m_Ca) * (1 / Z_sun_Ca) */
      const float Mfrac =
          table->Ca_over_Si_ratio_in_solar * table->atomicmass[element_Ca] *
          table->atomicmass_inv[element_Si] * table->Abundances[indxCa] *
          table->Abundances_inv[indx1d] *
          Z_mass_frac[element_from_table_to_code(element_Si)];

      totmass += Mfrac;
      metalmass += Mfrac;

    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Invalid element!");
#endif
    }
  }

  /* Get total metallicity [Z/Z_sun] from the particle data */
  const float Z_total =
      (float)chemistry_get_total_metal_mass_fraction_for_cooling(p);
  float ZZsol = Z_total * table->Zsol_inv[0];
  if (ZZsol <= 0.0f) ZZsol = FLT_MIN;
  const float logZZsol = log10f(ZZsol);

  /* All other elements (element_OA): scale with metallicity */
  ratio_solar[element_OA] = ZZsol;

  /* Get index and offset from the metallicity table corresponding to this
   * metallicity */
  int met_index;
  float d_met;

  cooling_get_index_1d(table->Metallicity, colibre_cooling_N_metallicity,
                       logZZsol, &met_index, &d_met);

  /* At this point ratio_solar is (nx/nH) / (nx/nH)_sol.
   * To multiply with the tables, we want the individual abundance ratio
   * relative to what is used in the tables for each metallicity */

  /* For example: for a metallicity of 1 per cent solar, the metallicity bin
   * for logZZsol = -2 has already the reduced cooling rates for each element;
   * it should therefore NOT be multiplied by 0.01 again.
   *
   * BUT: if e.g. Carbon is twice as abundant as the solar abundance ratio,
   * i.e. nC / nH = 0.02 * (nC/nH)_sol for the overall metallicity of 0.01,
   * the Carbon cooling rate is multiplied by 2
   *
   * We only do this if we are not in the primodial metallicity bin in which
   * case the ratio to solar should be 0.
   */

  for (int i = 0; i < colibre_cooling_N_elementtypes; i++) {

    /* Are we considering a metal and are *not* in the primodial metallicity
     * bin? Or are we looking at H or He? */
    if ((met_index > 0) || (i == element_H) || (i == element_He)) {

      /* Get min/max abundances */
      const float log_nx_nH_min =
          table->LogAbundances[cooling_row_major_index_2d(
              met_index, i, colibre_cooling_N_metallicity,
              colibre_cooling_N_elementtypes)];

      const float log_nx_nH_max =
          table->LogAbundances[cooling_row_major_index_2d(
              met_index + 1, i, colibre_cooling_N_metallicity,
              colibre_cooling_N_elementtypes)];

      /* Get solar abundances */
      const float log_nx_nH_sol =
          table->LogAbundances[cooling_row_major_index_2d(
              table->indxZsol, i, colibre_cooling_N_metallicity,
              colibre_cooling_N_elementtypes)];

      /* Interpolate ! (linearly in log-space) */
      const float log_nx_nH =
          (log_nx_nH_min * (1.f - d_met) + log_nx_nH_max * d_met) -
          log_nx_nH_sol;

      ratio_solar[i] *= exp10f(-log_nx_nH);

    } else {

      /* Primordial bin --> Z/Z_sun is 0 for that element */
      ratio_solar[i] = 0.0;
    }
  }

  /* at this point ratio_solar is (nx/nH) / (nx/nH)_table [Z],
   * the metallicity dependent abundance ratio for solar abundances.
   * We also return the total metallicity */

  return logZZsol;
}

#endif  // SWIFT_COLIBRE_TABLES_RESTRICT_H
