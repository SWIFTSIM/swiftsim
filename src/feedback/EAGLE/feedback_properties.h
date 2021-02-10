/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_FEEDBACK_PROPERTIES_H
#define SWIFT_EAGLE_FEEDBACK_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "feedback.h"
#include "chemistry.h"
#include "hydro_properties.h"

/**
 * @brief Modes of energy injection for SNII feedback
 */
enum SNII_feedback_models {
  SNII_random_ngb_model,       /*< Random neighbour model for SNII feedback */
  SNII_isotropic_model,        /*< Isotropic model of SNII feedback */
  SNII_minimum_distance_model, /*< Minimum-distance model of SNII feedback */
  SNII_minimum_density_model   /*< Minimum-density model of SNII feedback */
};

/**
 * @brief Type of high-density oversampling criterion
 */
enum SNII_oversampling_criterion {
  eagle_SNII_timescale_none,        /*<! No additional criterion */
  eagle_SNII_timescale_gasconsum,   /*<! Gas consumption timescale */
  eagle_SNII_timescale_freefall     /*<! Free-fall timescale */
};

/**
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {

  /*! Yield table mass bins */
  double *mass;

  /*! Yield table metallicity (metal mass fractions) bins */
  double *metallicity;

  /*! Array to store yield table (individual metals produced by the star)
     resampled by IMF mass bins */
  double *yield_IMF_resampled;

  /*! Array to store yield table being read in */
  double *yield;

  /*! Array to store table of ejecta (metals alredy in the stars that are
    ejected) resampled by IMF mass bins */
  double *ejecta_IMF_resampled;

  /*! Array to store table of ejecta being read in */
  double *ejecta;

  /*! Array to store table of total mass released ( metals produced by the star)
    resampled by IMF mass bins */
  double *total_metals_IMF_resampled;

  /* Array to store table of total mass released being read in */
  double *total_metals;
};

/**
 * @brief Stores tables to determine stellar lifetimes. Used for calculation of
 * IMF
 */
struct lifetime_table {

  /* table of masses */
  double *mass;

  /* table of metallicities */
  double *metallicity;

  /* table of lifetimes depending on mass an metallicity */
  double **dyingtime;
};

/**
 * @brief Functional form of the SNIa delay time distribution.
 */
enum eagle_feedback_SNIa_DTD {

  /*! Power-law with slope -1 */
  eagle_feedback_SNIa_DTD_power_law = 1,

  /*! Exponential model (EAGLE default) */
  eagle_feedback_SNIa_DTD_exponential = 2
};

/**
 * @brief Properties of the EAGLE feedback model.
 */
struct feedback_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we doing AGB enrichment? */
  int with_AGB_enrichment;

  /*! Are we doing SNII enrichment? */
  int with_SNII_enrichment;

  /*! Are we doing SNIa enrichment? */
  int with_SNIa_enrichment;

  /*! Are we doing SNII feedback? */
  int with_SNII_feedback;

  /*! Are we doing SNIa feedback? */
  int with_SNIa_feedback;

  /* ------------ Yield tables    ----------------- */

  /* Yield tables for AGB and SNII  */
  struct yield_table yield_AGB;
  struct yield_table yield_SNII;

  /* Arrays of yield tables for SNIa */
  double *yield_SNIa_IMF_resampled;
  double yield_SNIa_total_metals_IMF_resampled;
  double *yields_SNIa;

  /* Arrays for names of elements being tracked for each enrichment channel */
  char **SNIa_element_names;
  char **SNII_element_names;
  char **AGB_element_names;

  /* Array of mass bins for yield calculations */
  double *yield_mass_bins;

  /* Location of yield tables */
  char yield_table_path[200];

  /* ------------- Lifetime tracks   --------------- */

  /* Table of lifetime values */
  struct lifetime_table lifetimes;

  /* ------------- SNII parameters    --------------- */

  /* Array of adjustment factors for SNII  */
  float SNII_yield_factor[chemistry_element_count];

  /* ------------- SNIa parameters    --------------- */

  /* What delay time distribution are we using? */
  enum eagle_feedback_SNIa_DTD SNIa_DTD;

  /*! Normalisation of the SNIa DTD in the exponential model */
  float SNIa_DTD_exp_norm;

  /*! Time-scale of the SNIa decay function in the exponential model in
   * Giga-years */
  float SNIa_DTD_exp_timescale_Gyr;

  /*! Inverse of time-scale of the SNIa decay function in the exponential model
   * in Giga-years */
  float SNIa_DTD_exp_timescale_Gyr_inv;

  /*! Normalisation of the SNIa DTD in the power-law model */
  float SNIa_DTD_power_law_norm;

  /*! Stellar age below which no SNIa explode in Giga-years */
  float SNIa_DTD_delay_Gyr;

  /*! Energy released by one supernova type II in cgs units */
  double E_SNIa_cgs;

  /*! Energy released by one supernova type II in internal units */
  float E_SNIa;

  /* ------------- AGB parameters    ---------------- */

  /*! Specific kinetic energy injected from AGB ejectas (in internal units). */
  float AGB_ejecta_specific_kinetic_energy;

  /* ------------- Conversion factors --------------- */

  /*! Conversion factor from internal mass unit to solar mass */
  double mass_to_solar_mass;

  /*! Conversion factor from internal mass unit to solar mass */
  double solar_mass_to_mass;

  /*! Conversion factor from density in internal units to Hydrogen number
   * density in cgs */
  double rho_to_n_cgs;

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /* ------------- Parameters for IMF --------------- */

  /*! Array to store calculated IMF */
  double *imf;

  /*! Arrays to store IMF mass bins */
  double *imf_mass_bin;

  /*! Arrays to store IMF mass bins (log10)*/
  double *imf_mass_bin_log10;

  /*! Minimal stellar mass considered by the IMF (in solar masses) */
  double imf_min_mass_msun;

  /*! Maximal stellar mass considered by the IMF (in solar masses) */
  double imf_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered by the IMF (in solar masses)
   */
  double log10_imf_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered by the IMF (in solar masses)
   */
  double log10_imf_max_mass_msun;

  /* ------------ SNe feedback properties ------------ */

  /*! SNII feedback model: random, isotropic or minimum distance */
  enum SNII_feedback_models feedback_model;

  /*! Minimal stellar mass considered for SNII feedback (in solar masses) */
  double SNII_min_mass_msun;

  /*! Maximal stellar mass considered for SNII feedback (in solar masses) */
  double SNII_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered for SNII feedback (in solar
   * masses) */
  double log10_SNII_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered for SNII feedback (in solar
   * masses) */
  double log10_SNII_max_mass_msun;

  /*! Number of type II supernovae per solar mass */
  float num_SNII_per_msun;

  /*! Are we sampling the SNII life-times or using a fixed delay? */
  int SNII_sampled_delay;

  /*! Wind delay time for SNII when using a fixed delay */
  double SNII_wind_delay;

  /*! Use variable temperature increase? */
  int SNII_use_variable_delta_T;

  /*! Buffer factor for numerical efficiency temperature */
  double SNII_T_crit_factor;

  /*! Should we use the instantaneous or birth density for determining dT? */
  int SNII_use_instantaneous_density_for_dT;

  /*! Should we use the instantaneous or birth metallicity for finding dT? */
  int SNII_use_instantaneous_Z_for_dT;

  /*! Number of neighbours that should be heatable by SNII */
  double SNII_delta_T_num_ngb_to_heat;

  /*! Maximum temperature increase induced by SNII feedback [Kelvin] */
  double SNII_delta_T_max;

  /*! Switch to enable sub-critical SNII efficiency compensation */
  int SNII_with_energy_compensation;

  /*! Assumed power-law for sub-critical SNII efficiency scaling */
  double SNII_efficiency_zeta;

  /*! Maximum allowed increase in SNII energy to compensate numerical losses */
  double SNII_delta_T_omega_max;

  /*! Maximum allowed sampling suppression factor */
  float SNII_maximum_nu;

  /*! Normalisation SFR density for sampling reduction at high density */
  double SNII_gamma_star;

  /*! Switch to use adaptive sampling reduction based on smoothing length */
  int SNII_sampling_reduction_within_smoothing_length;

  /*! Switch to allow nu < 0, i.e. higher sampling at low density */
  int SNII_with_nu_below_one;

  /*! Switch to increase the sampling criterion according to a time scale */
  enum SNII_oversampling_criterion SNII_with_oversampling_timescale;

  /*! Set omega to enforce sampling, rather than energy conservation? */
  int SNII_omega_by_sampling;

  /* --- Parameters describing the dT_crit and dT_limit planes --- */

  /*! Solar metallicity value */
  float Z_Sol;

  /*! Minimum allowed dT_crit */
  float SNII_dTcrit_floor;

  /*! Density exponent for dT_crit plane */
  float SNII_dTcrit_exp_nH;

  /*! Metallicity exponent for dT_crit plane */
  float SNII_dTcrit_exp_Z;

  /*! Normalisation temperature for dT_crit plane */
  float SNII_dTcrit_norm;

  /*! Density exponent for dT_limit plane */
  float SNII_dTlimit_exp_nH;

  /*! Metallicity exponent for dT_limit plane */
  float SNII_dTlimit_exp_Z;

  /*! Minimum metallicity for dT_limit plane */
  float SNII_dT_Zmin;

  /*! Normalisation temperature for dT_limit plane */
  float SNII_dTlimit_norm;

  /*! Minimum allowed efficiency for compromise dT finding */
  float SNII_efficiency_eta_min;

  /*! (Constant) temperature increase induced by SNe feedback */
  double SNe_deltaT_desired;

  /*! Energy released by one supernova type II in cgs units */
  double E_SNII_cgs;

  /*! Energy released by one supernova type II in internal units */
  float E_SNII;

  /*! Minimal energy fraction for supernova type II feedback */
  double f_E_min;

  /*! Maximal energy fraction for supernova type II feedback */
  double f_E_max;

  /*! Pivot point for the metallicity dependance of the feedback energy fraction
   * model */
  double Z_0;

  /*! Pivot point for the density dependance of the feedback energy fraction
   * model */
  double n_0_cgs;

  /*! Slope of the density dependance of the feedback energy fraction model */
  double n_n;

  /*! Slope of the metallicity dependance of the feedback energy fraction model
   */
  double n_Z;

  /*! Are we using the birth density to compute f_th or the properties at
   * feedback time? */
  int use_birth_density_for_f_th;

  /*! Are we using the birth metallicity to compute f_th or the properties at
   * feedback time? */
  int use_birth_Z_for_f_th;

  /* ------------ Enrichment sampling properties ------------ */

  /*! Star age above which the enrichment will be downsampled (in internal
   * units) */
  double stellar_evolution_age_cut;

  /*! Number of time-steps in-between two enrichment events */
  int stellar_evolution_sampling_rate;
};

void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo);

#endif /* SWIFT_EAGLE_FEEDBACK_PROPERTIES_H */
