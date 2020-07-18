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
#ifndef SWIFT_COLIBRE_FEEDBACK_PROPERTIES_H
#define SWIFT_COLIBRE_FEEDBACK_PROPERTIES_H

#include "SNIa_DTD_struct.h"
#include "chemistry.h"
#include "hydro_properties.h"
//#include "dust_properties.h" <-- BREAKS CODE, WHY?
//#include "dust.h" */

/*! Number of elements to be read from the yield tables */
#define enrichment_of_N_elements_from_yield_tables 9

/**
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {

  /* Yield table mass bins */
  double *mass;

  /* Yield table metallicity bins */
  double *metallicity;

  /* Array to store yield table resampled by IMF mass bins */
  double *yield_IMF_resampled;

  /* Array to store yield table being read in */
  double *yield;

  /* Array to store table of ejecta resampled by IMF mass bins */
  double *ejecta_IMF_resampled;

  /* Array to store table of ejecta being read in */
  double *ejecta;

  /* Array to store table of total mass released resampled by IMF mass bins */
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

  /* table of lifetimes depending on mass and metallicity */
  double **dyingtime;
};

/**
 * @brief Properties of the COLIBRE feedback model.
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

  /*! Are we doing r-process enrichment? */
  int with_r_process_enrichment;

  /*! Are we doing HII regions? */
  int with_HIIRegions;

  /*! Are we doing Stellar winds? */
  int with_StellarWinds;

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
  float SNII_yield_factor[enrichment_of_N_elements_from_yield_tables];

  /* ------------- SNIa parameters    --------------- */

  /*! Energy released by one supernova type Ia in cgs units */
  double E_SNIa_cgs;

  /*! Energy released by one supernova type Ia in internal units */
  float E_SNIa;

  /* SNIa DTD struct with information about the DTD */
  struct SNIa_delay_time_distribution dtd_data;

  /* ----------- SNeIa feedback properties -------------- */

  /* Temperature increase induced by SNIe feedback */
  float SNIa_deltaT_desired;

  /* Energy fraction for supernova type Ia feedback */
  float SNIa_f_E;

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

  /*! Conversion factor for momentum to cgs */
  double Momentum_to_cgs;

  /*! Conversion factor from Myr to sec */
  double Myr_to_sec;

  /*! Conversion factor from sec to Myr */
  double sec_to_Myr;

  /* ------------- Parameters for IMF --------------- */

  /*! Array to store calculated IMF */
  float *imf;

  /*! Arrays to store IMF mass bins */
  float *imf_mass_bin;

  /*! Arrays to store IMF mass bins (log10)*/
  float *imf_mass_bin_log10;

  /*! Minimal stellar mass considered by the IMF (in solar masses) */
  float imf_min_mass_msun;

  /*! Maximal stellar mass considered by the IMF (in solar masses) */
  float imf_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered by the IMF (in solar masses)
   */
  float log10_imf_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered by the IMF (in solar masses)
   */
  float log10_imf_max_mass_msun;

  /* ------------ SNeII feedback properties ------------ */

  /*! Log 10 of the minimal stellar mass considered for SNII feedback (in solar
   * masses) */
  float log10_SNII_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered for SNII feedback (in solar
   * masses) */
  float log10_SNII_max_mass_msun;

  /*! Number of type II supernovae per solar mass */
  float num_SNII_per_msun;

  /*! Wind delay time for SNII */
  double SNII_wind_delay;

  /*! Temperature increase induced by SNIIe feedback */
  float SNII_deltaT_desired;

  /*! Fraction of SNII energy injected in kinetic form */
  float SNII_f_kinetic;

  /*! Kick velocity in SNII feedback in internal units */
  float SNII_delta_v;

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

  /* Timescale above which stars no longer inject momentum in Myr */
  double SW_max_age_Myr;

  /* Desired delta_v in km/s of particles suject to the wind. */
  /* higher values makes less likely to kick particles. */
  double delta_v;

  /* ------------ r-process enrichment properties ------------ */
  /* Number of neutron star mergers per unit of Msolar */
  double NSM_per_Msun;

  /* Amount of europium (in units of Msolar) relesed by NSM */
  double yield_Eu_from_NSM;

  /* Number of CEJSN per unit of Msolar */
  double CEJSN_per_Msun;

  /* Amount of europium (in units of Msolar) relesed by CEJSN */
  double yield_Eu_from_CEJSN;

  /* Number of collapsar per unit of Msolar */
  double collapsar_per_Msun;

  /* Amount of europium (in units of Msolar) relesed by collapsar */
  double yield_Eu_from_collapsar;

  /*! Log 10 of the maximal stellar mass that will turn into a collapsar (in
   * solar masses) */
  float log10_collapsar_max_mass_msun;

  /*! Log 10 of the minimal stellar mass that will turn into a collapsar (in
   * solar masses) */
  float log10_collapsar_min_mass_msun;

  /* ------------ Early feedback properties ------------ */

  /* Location of early feedback tables */
  char early_feedback_table_path[200];

  /* Maximum age in Myr of star particle to build HII region */
  float HIIregion_max_age_Myr;

  /* Time between rebuilding the HII region in Myr */
  float HIIregion_dt_Myr;

  /* Recombination coefficient in cgs units [cm3 s-1]*/
  float alpha_caseb_recomb;

  /*! Energy floor of the HII region for injection (internal units) */
  float HII_u;

  /* Number of age bins */
  int HII_nr_agebins;

  /* Number of metallicity bins */
  int HII_nr_metbins;

  /* Metallicity bins (log Z, metal mass fractions) from BPASS */
  float *HII_log10_Zbins;

  /* Age bins (star age in Myr) */
  float *HII_agebins;

  /* Cumululative number of ionizing photons per g stellar mass
   * dimension [HII_nr_metbins, HII_nr_agebins] */
  float *HII_log10_Qcum;

  /* Cumulative momentum input per g stellar mass from stellar winds
   * dimension [HII_nr_metbins, HII_nr_agebins] */
  float *SW_log10_Pcum;

  /* Minimum metallicity in the early feedback tables */
  double Zmin_early_fb;

  /* Maximum metallicity in the early feedback tables */
  double Zmax_early_fb;

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
                         const struct cosmology *cosmo/*, struct dustevo_props* dp*/);

#endif /* SWIFT_COLIBRE_FEEDBACK_PROPERTIES_H */
