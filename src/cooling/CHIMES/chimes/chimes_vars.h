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

#include <stdio.h>
#include <stdlib.h>

/*************************
 ** SWIFT-specific code **
 *************************/ 
#include "cooling/CHIMES/colibre_tables.h" 
/** End SWIFT-specific **/  

/*!< Electron mass. Units: g. */ 
#define ELECTRON_MASS 9.1093829e-28 

/*!< Proton mass. Units: g. */ 
#define PROTON_MASS 1.6726218e-24 

/*!< Pi. */ 
#define PI            3.1415927 

/*!< Speed of light. Units: cm / s. */ 
#define LIGHTSPEED    3.0e10 

/*!< Boltzmann constant. Units: erg / K. */ 
#define BOLTZMANNCGS     1.3806e-16 

/*!< Boltzmann constant. Units: eV / K. */ 
#define BOLTZMANN_EVK   8.6173324e-5 

/*!< For dust processes involving G0, e.g. photoelectric heating, we attenuate G0 by exp(- G0_gamma * Av). */ 
#define G0_GAMMA 2.77 

/*!< The total number of species in the full network. */ 
#define TOTSIZE	  157 

/*!< Total number of molecules in the full network. */ 
#define MOLSIZE   20 

/*!< Effective dust absorption and scattering cross section, to calculate Av from N_Htot. Units: cm^2 / H atom. */ 
#define DUSTEFFSIZE 4.0e-22 

/*!< Dust surface area for H2 formation. Units: cm^2 / H atom. */ 
#define DUST_CROSS_SECTION 1.0e-21 

/*!< All reaction and cooling rates are limited to this maximum temperature. */ 
#define MAX_TEMPERATURE_FOR_RATES  2.0e9 

/*!< Maximum number of UV spectra. */ 
#define MAX_UV_SPECTRA 20 

/*!< Exclude species when their corresponding element abundance is below this threshold. */ 
#define METALS_MINIMUM_THRESHOLD 1.0e-30

/*!< The maximum number of steps in the CVODE solver. */ 
#define MAXSTEPS 1e5 

/*!< Defines whether to use single or double precision throughout the chemistry solver. */ 
typedef double ChimesFloat; 

/** 
 * Structure containing the variables that 
 * are specific to each gas particle/cell. 
 */ 
struct gasVariables
{
  /* NOTE: all abundances are in
   * the form ni / nHtot, i.e. the
   * ratio of number densities of
   * species i to total hydrogen. */
  ChimesFloat element_abundances[10];    /*!< Element abundances, in the order: He, C, N, O, Ne, Mg, Si, S, Ca, Fe. */
  ChimesFloat nH_tot;                    /*!< Hydrogen number density. Units: cm^-3. */
  ChimesFloat temperature;               /*!< Gas temperature. Units: K. */
  ChimesFloat TempFloor;                 /*!< Minimum temperature. Units: K. */ 
  ChimesFloat divVel;                    /*!< Velocity divergence. Units: s^-1 */
  ChimesFloat doppler_broad;             /*!< Doppler broadening from turbulence only (thermal broadening is added separately). Units: km s^-1. */
  ChimesFloat *isotropic_photon_density; /*!< Isotropic photon density (i.e. flux divided by the speed of light) for each spectrum. Units: cm^-3 */
  ChimesFloat *G0_parameter;             /*!< Strength of the 6-13.6 eV FUV band for each spectrum. Units: Habing. */ 
  ChimesFloat *H2_dissocJ;               /*!< n / (isotropic_photon_density * c) for each spectrum, where n is photon number density in the 12.24-13.51 eV band. */
  ChimesFloat cr_rate;                   /*!< HI ionisation rate from cosmic rays. Units: s^-1. */ 
  ChimesFloat metallicity;               /*!< Metallicity relative to solar (where Z_sol = 0.0129). */
  ChimesFloat dust_ratio;                /*!< Dust-to-gas ratio relative to the Milky Way. */ 
  ChimesFloat cell_size;                 /*!< Shielding length. Units: cm. */
  ChimesFloat hydro_timestep;            /*!< Total time over which to integrate the chemistry and cooling. Units: s. */ 
  int ForceEqOn;                         /*!< 0 - Evolve chemistry in non-eq; 1 - set abundances to eqm from tables. */ 
  int ThermEvolOn;                       /*!< 0 - Hold temperature fixed; 1 - evolve temperature. */ 
  int InitIonState;                      /*!< Sets initial ionisation state if using #initialise_gas_abundances(). */ 
  ChimesFloat constant_heating_rate;     /*!< Extra heating term to add to the radiative cooling rates (positive for heating). Units: erg s^-1 cm^-3. */
  ChimesFloat *abundances;               /*!<  Abundance array, defined for species i as n_i / n_H. */ 

/*************************
 ** SWIFT-specific code **
 *************************/ 
  float abundance_ratio[colibre_cooling_N_elementtypes];  /*!< Element abundance ratios for the Colibre cooling tables, used in hybrid cooling. */ 
/** End SWIFT-specific **/ 
};

/**  
 * Structure containing the global variables 
 * that govern the behaviour of the chemistry 
 * solver. 
 */ 
struct globalVariables
{
  char MainDataTablePath[500];                 /*!< Path to the chimes_main_data.hdf5 data file. */ 
  char PhotoIonTablePath[MAX_UV_SPECTRA][500]; /*!< Array of strings containing the paths to the cross sections tables, one for each UV spectrum. */ 
  char EqAbundanceTablePath[500];              /*!< Path to the equilibrium abundance table. */ 
  int cellSelfShieldingOn;                     /*!< 0 - switch off self-shielding; 1 - switch on self-shielding. */ 
  int N_spectra;                               /*!< The number of UV spectra. */ 
  int redshift_dependent_UVB_index;            /*!< Specifies which of the UV spectra corresponds to the redshift-dependent UVB. */ 
  int use_redshift_dependent_eqm_tables;       /*!< Use redshift-dependent eqm abundance tables if this flag is set to 1. */ 
  ChimesFloat redshift;                        /*!< Current redshift. */ 
  ChimesFloat reionisation_redshift;           /*!< Before this redshift, the redshift-dependent UVB is set to zero. */ 
  int StaticMolCooling;                        /*!< For CO and H2O cooling, use: 0 - thermal velocity dispersion; 1 - velocity divergence. */ 
  ChimesFloat T_mol;                           /*!< Disable molecules above this temperature. */ 
  ChimesFloat grain_temperature;               /*!< Grain temperature used in the H2 formation rate on dust. */ 
  ChimesFloat cmb_temperature;                 /*!< CMB temperature used for Compton cooling. */ 
  ChimesFloat relativeTolerance;               /*!< Relative tolerance used in the CVODE integration. */ 
  ChimesFloat absoluteTolerance;               /*!< Absolute tolerance used in the CVODE integration. */ 
  ChimesFloat thermalAbsoluteTolerance;        /*!< Absolute tolerance used for integrating the temperature in CVODE. */ 
  ChimesFloat explicitTolerance;               /*!< Tolerance below which we will use the explicit solution. */ 
  int element_included[9];                     /*!< Array of flags to exclude (0) or include (1) each metal element. */ 
  int speciesIndices[TOTSIZE];                 /*!< Maps the position of each species in the abundance array. */ 
  int totalNumberOfSpecies;                    /*!< Total number of species included in the network. */ 
  int scale_metal_tolerances;                  /*!< Scale the absolute tolerances by the corresponding element abundance. */ 
  int chimes_debug;                            /*!< If set to 1, print gasVariables if CVODE returns an error or warning message. */ 

/*************************
 ** SWIFT-specific code **
 *************************/ 
  const struct colibre_cooling_tables *colibre_table;  /*!< Colibre cooling table, for hybrid cooling. */ 
  int hybrid_cooling_mode;                             /*!< 0 - do not use hybrid cooling; 1 - use hybrid cooling. */ 
/** End SWIFT-specific **/ 
}; 

/** 
 * Structure containing information 
 * about each non-equilibrium species 
 */
struct Species_Structure 
{
  int include_species;           /*!< Flag indicating whether the species is included in the network. */ 
  ChimesFloat element_abundance; /*!< Abundance of the corresponding element. */ 
  ChimesFloat creation_rate;     /*!< Sum of the rates of all reactions that create this species. */ 
  ChimesFloat destruction_rate;  /*!< Sum of the rates of all reactions that destroy this species. */ 
};

/** 
 * Structure containing the table bins 
 * for all of the rate tables.
 */ 
extern struct chimes_table_bins_struct 
{ 
  int N_Temperatures;                        /*!< Number of temperature bins. */              
  ChimesFloat *Temperatures;                 /*!< Temperature array. */ 
  int N_Dust_Temperatures;                   /*!< Number of dust temperature bins. */ 
  ChimesFloat *Dust_Temperatures;            /*!< Dust temperature array. */ 
  int N_Psi;                                 /*!< Number of bins for the Psi parameter. */ 
  ChimesFloat *Psi;                          /*!< Psi parameter array. */ 
  int N_secondary_cosmic_ray_xHII;           /*!< Number of xHII bins for secondary cosmic ray reactions. */ 
  ChimesFloat *secondary_cosmic_ray_xHII;    /*!< xHII array for secondary cosmic ray reactions. */ 
  int N_Column_densities;                    /*!< Number of column density bins. */ 
  ChimesFloat *Column_densities;             /*!< Column density array. */ 
  int N_H2self_column_densities;             /*!< Number of H2 column density bins for H2 self-shielding. */ 
  ChimesFloat *H2self_column_densities;      /*!< H2 column density array for H2 self-shielding. */ 
  int N_b_turbulence;                        /*!< Number of turbulent doppler broadening bins for H2 self-shielding. */ 
  ChimesFloat *b_turbulence;                 /*!< Turbulent doppler broadening array for H2 self-shielding. */ 
  int N_COself_column_densities;             /*!< Number of CO column density bins for CO self-shielding. */  
  ChimesFloat *COself_column_densities;      /*!< CO column density array for CO self-shielding. */ 
  int N_H2CO_column_densities;               /*!< Number of H2 column density bins for CO shielding. */ 
  ChimesFloat *H2CO_column_densities;        /*!< H2 column density array for CO shielding. */ 
  int N_mol_cool_Temperatures;               /*!< Number of temperature bins for molecular cooling. */ 
  ChimesFloat *mol_cool_Temperatures;        /*!< Temperature array for molecular cooling. */ 
  int N_CO_cool_rot_ColumnDensities;         /*!< Number of CO column density bins for CO rotational cooling. */ 
  ChimesFloat *CO_cool_rot_ColumnDensities;  /*!< CO column density array for CO rotational cooling. */ 
  int N_CO_cool_vib_ColumnDensities;         /*!< Number of CO column density bins for CO vibrational cooling. */ 
  ChimesFloat *CO_cool_vib_ColumnDensities;  /*!< CO column density array for CO vibrational cooling. */ 
  int N_H2O_cool_hiT_Temperatures;           /*!< Number of temperature bins for H2O cooling at high temperatures. */ 
  ChimesFloat *H2O_cool_hiT_Temperatures;    /*!< Temperature array for H2O cooling at high temperatures. */ 
  int N_H2O_cool_lowT_Temperatures;          /*!< Number of temperature bins for H2O cooling at low temperatures. */ 
  ChimesFloat *H2O_cool_lowT_Temperatures;   /*!< Temperature array for H2O cooling at low temperatures. */ 
  int N_H2O_cool_rot_ColumnDensities;        /*!< Number of H2O column density bins for H2O rotational cooling. */ 
  ChimesFloat *H2O_cool_rot_ColumnDensities; /*!< H2O column density array for H2O rotational cooling. */ 
  int N_H2O_cool_vib_ColumnDensities;        /*!< Number of H2O column density bins for H2O vibrational cooling. */ 
  ChimesFloat *H2O_cool_vib_ColumnDensities; /*!< H2O column density array for H2O vibrational cooling. */ 
  int N_cool_2d_Temperatures;                /*!< Number of temperature bins for the 2-d cooling channels. */ 
  ChimesFloat *cool_2d_Temperatures;         /*!< Temperature array for the 2-d cooling channels. */ 
  int N_cool_hiT_2d_Temperatures;            /*!< Number of temperature bins for the 2-d cooling channels at high temperature. */ 
  ChimesFloat *cool_hiT_2d_Temperatures;     /*!< Temperature array for the 2-d cooling channels at high temperature. */ 
  int N_cool_2d_ElectronDensities;           /*!< Number of electron density bins for the 2-d cooling channels. */ 
  ChimesFloat *cool_2d_ElectronDensities;    /*!< Electron density array for the 2-d cooling channels. */ 
  int N_cool_4d_Temperatures;                /*!< Number of temperature bins for the 4-d cooling channels. */ 
  ChimesFloat *cool_4d_Temperatures;         /*!< Temperature array for the 4-d cooling channels. */ 
  int N_cool_hiT_4d_Temperatures;            /*!< Number of temperature bins for the 4-d cooling channels at high temperature. */ 
  ChimesFloat *cool_hiT_4d_Temperatures;     /*!< Temperature array for the 4-d cooling channels at high temperature. */ 
  int N_cool_4d_HIDensities;                 /*!< Number of HI density bins for the 4-d cooling channels. */ 
  ChimesFloat *cool_4d_HIDensities;          /*!< HI density array for the 4-d cooling channels. */ 
  int N_cool_4d_ElectronDensities;           /*!< Number of electron density bins for the 4-d cooling channels. */ 
  ChimesFloat *cool_4d_ElectronDensities;    /*!< Electron density array for the 4-d cooling channels. */ 
  int N_cool_4d_HIIDensities;                /*!< Number of HII density bins for the 4-d cooling channels. */ 
  ChimesFloat *cool_4d_HIIDensities;         /*!< HII density array for the 4-d cooling channels. */ 
} chimes_table_bins; 

/** 
 * Structure containing the rates for 
 * the T-dependent reaction group.
 */ 
extern struct chimes_T_dependent_struct 
{ 
  int N_reactions[2];                          /*!< Number of reactions, excluding and including molecules. */ 
  int **reactants;                             /*!< Reactants for each reaction. */ 
  int **products;                              /*!< Products for each reaction. */ 
  int **element_incl;                          /*!< Flags indicating which elements are required for each reaction. */ 
  int *molecular_flag;                         /*!< Flags indicating whether each reaction involves molecules. */ 
  int H2_collis_dissoc_heating_reaction_index; /*!< Index pointing to the H2 collisional dissociation reaction. Needed for the associated heating. */ 
  int H2_form_heating_reaction_index;          /*!< Index pointing to the H2 gas-phase formation reaction. Needed for the associated heating. */ 
  ChimesFloat **rates;                         /*!< Arrays containing the rate coefficients as a function of temperature for each reaction. */ 
} chimes_table_T_dependent; 

/** 
 * Structure containing the rates for 
 * the constant reaction group.
 */ 
extern struct chimes_constant_struct 
{ 
  int N_reactions[2];                 /*!< Number of reactions, excluding and including molecules. */ 
  int **reactants;                    /*!< Reactants for each reaction. */ 
  int **products;                     /*!< Products for each reaction. */ 
  int **element_incl;                 /*!< Flags indicating which elements are required for each reaction. */ 
  int *molecular_flag;                /*!< Flags indicating whether each reaction involves molecules. */ 
  int H2_form_heating_reaction_index; /*!< Index pointing to the H2 gas-phase formation reaction. Needed for the associated heating. */ 
  ChimesFloat *rates;                 /*!< Array containing the constant rate coefficients for each reaction. */ 
} chimes_table_constant; 

/** 
 * Structure containing the rates for 
 * the recombination_AB reaction group.
 */ 
extern struct chimes_recombination_AB_struct 
{ 
  int N_reactions[2];            /*!< Number of reactions, excluding and including molecules. */ 
  int **reactants;               /*!< Reactants for each reaction. */ 
  int *products;                 /*!< Product for each reaction. */ 
  int **element_incl;            /*!< Flags indicating which elements are required for each reaction. */ 
  int *molecular_flag;           /*!< Flags indicating whether each reaction involves molecules. */ 
  ChimesFloat ***rates;          /*!< Arrays containing the rate coefficients as a function of temperature, for case A and case B recombinations separately.*/ 
} chimes_table_recombination_AB; 

/** 
 * Structure containing the rates for 
 * the grain_recombination reaction group.
 */ 
extern struct chimes_grain_recombination_struct 
{ 
  int N_reactions[2];   /*!< Number of reactions, excluding and including molecules. */ 
  int **reactants;      /*!< Reactants for each reaction. */ 
  int *products;        /*!< Products for each reaction. */ 
  int **element_incl;   /*!< Flags indicating which elements are required for each reaction. */ 
  ChimesFloat ***rates; /*!< Arrays containing the rate coefficients as a function of temperature and Psi for each reaction. */ 
} chimes_table_grain_recombination; 

/** 
 * Structure containing the rates for 
 * the cosmic_ray reaction group.
 */ 
extern struct chimes_cosmic_ray_struct 
{ 
  int N_reactions[2];            /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;                /*!< Reactants for each reaction. */ 
  int **products;                /*!< Products for each reaction. */ 
  int **element_incl;            /*!< Flags indicating which elements are required for each reaction. */ 
  int *molecular_flag;           /*!< Flags indicating whether each reaction involves molecules. */ 
  int *secondary_base_reaction;  /*!< Indices of the HI and HeI CR reactions. Needed for secondary ionisations. */ 
  ChimesFloat **secondary_ratio; /*!< Ratio of secondary to primary CR ionisation of HI and HeI, as a function of xHII. */ 
  ChimesFloat *rates;            /*!< Array containing the CR rates relative to HI for each reaction. */ 
} chimes_table_cosmic_ray; 

/** 
 * Structure containing the rates for 
 * the CO_cosmic_ray reaction group.
 */ 
extern struct chimes_CO_cosmic_ray_struct 
{ 
  int N_reactions[2];  /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;      /*!< Reactant for each reaction. */ 
  int **products;      /*!< Products for each reaction. */ 
  int **element_incl;  /*!< Flags indicating which elements are required for each reaction. */ 
  ChimesFloat **rates; /*!< Arrays containing the rate as a function of temperature for each reaction. */ 
} chimes_table_CO_cosmic_ray; 

/** 
 * Structure containing the rates for 
 * the H2_dust_formation (on dust grains) 
 * reaction group.
 */ 
extern struct chimes_H2_dust_formation_struct 
{ 
  int *reactants;      /*!< Reactants. */ 
  int *products;       /*!< Products. */ 
  ChimesFloat **rates; /*!< Array of the rate coefficient as a function of gas and dust temperature. */ 
} chimes_table_H2_dust_formation; 

/** 
 * Structure containing the rates for 
 * the H2_collis_dissoc reaction group.
 */ 
extern struct chimes_H2_collis_dissoc_struct 
{ 
  int N_reactions[2];               /*!< Number of reactions, excluding and including molecules. */ 
  int **reactants;                  /*!< Reactants for each reaction. */ 
  int **products;                   /*!< Products for each reaction. */ 
  int Heating_reaction_index;       /*!< Index of the collisional dissociation reaction that contributes to heating. */ 
  ChimesFloat **k0;                 /*!< Arrays of the low-density rate coefficients as a function of temperature for each reaction. */ 
  ChimesFloat **kLTE;               /*!< Arrays of the LTE rate coefficients as a function of temperature for each reaction. */ 
  ChimesFloat *critical_density_H;  /*!< Array of the critical density for HI collisions as a function of temperature. */ 
  ChimesFloat *critical_density_H2; /*!< Array of the critical density for H2 collisions as a function of temperature. */ 
  ChimesFloat *critical_density_He; /*!< Array of the critical density for HeI collisions as a function of temperature. */ 
} chimes_table_H2_collis_dissoc; 

/** 
 * Structure containing the rates for 
 * the photoion_fuv reaction group.
 */ 
extern struct chimes_photoion_fuv_struct 
{ 
  int N_reactions[2];        /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;            /*!< Reactant for each reaction. */ 
  int **products;            /*!< Products for each reaction. */ 
  int **element_incl;        /*!< Flags indicating which elements are required for each reaction. */ 
  ChimesFloat *gamma;        /*!< Array of the gamma dust-shielding parameter for each reaction. */ 
  ChimesFloat **sigmaPhot;   /*!< Arrays of the cross sections for each reaction and each spectrum. Units: cm^-2 */ 
  ChimesFloat **epsilonPhot; /*!< Arrays of the average excess photon energy for each reaction and each spectrum. Units: erg */ 
} chimes_table_photoion_fuv; 

/** 
 * Structure containing the rates for 
 * the photoion_euv reaction group.
 */ 
extern struct chimes_photoion_euv_struct 
{ 
  int N_reactions[2];               /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;                   /*!< Reactant for each reaction. */ 
  int **products;                   /*!< Products for each reaction. */ 
  int **element_incl;               /*!< Flags indicating which elements are required for each reaction. */ 
  int *molecular_flag;              /*!< Flags indicating whether each reaction involves molecules. */ 
  ChimesFloat *E_thresh;            /*!< Ionisation energy for each reaction. */ 
  ChimesFloat **sigmaPhot;          /*!< Arrays of the cross sections for each reaction and each spectrum. Units: cm^-2 */ 
  ChimesFloat ****shieldFactor_1D;  /*!< Arrays of the shield factor components in 1D (as a function of 1 column density) for each reaction. */ 
  ChimesFloat *****shieldFactor_2D; /*!< Arrays of the shield factor components in 2D (as a function of 2 column densities) for each reaction. */ 
} chimes_table_photoion_euv; 

/** 
 * Structure containing the rates for 
 * the photoion_auger_fuv reaction group.
 */ 
extern struct chimes_photoion_auger_fuv_struct 
{ 
  int N_reactions[2];       /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;           /*!< Reactant for each reaction. */ 
  int **products;           /*!< Products for each reaction. */ 
  int **element_incl;       /*!< Flags indicating which elements are required for each reaction. */ 
  int *base_reaction;       /*!< Index of the corresponding base reaction, for each auger reaction. */ 
  int *number_of_electrons; /*!< Number of electrons released for each reaction. */ 
  ChimesFloat **sigmaPhot;  /*!< Arrays of the cross sections for each reaction and each spectrum. Units: cm^-2 */ 
} chimes_table_photoion_auger_fuv; 

/** 
 * Structure containing the rates for 
 * the photoion_auger_euv reaction group.
 */ 
extern struct chimes_photoion_auger_euv_struct 
{ 
  int N_reactions[2];       /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;           /*!< Reactant for each reaction. */ 
  int **products;           /*!< Products for each reaction. */ 
  int **element_incl;       /*!< Flags indicating which elements are required for each reaction. */ 
  int *base_reaction;       /*!< Index of the corresponding base reaction, for each auger reaction. */ 
  int *number_of_electrons; /*!< Number of electrons released for each reaction. */ 
  ChimesFloat **sigmaPhot;  /*!< Arrays of the cross sections for each reaction and each spectrum. Units: cm^-2 */ 
} chimes_table_photoion_auger_euv; 

/** 
 * Structure containing the rates for 
 * the photodissoc_group1 reaction group.
 */ 
extern struct chimes_photodissoc_group1_struct 
{ 
  int N_reactions[2];  /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;      /*!< Reactant for each reaction. */ 
  int **products;      /*!< Products for each reaction. */ 
  int **element_incl;  /*!< Flags indicating which elements are required for each reaction. */ 
  int *molecular_flag; /*!< Flags indicating whether each reaction involves molecules. */ 
  ChimesFloat *gamma;  /*!< Array of the gamma dust-shielding parameter for each reaction. */ 
  ChimesFloat *rates;  /*!< Array of the rates for each reaction. */ 
} chimes_table_photodissoc_group1; 

/** 
 * Structure containing the rates for 
 * the photodissoc_group2 reaction group.
 */ 
extern struct chimes_photodissoc_group2_struct 
{ 
  int N_reactions[2];       /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;           /*!< Reactant for each reaction. */ 
  int **products;           /*!< Products for each reaction. */ 
  int **element_incl;       /*!< Flags indicating which elements are required for each reaction. */ 
  ChimesFloat *gamma_coeff; /*!< Parameters to calculate gamma for all reactions in this group. */ 
  ChimesFloat *rates;       /*!< Array of the rates for each reaction. */ 
} chimes_table_photodissoc_group2; 

/** 
 * Structure containing the rates for 
 * the H2_photodissoc reaction group.
 */ 
extern struct chimes_H2_photodissoc_struct 
{ 
  int N_reactions[2];             /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;                 /*!< Reactant for each reaction. */ 
  int **products;                 /*!< Products for each reaction. */ 
  ChimesFloat *gamma;             /*!< Array of the gamma dust-shielding parameter for each reaction. */ 
  ChimesFloat *rates;             /*!< Array of the rates for each reaction. */ 
  ChimesFloat ****self_shielding; /*!< Array of the H2 self-shielding factor as a function of temperature, H2 column density and turbulent doppler broadening parameter. */
} chimes_table_H2_photodissoc; 

/** 
 * Structure containing the rates for 
 * the CO_photodissoc reaction group.
 */ 
extern struct chimes_CO_photodissoc_struct 
{ 
  int N_reactions[2];            /*!< Number of reactions, excluding and including molecules. */ 
  int *reactants;                /*!< Reactant for each reaction. */ 
  int **products;                /*!< Products for each reaction. */ 
  int **element_incl;            /*!< Flags indicating which elements are required for each reaction. */ 
  ChimesFloat *gamma;            /*!< Array of the gamma dust-shielding parameter for each reaction. */ 
  ChimesFloat *rates;            /*!< Array of the rates for each reaction. */ 
  ChimesFloat ***self_shielding; /*!< Array of the CO self-shielding factor as a function CO and H2 column densities. */
} chimes_table_CO_photodissoc; 

/** 
 * Structure containing the 
 * spectra information. 
 */ 
extern struct chimes_spectra_struct 
{ 
  ChimesFloat *isotropic_photon_density; /*!< Isotropic photon density (i.e. flux divided by the speed of light) for each spectrum. Units: cm^-3 */
  ChimesFloat *G0_parameter;             /*!< Strength of the 6-13.6 eV FUV band for each spectrum. Units: Habing. */ 
  ChimesFloat *H2_dissocJ;               /*!< n / (isotropic_photon_density * c) for each spectrum, where n is photon number density in the 12.24-13.51 eV band. */
} chimes_table_spectra; 

/** 
 * Structure containing the cooling 
 * and heating rates. 
 */ 
extern struct chimes_cooling_struct 
{
  int N_coolants;                             /*!< Number of standard cooling channels. */ 
  int N_coolants_2d;                          /*!< Number of 2-d cooling channels. */ 
  int N_coolants_4d;                          /*!< Number of 4-d cooling channels. */ 
  int *coolants;                              /*!< The species corresponding to each standard cooling channel. */ 
  int *coolants_2d;                           /*!< The species corresponding to each 2-d cooling channel. */ 
  int *coolants_4d;                           /*!< The species corresponding to each 4-d cooling channel. */ 
  ChimesFloat **rates;                        /*!< Arrays of the rates for each standard cooling channel as a function of temperature. */ 
  ChimesFloat ***rates_2d;                    /*!< Arrays of the low-temperature rates for each 2-d cooling channel as a function of temperature and ne. */ 
  ChimesFloat *****rates_4d;                  /*!< Arrays of the low-temperature rates for each 4-d cooling channel as a function of temperature, ne, nHI and nHII. */ 
  ChimesFloat **rates_hiT_2d;                 /*!< Arrays of the high-temperature rates for each 2-d cooling channel as a function of temperature. */ 
  ChimesFloat **rates_hiT_4d;                 /*!< Arrays of the high-temperature rates for each 4-d cooling channel as a function of temperature. */ 
  ChimesFloat **photoelectric_heating;        /*!< Photoelectric heating rate as a function of temperature and Psi. */ 
  ChimesFloat *gas_grain_transfer;            /*!< Gas-dust transfer cooling rate as a function of temperature. */ 
  ChimesFloat **grain_recombination;          /*!< Grain recombination cooling rate as a function of temperature and Psi. */ 
  ChimesFloat *H2_cool_lowDens_H2;            /*!< Low-density H2 rovibrational cooling rate from H2 collisions as a function of temperature. */ 
  ChimesFloat *H2_cool_lowDens_HI;            /*!< Low-density H2 rovibrational cooling rate from HI collisions as a function of temperature. */ 
  ChimesFloat *H2_cool_lowDens_HII;           /*!< Low-density H2 rovibrational cooling rate from HII collisions as a function of temperature. */ 
  ChimesFloat *H2_cool_lowDens_HeI;           /*!< Low-density H2 rovibrational cooling rate from HeI collisions as a function of temperature. */ 
  ChimesFloat *H2_cool_lowDens_elec;          /*!< Low-density H2 rovibrational cooling rate from electron collisions as a function of temperature. */ 
  ChimesFloat *H2_cool_LTE;                   /*!< LTE H2 rovibrational cooling rate as a function of temperature. */ 
  ChimesFloat *CO_cool_rot_L0;                /*!< CO rotational cooling L0 as a function of temperature. */ 
  ChimesFloat **CO_cool_rot_Llte;             /*!< CO rotational cooling Llte as a function of temperature and CO column density. */ 
  ChimesFloat **CO_cool_rot_nhalf;            /*!< CO rotational cooling nhalf as a function of temperature and CO column density. */ 
  ChimesFloat **CO_cool_rot_a;                /*!< CO rotational cooling a as a function of temperature and CO column density. */ 
  ChimesFloat *CO_cool_vib_L0;                /*!< CO vibrational cooling L0 as a function of temperature. */ 
  ChimesFloat **CO_cool_vib_Llte;             /*!< CO vibrational cooling Llte as a function of temperature and CO column density. */ 
  ChimesFloat *H2O_cool_rot_hiT_L0;           /*!< H2O rotational cooling high-temperature L0 as a function of temperature. */ 
  ChimesFloat **H2O_cool_rot_hiT_Llte;        /*!< H2O rotational cooling high-temperature Llte as a function of temperature and H2O column density. */ 
  ChimesFloat **H2O_cool_rot_hiT_nhalf;       /*!< H2O rotational cooling high-temperature nhalf as a function of temperature and H2O column density. */ 
  ChimesFloat **H2O_cool_rot_hiT_a;           /*!< H2O rotational cooling high-temperature a as a function of temperature and H2O column density. */ 
  ChimesFloat *H2Oortho_cool_rot_lowT_L0;     /*!< ortho-H2O rotational cooling low-temperature L0 as a function of temperature. */ 
  ChimesFloat **H2Oortho_cool_rot_lowT_Llte;  /*!< ortho-H2O rotational cooling low-temperature Llte as a function of temperature and H2O column density. */ 
  ChimesFloat **H2Oortho_cool_rot_lowT_nhalf; /*!< ortho-H2O rotational cooling low-temperature nhalf as a function of temperature and H2O column density. */ 
  ChimesFloat **H2Oortho_cool_rot_lowT_a;     /*!< ortho-H2O rotational cooling low-temperature a as a function of temperature and H2O column density. */ 
  ChimesFloat *H2Opara_cool_rot_lowT_L0;      /*!< para-H2O rotational cooling low-temperature L0 as a function of temperature. */ 
  ChimesFloat **H2Opara_cool_rot_lowT_Llte;   /*!< para-H2O rotational cooling low-temperature Llte as a function of temperature and H2O column density. */ 
  ChimesFloat **H2Opara_cool_rot_lowT_nhalf;  /*!< para-H2O rotational cooling low-temperature nhalf as a function of temperature and H2O column density. */ 
  ChimesFloat **H2Opara_cool_rot_lowT_a;      /*!< para-H2O rotational cooling low-temperature a as a function of temperature and H2O column density. */ 
  ChimesFloat *H2O_cool_vib_L0;               /*!< H2O vibrational cooling L0 as a function of temperature. */ 
  ChimesFloat **H2O_cool_vib_Llte;            /*!< H2O vibrational cooling Llte as a function of temperature and H2O column density. */ 
} chimes_table_cooling; 

/** 
 * Structure containing the 
 * equilibrium abundance tables. 
 */ 
extern struct chimes_eqm_abundances_struct
{
  int N_Temperatures;          /*!< Number of temperature bins. */ 
  int N_Densities;             /*!< Number of density bins. */ 
  int N_Metallicities;         /*!< Number of metallicity bins. */ 
  ChimesFloat *Temperatures;   /*!< Temperature array. */ 
  ChimesFloat *Densities;      /*!< Density array. */ 
  ChimesFloat *Metallicities;  /*!< Metallicity array. */ 
  ChimesFloat ****Abundances;  /*!< Equilibrium abundance array. */ 
} chimes_table_eqm_abundances; 

/** 
 * Structure containing information needed 
 * to interpolate the redshift-dependent UVB. 
 */ 
extern struct chimes_redshift_dependent_UVB_struct 
{
  int N_redshifts;                                       /*!< Number of redshift bins. */ 
  int z_index_low;                                       /*!< Index of the lower redshift bin that has been loaded. */ 
  int z_index_hi;                                        /*!< Index of the higher redshift bin that has been loaded. */ 
  ChimesFloat *redshift_bins;                            /*!< Redshift array. */ 
  ChimesFloat **photoion_fuv_sigmaPhot;                  /*!< Cross sections for the photoion_fuv group from the two redshift bins. */ 
  ChimesFloat **photoion_fuv_epsilonPhot;                /*!< Average excess photon energies for the photoion_euv group from the two redshift bins. */ 
  ChimesFloat **photoion_euv_sigmaPhot;                  /*!< Cross sections for the photoion_euv group from the two redshift bins. */ 
  ChimesFloat ****photoion_euv_shieldFactor_1D;          /*!< 1D shield factor components for the photoion_euv group from the two redshift bins. */ 
  ChimesFloat *****photoion_euv_shieldFactor_2D;         /*!< 2D shield factor components for the photoion_euv group from the two redshift bins. */ 
  ChimesFloat **photoion_auger_fuv_sigmaPhot;            /*!< Cross sections for the photoion_auger_fuv group from the two redshift bins. */ 
  ChimesFloat **photoion_auger_euv_sigmaPhot;            /*!< Cross sections for the photoion_auger_euv group from the two redshift bins. */ 
  ChimesFloat isotropic_photon_density[2];               /*!< Isotropic photon density from the two redshift bins. */ 
  ChimesFloat G0_parameter[2];                           /*!< Strength of the 6-13.6 FUV band from the two redshift bins. */ 
  ChimesFloat H2_dissocJ[2];                             /*!< n[12.24-13.51 eV] / (isotropic_photon_density * c) from the two redshift bins. */ 
  int **photoion_fuv_element_incl;                       /*!< Flags indicating which elements are required for each reaction in the photoion_fuv group. */ 
  int **photoion_euv_element_incl;                       /*!< Flags indicating which elements are required for each reaction in the photoion_euv group. */ 
  int **photoion_auger_fuv_element_incl;                 /*!< Flags indicating which elements are required for each reaction in the photoion_auger_fuv group. */ 
  int **photoion_auger_euv_element_incl;                 /*!< Flags indicating which elements are required for each reaction in the photoion_auger_euv group. */ 
  struct chimes_eqm_abundances_struct eqm_abundances[2]; /*!< Equilibrium abundance tables from the two redshift bins. */ 
} chimes_table_redshift_dependent_UVB; 


/** 
 * Structure containing the input data 
 * that will be sent to the CVODE solver. 
 */ 
struct UserData
{
  struct gasVariables *myGasVars;                           /*<! The #gasVariables struct. */ 
  struct globalVariables *myGlobalVars;                     /*!< The #globalVariables struct. */ 
  struct Species_Structure *species;                        /*!< The #Species_Structure struct. */ 
  struct chimes_current_rates_struct *chimes_current_rates; /*!< The #chimes_current_rates_struct struct. */ 
  void *cvode_mem;                                          /*!< Pointer to the CVODE memory. */ 
  ChimesFloat HI_column;                                    /*!< HI column density. */ 
  ChimesFloat H2_column;                                    /*!< H2 column density. */ 
  ChimesFloat HeI_column;                                   /*!< HeI column density. */ 
  ChimesFloat HeII_column;                                  /*!< HeII column density. */ 
  ChimesFloat CO_column;                                    /*!< CO column density. */ 
  ChimesFloat H2O_column;                                   /*!< H2O column density. */ 
  ChimesFloat OH_column;                                    /*!< OH column density. */ 
  ChimesFloat extinction;                                   /*!< Dust extinction. */ 
  int network_size;                                         /*!< Size of the network. */ 
  int mol_flag_index;                                       /*!< Flag for whether to exclude (0) or include (1) molecules. */  
  int case_AB_index[2];                                     /*!< Flags to use Case A (0) or Case B (1) recombination for HI and HeI. */ 
};

/** 
 * Structure containing the current rates 
 * and rate coefficients for each reaction. 
 */ 
struct chimes_current_rates_struct
{
  ChimesFloat *data_buffer;                          /*!< Memory buffer for all 1-d arrays. */ 
  ChimesFloat **data_buffer_2d;                      /*!< Memory buffer for all 2-d arrays. */ 
  ChimesFloat *T_dependent_rate_coefficient;         /*!< Current rate coefficients for the T_dependent group. */ 
  ChimesFloat *T_dependent_rate;                     /*!< Current rates for the T_dependent group. */ 
  ChimesFloat *constant_rate;                        /*!< Current rates for the constant group. */ 
  ChimesFloat *recombination_AB_rate_coefficient;    /*!< Current rate coefficients for the recombination_AB group. */ 
  ChimesFloat *recombination_AB_rate;                /*!< Current rates for the recombination_AB group. */ 
  ChimesFloat *grain_recombination_rate_coefficient; /*!< Current rate coefficients for the grain_recombination group. */ 
  ChimesFloat *grain_recombination_rate;             /*!< Current rates for the grain_recombination group. */ 
  ChimesFloat *cosmic_ray_rate;                      /*!< Current rates for the cosmic_ray group. */ 
  ChimesFloat *CO_cosmic_ray_rate_coefficient;       /*!< Current rate coefficients for the CO_cosmic_ray group. */ 
  ChimesFloat *CO_cosmic_ray_rate;                   /*!< Current rates for the CO_cosmic_ray group. */ 
  ChimesFloat H2_dust_formation_rate_coefficient;    /*!< Current rate coefficients for the H2_dust_formation group. */ 
  ChimesFloat H2_dust_formation_rate;                /*!< Current rates for the H2_dust_formation group. */ 
  ChimesFloat *H2_collis_dissoc_rate_coefficient;    /*!< Current rate coefficients for the H2_collis_dissoc group. */ 
  ChimesFloat *H2_collis_dissoc_rate;                /*!< Current rates for the H2_collis_dissoc group. */ 
  ChimesFloat H2_collis_dissoc_crit_H;               /*!< Current critical densities from HI collisions for the H2_collis_dissoc group. */ 
  ChimesFloat H2_collis_dissoc_crit_H2;              /*!< Current critical densities from H2 collisions for the H2_collis_dissoc group. */ 
  ChimesFloat H2_collis_dissoc_crit_He;              /*!< Current critical densities from HeI collisions for the H2_collis_dissoc group. */ 
  ChimesFloat *H2_collis_dissoc_log_k0;              /*!< Current low-density rate coefficients for the H2_collis_dissoc group. */ 
  ChimesFloat *H2_collis_dissoc_log_kLTE;            /*!< Current LTE rate coefficients for the H2_collis_dissoc group. */ 
  ChimesFloat *photoion_fuv_shield_factor;           /*!< Current shield factors for the photoion_fuv group. */ 
  ChimesFloat *photoion_fuv_rate_coefficient;        /*!< Current rate coefficients for the photoion_fuv group. */ 
  ChimesFloat *photoion_fuv_rate;                    /*!< Current rates for the photoion_fuv group. */ 
  ChimesFloat *photoion_fuv_heat_rate;               /*!< Current heating rates for the photoion_fuv group. */ 
  ChimesFloat **photoion_euv_shield_factor;          /*!< Current shield factors for the photoion_euv group. */ 
  ChimesFloat *photoion_euv_rate_coefficient;        /*!< Current rate coefficients for the photoion_euv group. */ 
  ChimesFloat *photoion_euv_rate;                    /*!< Current rates for the photoion_euv group. */ 
  ChimesFloat **photoion_euv_epsilon;                /*!< Current excess photon energies for the photoion_euv group. */ 
  ChimesFloat *photoion_euv_heat_rate;               /*!< Current heating rates for the photoion_euv group. */ 
  ChimesFloat *photoion_auger_fuv_rate_coefficient;  /*!< Current rate coefficients for the photoion_auger_fuv group. */ 
  ChimesFloat *photoion_auger_fuv_rate;              /*!< Current rates for the photoion_auger_fuv group. */ 
  ChimesFloat *photoion_auger_euv_rate_coefficient;  /*!< Current rate coefficients for the photoion_auger_euv group. */ 
  ChimesFloat *photoion_auger_euv_rate;              /*!< Current rates for the photoion_auger_euv group. */ 
  ChimesFloat *photodissoc_group1_shield_factor;     /*!< Current shield factors for the photodissoc_group1 group. */ 
  ChimesFloat *photodissoc_group1_rate_coefficient;  /*!< Current rate coefficients for the photodissoc_group1 group. */ 
  ChimesFloat *photodissoc_group1_rate;              /*!< Current rates for the photodissoc_group1 group. */ 
  ChimesFloat photodissoc_group2_shield_factor;      /*!< Current shield factors for the photodissoc_group2 group. */ 
  ChimesFloat *photodissoc_group2_rate_coefficient;  /*!< Current rate coefficients for the photodissoc_group2 group. */ 
  ChimesFloat *photodissoc_group2_rate;              /*!< Current rates for the photodissoc_group2 group. */ 
  ChimesFloat *H2_photodissoc_shield_factor;         /*!< Current shield factors for the H2_photodissoc group. */ 
  ChimesFloat *H2_photodissoc_rate_coefficient;      /*!< Current rate coefficients for the H2_photodissoc group. */ 
  ChimesFloat *H2_photodissoc_rate;                  /*!< Current rates for the H2_photodissoc group. */ 
  ChimesFloat *CO_photodissoc_shield_factor;         /*!< Current shield factors for the CO_photodissoc group. */ 
  ChimesFloat *CO_photodissoc_rate_coefficient;      /*!< Current rate coefficients for the CO_photodissoc group. */ 
  ChimesFloat *CO_photodissoc_rate;                  /*!< Current rates for the CO_photodissoc group. */ 
  ChimesFloat *cooling_rate;                         /*!< Current rates for the standard cooling channels. */ 
  ChimesFloat *cooling_rate_2d;                      /*!< Current rates for the 2-d cooling channels. */ 
  ChimesFloat *cooling_rate_4d;                      /*!< Current rates for the 4-d cooling channels. */ 
}; 

/** 
 * Enumeration of the indices for 
 * all species in the network. 
 */ 
enum 
  {
    sp_elec,		/* 0 */
    sp_HI,		/* 1 */
    sp_HII,		/* 2 */
    sp_Hm,		/* 3 */
    sp_HeI,		/* 4 */
    sp_HeII,		/* 5 */
    sp_HeIII,		/* 6 */
    sp_CI,		/* 7 */
    sp_CII,		/* 8 */
    sp_CIII,		/* 9 */
    sp_CIV,		/* 10 */
    sp_CV,		/* 11 */
    sp_CVI,		/* 12 */
    sp_CVII,		/* 13 */
    sp_Cm,		/* 14 */
    sp_NI,		/* 15 */
    sp_NII,		/* 16 */
    sp_NIII,		/* 17 */
    sp_NIV,		/* 18 */
    sp_NV,		/* 19 */
    sp_NVI,		/* 20 */
    sp_NVII,		/* 21 */
    sp_NVIII,		/* 22 */
    sp_OI,		/* 23 */
    sp_OII,		/* 24 */
    sp_OIII,		/* 25 */
    sp_OIV, 		/* 26 */
    sp_OV,		/* 27 */
    sp_OVI,		/* 28 */
    sp_OVII,		/* 29 */
    sp_OVIII,		/* 30 */
    sp_OIX,		/* 31 */
    sp_Om,		/* 32 */
    sp_NeI,		/* 33 */
    sp_NeII,		/* 34 */
    sp_NeIII,		/* 35 */
    sp_NeIV,		/* 36 */
    sp_NeV,		/* 37 */
    sp_NeVI,		/* 38 */
    sp_NeVII,		/* 39 */
    sp_NeVIII,		/* 40 */
    sp_NeIX,		/* 41 */
    sp_NeX,		/* 42 */
    sp_NeXI,		/* 43 */
    sp_MgI,		/* 44 */
    sp_MgII,		/* 45 */
    sp_MgIII,		/* 46 */
    sp_MgIV,		/* 47 */
    sp_MgV,		/* 48 */
    sp_MgVI,		/* 49 */
    sp_MgVII,		/* 50 */
    sp_MgVIII,		/* 51 */
    sp_MgIX,		/* 52 */
    sp_MgX,		/* 53 */
    sp_MgXI,		/* 54 */
    sp_MgXII,		/* 55 */
    sp_MgXIII,		/* 56 */
    sp_SiI,		/* 57 */
    sp_SiII,		/* 58 */	
    sp_SiIII,		/* 59 */
    sp_SiIV,		/* 60 */
    sp_SiV,		/* 61 */
    sp_SiVI,		/* 62 */
    sp_SiVII,		/* 63 */
    sp_SiVIII,		/* 64 */
    sp_SiIX,		/* 65 */
    sp_SiX,		/* 66 */
    sp_SiXI,		/* 67 */
    sp_SiXII,		/* 68 */
    sp_SiXIII,		/* 69 */
    sp_SiXIV,		/* 70 */
    sp_SiXV,		/* 71 */
    sp_SI,		/* 72 */
    sp_SII,		/* 73 */
    sp_SIII,		/* 74 */
    sp_SIV,		/* 75 */
    sp_SV,		/* 76 */
    sp_SVI,		/* 77 */
    sp_SVII,		/* 78 */
    sp_SVIII,		/* 79 */
    sp_SIX,		/* 80 */
    sp_SX,		/* 81 */
    sp_SXI,		/* 82 */
    sp_SXII,		/* 83 */
    sp_SXIII,		/* 84 */
    sp_SXIV,		/* 85 */
    sp_SXV,		/* 86 */
    sp_SXVI,		/* 87 */
    sp_SXVII,		/* 88 */
    sp_CaI,		/* 89 */
    sp_CaII,		/* 90 */
    sp_CaIII,		/* 91 */
    sp_CaIV,		/* 92 */
    sp_CaV,		/* 93 */
    sp_CaVI,		/* 94 */
    sp_CaVII,		/* 95 */
    sp_CaVIII,		/* 96 */
    sp_CaIX,		/* 97 */
    sp_CaX,		/* 98 */
    sp_CaXI,		/* 99 */
    sp_CaXII,		/* 100 */
    sp_CaXIII,		/* 101 */
    sp_CaXIV,		/* 102 */
    sp_CaXV,		/* 103 */
    sp_CaXVI,		/* 104 */
    sp_CaXVII,		/* 105 */
    sp_CaXVIII,	        /* 106 */
    sp_CaXIX,		/* 107 */
    sp_CaXX,		/* 108 */
    sp_CaXXI,		/* 109 */
    sp_FeI,		/* 110 */
    sp_FeII,		/* 111 */
    sp_FeIII,		/* 112 */
    sp_FeIV,		/* 113 */
    sp_FeV,		/* 114 */
    sp_FeVI,		/* 115 */
    sp_FeVII,		/* 116 */
    sp_FeVIII,		/* 117 */
    sp_FeIX,		/* 118 */
    sp_FeX,		/* 119 */
    sp_FeXI,		/* 120 */
    sp_FeXII,		/* 121 */
    sp_FeXIII,		/* 122 */
    sp_FeXIV,		/* 123 */
    sp_FeXV,		/* 124 */
    sp_FeXVI,		/* 125 */
    sp_FeXVII,		/* 126 */
    sp_FeXVIII,	        /* 127 */
    sp_FeXIX,		/* 128 */
    sp_FeXX,		/* 129 */
    sp_FeXXI,		/* 130 */
    sp_FeXXII,		/* 131 */
    sp_FeXXIII,	        /* 132 */
    sp_FeXXIV,		/* 133 */
    sp_FeXXV,		/* 134 */
    sp_FeXXVI,		/* 135 */
    sp_FeXXVII,	        /* 136 */
    sp_H2,		/* 137 */
    sp_H2p,		/* 138 */
    sp_H3p,		/* 139 */
    sp_OH,		/* 140 */
    sp_H2O,		/* 141 */
    sp_C2,		/* 142 */
    sp_O2,		/* 143 */
    sp_HCOp,		/* 144 */
    sp_CH,		/* 145 */
    sp_CH2,		/* 146 */
    sp_CH3p,		/* 147 */
    sp_CO,		/* 148 */
    sp_CHp,		/* 149 */
    sp_CH2p,		/* 150 */
    sp_OHp,		/* 151 */
    sp_H2Op,		/* 152 */
    sp_H3Op,		/* 153 */
    sp_COp,		/* 154 */
    sp_HOCp,		/* 155 */
    sp_O2p		/* 156 */
  };
