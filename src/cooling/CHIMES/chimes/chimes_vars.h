#include <stdio.h>
#include <stdlib.h>

#define ELECTRON_MASS 9.1093829e-28
#define PROTON_MASS 1.6726218e-24
#define PI            3.1415927
#define LIGHTSPEED    3.0e10		/* In cm/s */
#define BOLTZMANNCGS     1.3806e-16	/* In ergs/K */
#define BOLTZMANN_EVK   8.6173324e-5	/*In eV/K*/
#define G0_GAMMA 2.77  /* For dust processes involving G0, e.g. photoelectric heating, we attenuate
			* G0 by exp(- G0_gamma * Av) */

#define TOTSIZE	  157    /* The total number of species in the full network */
#define MOLSIZE   20     /* Total number of molecules in the full network */
#define DUSTEFFSIZE 4.0e-22	/* in cm^2 / H atom. Used to convert between Av & N_Htot */
#define DUST_CROSS_SECTION 1.0e-21  /* in cm^2 / H atom. Used in H2 formation on dust. 
				     * NOTE: NOT the same as DUSTEFFSIZE. The latter is
				     * Av / NH and thus includes absorption AND scattering. */
#define MAX_TEMPERATURE_FOR_RATES  2.0e9   /* If the temperature in gasVars exceeds this value, 
                                            * we limit the temperature used to calculate reaction 
					    * and cooling rates (in set_rates.c and cooling.c) to 
                                            * this value. */ 
#define MAX_UV_SPECTRA 20 

#define METALS_MINIMUM_THRESHOLD 1.0e-30 

#define MAXSTEPS 1e5      /* The maximum number of steps in the CVODE solver. */ 

typedef double ChimesFloat; 

/* This structure contains variables that 
 * are specific to each gas particle/cell. */ 
struct gasVariables
{
  /* NOTE: all abundances are in
   * the form ni / nHtot, i.e. the
   * ratio of number densities of
   * species i to total hydrogen. */
  ChimesFloat element_abundances[10];  /* In the order He, C, N, O, Ne, Mg, Si, S, Ca, Fe */
  ChimesFloat nH_tot;                  /* cm^-3 */
  ChimesFloat temperature;             /* K */
  ChimesFloat TempFloor;
  ChimesFloat divVel;                  /* s^-1 */
  ChimesFloat doppler_broad;           /* km s^-1. NOTE: this is ONLY from turbulence; thermal broadening is added later. */
  ChimesFloat *isotropic_photon_density;  /* cm^-3 */
  ChimesFloat *G0_parameter;
  ChimesFloat *H2_dissocJ;             /* This is n / (isotropic_photon_density * c),
				   * where n is photon number density in the 
				   * band 12.24 eV to 13.51 eV (all in cgs units). */
  ChimesFloat cr_rate;
  ChimesFloat metallicity;             /* Z / Z_sol */
  ChimesFloat dust_ratio;              /* Relative to Milky Way */
  ChimesFloat cell_size;               /* cm; use kernel smoothing length in SPH */
  ChimesFloat hydro_timestep;          /* s */
  int ForceEqOn;
  int ThermEvolOn;
  ChimesFloat constant_heating_rate;   /* erg s^-1 cm^-3 */
  ChimesFloat *abundances;             /* The size of this array will be set by init_chimes() */
};

/* This structure contains the global 
 * variables */ 
struct globalVariables
{
  /* The following are parameters that will
   * need to be read in from the parameter file. */
  char MainDataTablePath[500];
  char PhotoIonTablePath[MAX_UV_SPECTRA][500];
  char EqAbundanceTablePath[500];
  int cellSelfShieldingOn;
  int N_spectra;          /* The number of UV spectra. */ 
  int StaticMolCooling;
  ChimesFloat T_mol;
  int InitIonState;
  ChimesFloat grain_temperature;
  ChimesFloat cmb_temperature;
  ChimesFloat relativeTolerance;
  ChimesFloat absoluteTolerance;
  ChimesFloat thermalAbsoluteTolerance;
  ChimesFloat explicitTolerance; 
  int element_included[9];   /* Metals only */
  int speciesIndices[TOTSIZE];
  int totalNumberOfSpecies;
  int scale_metal_tolerances; 
  int update_colibre_ISRF; 
  int update_colibre_shielding; 
  ChimesFloat max_shielding_length_cgs; 
  ChimesFloat shielding_length_factor; 
  ChimesFloat colibre_cr_rate_0; 
}; 

/* The following structure contains
 * information about each non-eq species */
struct Species_Structure 
{
  int include_species; 
  ChimesFloat element_abundance; 
  ChimesFloat creation_rate; 
  ChimesFloat destruction_rate; 
};

/* The following data structures will contain 
 * the various data tables from chimes_main_data.hdf5. */ 
extern struct chimes_table_bins_struct 
{ 
  int N_Temperatures; 
  ChimesFloat *Temperatures; 
  int N_Dust_Temperatures; 
  ChimesFloat *Dust_Temperatures; 
  int N_Psi; 
  ChimesFloat *Psi; 
  int N_secondary_cosmic_ray_xHII; 
  ChimesFloat *secondary_cosmic_ray_xHII; 
  int N_Column_densities; 
  ChimesFloat *Column_densities; 
  int N_H2self_column_densities; 
  ChimesFloat *H2self_column_densities; 
  int N_b_turbulence; 
  ChimesFloat *b_turbulence; 
  int N_COself_column_densities; 
  ChimesFloat *COself_column_densities; 
  int N_H2CO_column_densities; 
  ChimesFloat *H2CO_column_densities; 
  int N_mol_cool_Temperatures; 
  ChimesFloat *mol_cool_Temperatures; 
  int N_CO_cool_rot_ColumnDensities; 
  ChimesFloat *CO_cool_rot_ColumnDensities; 
  int N_CO_cool_vib_ColumnDensities; 
  ChimesFloat *CO_cool_vib_ColumnDensities; 
  int N_H2O_cool_hiT_Temperatures; 
  ChimesFloat *H2O_cool_hiT_Temperatures; 
  int N_H2O_cool_lowT_Temperatures; 
  ChimesFloat *H2O_cool_lowT_Temperatures; 
  int N_H2O_cool_rot_ColumnDensities; 
  ChimesFloat *H2O_cool_rot_ColumnDensities; 
  int N_H2O_cool_vib_ColumnDensities; 
  ChimesFloat *H2O_cool_vib_ColumnDensities; 
  int N_cool_2d_Temperatures; 
  ChimesFloat *cool_2d_Temperatures; 
  int N_cool_hiT_2d_Temperatures; 
  ChimesFloat *cool_hiT_2d_Temperatures; 
  int N_cool_2d_ElectronDensities; 
  ChimesFloat *cool_2d_ElectronDensities; 
  int N_cool_4d_Temperatures; 
  ChimesFloat *cool_4d_Temperatures; 
  int N_cool_hiT_4d_Temperatures; 
  ChimesFloat *cool_hiT_4d_Temperatures; 
  int N_cool_4d_HIDensities; 
  ChimesFloat *cool_4d_HIDensities; 
  int N_cool_4d_ElectronDensities; 
  ChimesFloat *cool_4d_ElectronDensities; 
  int N_cool_4d_HIIDensities; 
  ChimesFloat *cool_4d_HIIDensities; 
} chimes_table_bins; 

extern struct chimes_T_dependent_struct 
{ 
  int N_reactions[2]; 
  int **reactants; 
  int **products; 
  int **element_incl; 
  int *molecular_flag; 
  int H2_collis_dissoc_heating_reaction_index; 
  int H2_form_heating_reaction_index; 
  ChimesFloat **rates; 
} chimes_table_T_dependent; 

extern struct chimes_constant_struct 
{ 
  int N_reactions[2]; 
  int **reactants; 
  int **products; 
  int **element_incl; 
  int *molecular_flag;
  int H2_form_heating_reaction_index; 
  ChimesFloat *rates; 
} chimes_table_constant; 

extern struct chimes_recombination_AB_struct 
{ 
  int N_reactions[2]; 
  int **reactants; 
  int *products; 
  int **element_incl; 
  int *molecular_flag; 
  ChimesFloat ***rates; 
} chimes_table_recombination_AB; 

extern struct chimes_grain_recombination_struct 
{ 
  int N_reactions[2]; 
  int **reactants; 
  int *products; 
  int **element_incl; 
  ChimesFloat ***rates; 
} chimes_table_grain_recombination; 

extern struct chimes_cosmic_ray_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  int *molecular_flag;
  int *secondary_base_reaction; 
  ChimesFloat **secondary_ratio; 
  ChimesFloat *rates; 
} chimes_table_cosmic_ray; 

extern struct chimes_CO_cosmic_ray_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  ChimesFloat **rates; 
} chimes_table_CO_cosmic_ray; 

extern struct chimes_H2_dust_formation_struct 
{ 
  int *reactants; 
  int *products; 
  ChimesFloat **rates; 
} chimes_table_H2_dust_formation; 

extern struct chimes_H2_collis_dissoc_struct 
{ 
  int N_reactions[2]; 
  int **reactants; 
  int **products; 
  int Heating_reaction_index; 
  ChimesFloat **k0; 
  ChimesFloat **kLTE; 
  ChimesFloat *critical_density_H; 
  ChimesFloat *critical_density_H2; 
  ChimesFloat *critical_density_He; 
} chimes_table_H2_collis_dissoc; 

extern struct chimes_photoion_fuv_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  ChimesFloat *gamma;   
  ChimesFloat **sigmaPhot;                 // cm^-2 
  ChimesFloat **epsilonPhot;                // erg 
} chimes_table_photoion_fuv; 

extern struct chimes_photoion_euv_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  int *molecular_flag; 
  ChimesFloat *E_thresh; 
  ChimesFloat **sigmaPhot;                  // cm^-2 
  ChimesFloat ****shieldFactor_1D; 
  ChimesFloat *****shieldFactor_2D; 
} chimes_table_photoion_euv; 

extern struct chimes_photoion_auger_fuv_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  int *base_reaction; 
  int *number_of_electrons; 
  ChimesFloat **sigmaPhot;                  // cm^-2 
} chimes_table_photoion_auger_fuv; 

extern struct chimes_photoion_auger_euv_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  int *base_reaction; 
  int *number_of_electrons; 
  ChimesFloat **sigmaPhot;                  // cm^-2 
} chimes_table_photoion_auger_euv; 

extern struct chimes_photodissoc_group1_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  int *molecular_flag; 
  ChimesFloat *gamma;  
  ChimesFloat *rates; 
} chimes_table_photodissoc_group1; 

extern struct chimes_photodissoc_group2_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  ChimesFloat *gamma_coeff;  
  ChimesFloat *rates; 
} chimes_table_photodissoc_group2; 

extern struct chimes_H2_photodissoc_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  ChimesFloat *gamma;  
  ChimesFloat *rates; 
  ChimesFloat ****self_shielding; 
} chimes_table_H2_photodissoc; 

extern struct chimes_CO_photodissoc_struct 
{ 
  int N_reactions[2]; 
  int *reactants; 
  int **products; 
  int **element_incl; 
  ChimesFloat *gamma;  
  ChimesFloat *rates; 
  ChimesFloat ***self_shielding; 
} chimes_table_CO_photodissoc; 

extern struct chimes_spectra_struct 
{ 
  ChimesFloat *isotropic_photon_density;
  ChimesFloat *G0_parameter; 
  ChimesFloat *H2_dissocJ; 
} chimes_table_spectra; 

extern struct chimes_cooling_struct 
{
  int N_coolants; 
  int N_coolants_2d; 
  int N_coolants_4d; 
  int *coolants; 
  int *coolants_2d; 
  int *coolants_4d; 
  ChimesFloat **rates; 
  ChimesFloat ***rates_2d; 
  ChimesFloat *****rates_4d; 
  ChimesFloat **rates_hiT_2d; 
  ChimesFloat **rates_hiT_4d; 
  ChimesFloat **photoelectric_heating; 
  ChimesFloat *gas_grain_transfer; 
  ChimesFloat **grain_recombination; 
  ChimesFloat *H2_cool_lowDens_H2; 
  ChimesFloat *H2_cool_lowDens_HI; 
  ChimesFloat *H2_cool_lowDens_HII; 
  ChimesFloat *H2_cool_lowDens_HeI; 
  ChimesFloat *H2_cool_lowDens_elec; 
  ChimesFloat *H2_cool_LTE; 
  ChimesFloat *CO_cool_rot_L0; 
  ChimesFloat **CO_cool_rot_Llte; 
  ChimesFloat **CO_cool_rot_nhalf; 
  ChimesFloat **CO_cool_rot_a; 
  ChimesFloat *CO_cool_vib_L0; 
  ChimesFloat **CO_cool_vib_Llte; 
  ChimesFloat *H2O_cool_rot_hiT_L0; 
  ChimesFloat **H2O_cool_rot_hiT_Llte; 
  ChimesFloat **H2O_cool_rot_hiT_nhalf; 
  ChimesFloat **H2O_cool_rot_hiT_a; 
  ChimesFloat *H2Oortho_cool_rot_lowT_L0; 
  ChimesFloat **H2Oortho_cool_rot_lowT_Llte; 
  ChimesFloat **H2Oortho_cool_rot_lowT_nhalf; 
  ChimesFloat **H2Oortho_cool_rot_lowT_a; 
  ChimesFloat *H2Opara_cool_rot_lowT_L0; 
  ChimesFloat **H2Opara_cool_rot_lowT_Llte; 
  ChimesFloat **H2Opara_cool_rot_lowT_nhalf; 
  ChimesFloat **H2Opara_cool_rot_lowT_a; 
  ChimesFloat *H2O_cool_vib_L0; 
  ChimesFloat **H2O_cool_vib_Llte; 
} chimes_table_cooling; 

extern struct chimes_eqm_abundances_struct
{
  int N_Temperatures;
  int N_Densities;
  int N_Metallicities;
  ChimesFloat *Temperatures;
  ChimesFloat *Densities;
  ChimesFloat *Metallicities;
  ChimesFloat ****Abundances;
} chimes_table_eqm_abundances; 

struct UserData
{
  struct gasVariables *myGasVars;
  struct globalVariables *myGlobalVars;
  struct Species_Structure *species;
  struct chimes_current_rates_struct *chimes_current_rates; 
  void *cvode_mem;
  ChimesFloat HI_column;
  ChimesFloat H2_column;
  ChimesFloat HeI_column;
  ChimesFloat HeII_column;
  ChimesFloat CO_column;
  ChimesFloat H2O_column;
  ChimesFloat OH_column;
  ChimesFloat extinction;
  int network_size;
  int mol_flag_index; 
  int case_AB_index[2]; 
};

struct chimes_current_rates_struct
{
  ChimesFloat *data_buffer; 
  ChimesFloat **data_buffer_2d; 
  ChimesFloat *T_dependent_rate_coefficient; 
  ChimesFloat *T_dependent_rate; 
  ChimesFloat *constant_rate; 
  ChimesFloat *recombination_AB_rate_coefficient; 
  ChimesFloat *recombination_AB_rate; 
  ChimesFloat *grain_recombination_rate_coefficient; 
  ChimesFloat *grain_recombination_rate; 
  ChimesFloat *cosmic_ray_rate; 
  ChimesFloat *CO_cosmic_ray_rate_coefficient; 
  ChimesFloat *CO_cosmic_ray_rate; 
  ChimesFloat H2_dust_formation_rate_coefficient;
  ChimesFloat H2_dust_formation_rate; 
  ChimesFloat *H2_collis_dissoc_rate_coefficient; 
  ChimesFloat *H2_collis_dissoc_rate; 
  ChimesFloat H2_collis_dissoc_crit_H; 
  ChimesFloat H2_collis_dissoc_crit_H2; 
  ChimesFloat H2_collis_dissoc_crit_He; 
  ChimesFloat *H2_collis_dissoc_log_k0; 
  ChimesFloat *H2_collis_dissoc_log_kLTE; 
  ChimesFloat *photoion_fuv_shield_factor; 
  ChimesFloat *photoion_fuv_rate_coefficient; 
  ChimesFloat *photoion_fuv_rate; 
  ChimesFloat *photoion_fuv_heat_rate; 
  ChimesFloat **photoion_euv_shield_factor; 
  ChimesFloat *photoion_euv_rate_coefficient; 
  ChimesFloat *photoion_euv_rate; 
  ChimesFloat **photoion_euv_epsilon; 
  ChimesFloat *photoion_euv_heat_rate; 
  ChimesFloat *photoion_auger_fuv_rate_coefficient; 
  ChimesFloat *photoion_auger_fuv_rate; 
  ChimesFloat *photoion_auger_euv_rate_coefficient; 
  ChimesFloat *photoion_auger_euv_rate; 
  ChimesFloat *photodissoc_group1_shield_factor; 
  ChimesFloat *photodissoc_group1_rate_coefficient; 
  ChimesFloat *photodissoc_group1_rate;
  ChimesFloat photodissoc_group2_shield_factor; 
  ChimesFloat *photodissoc_group2_rate_coefficient; 
  ChimesFloat *photodissoc_group2_rate; 
  ChimesFloat *H2_photodissoc_shield_factor; 
  ChimesFloat *H2_photodissoc_rate_coefficient; 
  ChimesFloat *H2_photodissoc_rate; 
  ChimesFloat *CO_photodissoc_shield_factor; 
  ChimesFloat *CO_photodissoc_rate_coefficient; 
  ChimesFloat *CO_photodissoc_rate; 
  ChimesFloat *cooling_rate; 
  ChimesFloat *cooling_rate_2d; 
  ChimesFloat *cooling_rate_4d; 
}; 

enum 
  {
    elec,		/* 0 */
    HI,		/* 1 */
    HII,		/* 2 */
    Hm,		/* 3 */
    HeI,		/* 4 */
    HeII,		/* 5 */
    HeIII,		/* 6 */
    CI,		/* 7 */
    CII,		/* 8 */
    CIII,		/* 9 */
    CIV,		/* 10 */
    CV,		/* 11 */
    CVI,		/* 12 */
    CVII,		/* 13 */
    Cm,		/* 14 */
    NI,		/* 15 */
    NII,		/* 16 */
    NIII,		/* 17 */
    NIV,		/* 18 */
    NV,		/* 19 */
    NVI,		/* 20 */
    NVII,		/* 21 */
    NVIII,		/* 22 */
    OI,		/* 23 */
    OII,		/* 24 */
    OIII,		/* 25 */
    OIV, 		/* 26 */
    OV,		/* 27 */
    OVI,		/* 28 */
    OVII,		/* 29 */
    OVIII,		/* 30 */
    OIX,		/* 31 */
    Om,		/* 32 */
    NeI,		/* 33 */
    NeII,		/* 34 */
    NeIII,		/* 35 */
    NeIV,		/* 36 */
    NeV,		/* 37 */
    NeVI,		/* 38 */
    NeVII,		/* 39 */
    NeVIII,		/* 40 */
    NeIX,		/* 41 */
    NeX,		/* 42 */
    NeXI,		/* 43 */
    MgI,		/* 44 */
    MgII,		/* 45 */
    MgIII,		/* 46 */
    MgIV,		/* 47 */
    MgV,		/* 48 */
    MgVI,		/* 49 */
    MgVII,		/* 50 */
    MgVIII,		/* 51 */
    MgIX,		/* 52 */
    MgX,		/* 53 */
    MgXI,		/* 54 */
    MgXII,		/* 55 */
    MgXIII,		/* 56 */
    SiI,		/* 57 */
    SiII,		/* 58 */	
    SiIII,		/* 59 */
    SiIV,		/* 60 */
    SiV,		/* 61 */
    SiVI,		/* 62 */
    SiVII,		/* 63 */
    SiVIII,		/* 64 */
    SiIX,		/* 65 */
    SiX,		/* 66 */
    SiXI,		/* 67 */
    SiXII,		/* 68 */
    SiXIII,		/* 69 */
    SiXIV,		/* 70 */
    SiXV,		/* 71 */
    SI,		/* 72 */
    SII,		/* 73 */
    SIII,		/* 74 */
    SIV,		/* 75 */
    SV,		/* 76 */
    SVI,		/* 77 */
    SVII,		/* 78 */
    SVIII,		/* 79 */
    SIX,		/* 80 */
    SX,		/* 81 */
    SXI,		/* 82 */
    SXII,		/* 83 */
    SXIII,		/* 84 */
    SXIV,		/* 85 */
    SXV,		/* 86 */
    SXVI,		/* 87 */
    SXVII,		/* 88 */
    CaI,		/* 89 */
    CaII,		/* 90 */
    CaIII,		/* 91 */
    CaIV,		/* 92 */
    CaV,		/* 93 */
    CaVI,		/* 94 */
    CaVII,		/* 95 */
    CaVIII,		/* 96 */
    CaIX,		/* 97 */
    CaX,		/* 98 */
    CaXI,		/* 99 */
    CaXII,		/* 100 */
    CaXIII,		/* 101 */
    CaXIV,		/* 102 */
    CaXV,		/* 103 */
    CaXVI,		/* 104 */
    CaXVII,		/* 105 */
    CaXVIII,	/* 106 */
    CaXIX,		/* 107 */
    CaXX,		/* 108 */
    CaXXI,		/* 109 */
    FeI,		/* 110 */
    FeII,		/* 111 */
    FeIII,		/* 112 */
    FeIV,		/* 113 */
    FeV,		/* 114 */
    FeVI,		/* 115 */
    FeVII,		/* 116 */
    FeVIII,		/* 117 */
    FeIX,		/* 118 */
    FeX,		/* 119 */
    FeXI,		/* 120 */
    FeXII,		/* 121 */
    FeXIII,		/* 122 */
    FeXIV,		/* 123 */
    FeXV,		/* 124 */
    FeXVI,		/* 125 */
    FeXVII,		/* 126 */
    FeXVIII,	/* 127 */
    FeXIX,		/* 128 */
    FeXX,		/* 129 */
    FeXXI,		/* 130 */
    FeXXII,		/* 131 */
    FeXXIII,	/* 132 */
    FeXXIV,		/* 133 */
    FeXXV,		/* 134 */
    FeXXVI,		/* 135 */
    FeXXVII,	/* 136 */
    H2,		/* 137 */
    H2p,		/* 138 */
    H3p,		/* 139 */
    OH,		/* 140 */
    H2O,		/* 141 */
    C2,		/* 142 */
    O2,		/* 143 */
    HCOp,		/* 144 */
    CH,		/* 145 */
    CH2,		/* 146 */
    CH3p,		/* 147 */
    CO,		/* 148 */
    CHp,		/* 149 */
    CH2p,		/* 150 */
    OHp,		/* 151 */
    H2Op,		/* 152 */
    H3Op,		/* 153 */
    COp,		/* 154 */
    HOCp,		/* 155 */
    O2p		/* 156 */
  };
