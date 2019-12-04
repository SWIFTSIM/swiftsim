#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h> 
#include <sundials/sundials_dense.h>
#include <hdf5.h>

#ifndef chimes_max
#define chimes_max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef chimes_min
#define chimes_min(a,b) ((a) < (b) ? (a) : (b))
#endif

// chimes_cooling.c 
ChimesFloat calculate_mean_molecular_weight(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars); 
ChimesFloat calculate_total_cooling_rate(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data, int mode); 
ChimesFloat calculate_total_number_density(ChimesFloat *my_abundances, ChimesFloat nH, struct globalVariables *myGlobalVars); 
ChimesFloat compton_cooling(ChimesFloat T, ChimesFloat Tcmb, ChimesFloat xe, ChimesFloat nH); 
void do_equilibrium_cooling(struct UserData data); 
ChimesFloat OH_rotational_cooling(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data); 
void update_cooling_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data); 

// chimes.c 
void chimes_network(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void set_equilibrium_abundances_from_tables(struct UserData data);

// init_chimes.c 
void allocate_gas_abundances_memory(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void free_gas_abundances_memory(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars); 
int compare_element_incl_arrays(int *reaction_array, int *network_array);
void allocate_eqm_table_memory(char *filename, struct chimes_eqm_abundances_struct *my_eqm_abundances, struct globalVariables *myGlobalVars); 
void load_eqm_table(char *filename, struct chimes_eqm_abundances_struct *my_eqm_abundances, struct globalVariables *myGlobalVars); 
void init_chimes(struct globalVariables *myGlobalVars); 
void initialise_gas_abundances(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void initialise_main_data(struct globalVariables *myGlobalVars);  
void read_cross_sections_tables(struct chimes_table_bins_struct *my_table_bins, struct chimes_photoion_fuv_struct *my_photoion_fuv, struct chimes_photoion_euv_struct *my_photoion_euv, struct chimes_photoion_auger_fuv_struct *my_photoion_auger_fuv, struct chimes_photoion_auger_euv_struct *my_photoion_auger_euv, struct chimes_spectra_struct *my_spectra, struct globalVariables *myGlobalVars); 
int set_species_index_array(struct globalVariables *myGlobalVariables);
void determine_current_rates_buffer_size(int *buffer_size, int *buffer_size_2D, struct globalVariables *myGlobalVars); 
void allocate_current_rates_memory(struct chimes_current_rates_struct *chimes_current_rates, struct globalVariables *myGlobalVars); 
void free_current_rates_memory(struct chimes_current_rates_struct *chimes_current_rates, struct globalVariables *myGlobalVars); 
void allocate_redshift_dependent_UVB_memory(struct globalVariables *myGlobalVars); 
void load_redshift_dependent_UVB(ChimesFloat redshift, int bin_index, struct globalVariables *myGlobalVars); 

// interpol.c 
void chimes_get_table_index(ChimesFloat *table, int ntable, ChimesFloat x, int *i, ChimesFloat *dx); 
ChimesFloat chimes_interpol_1d(ChimesFloat *table, int i, ChimesFloat dx); 
ChimesFloat chimes_interpol_2d(ChimesFloat **table, int i, int j, ChimesFloat dx, ChimesFloat dy); 
ChimesFloat chimes_interpol_3d(ChimesFloat ***table, int i, int j, int k, ChimesFloat dx, ChimesFloat dy, ChimesFloat dz); 
ChimesFloat chimes_interpol_4d(ChimesFloat ****table, int i, int j, int k, int l, ChimesFloat dx, ChimesFloat dy, ChimesFloat dz, ChimesFloat dw); 

// rate_equations.c 
void check_constraint_equations(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// update_rates.c 
void set_initial_rate_coefficients(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data); 
void set_species_structures(struct Species_Structure *mySpecies, struct gasVariables *myGasVars, int *total_network, int *nonmolecular_network, struct globalVariables *myGlobalVars);
void update_rate_coefficients(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data, int mode); 
void update_rate_vector(struct Species_Structure *mySpecies, struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data); 
void update_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data); 
void redshift_dependent_UVB_copy_lowz_to_hiz(struct globalVariables *myGlobalVars); 
void interpolate_redshift_dependent_UVB(struct globalVariables *myGlobalVars); 

