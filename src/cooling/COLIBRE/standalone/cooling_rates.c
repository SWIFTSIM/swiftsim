#include <stdio.h>
#include <math.h>

#include "physical_constants_cgs.h"
#include "cooling_tables.h"
#include "interpolate.h"


/**
 * @brief Computes the net cooling rate (cooling - heating) for a given element abundance ratio,
 * internal energy, redshift, and density. The unit of the net cooling rate is 
 * Lambda / nH**2 [erg cm^3 s-1] and all input values are in cgs. 
 * The Compton cooling is not taken from the tables but calculated analytically and added separately
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
 */

double colibre_cooling_rate(
    double log_u_cgs, double redshift, double n_H_cgs, float ZZsol,
    const float abundance_ratio[colibre_cooling_N_elementtypes], int n_H_index,
    float d_n_H, int met_index, float d_met, int red_index, float d_red,
    const struct cooling_function_data *restrict cooling){

  /* Get index of u along the internal energy axis */
    int U_index;
    float d_U;

    double cooling_rate, heating_rate, Compton_cooling_rate, temp, logtemp;
    double net_cooling_rate, electron_fraction;

    float weights_cooling[colibre_cooling_N_cooltypes-2];
    float weights_heating[colibre_cooling_N_heattypes-2];

    int i;

    /* set weights for cooling rates */
    for (i = 0; i < colibre_cooling_N_cooltypes-2; i++) {
  	if ( i <= colibre_cooling_N_elementtypes ) weights_cooling[i] = abundance_ratio[i];
	if ( i == cooltype_H2 )                    weights_cooling[i] = 1.;  /* use same H2 abundance as in tables */
	if ( i == cooltype_molecules )             weights_cooling[i] = 1.; 
	if ( i == cooltype_HD )                    weights_cooling[i] = 1.;  /* use same HD abundance as in tables */
	if ( i == cooltype_NetFFH )                weights_cooling[i] = 1.; 
	if ( i == cooltype_NetFFM )                weights_cooling[i] = 1.; 
	if ( i == cooltype_eeBrems )               weights_cooling[i] = 1.;  /* use same electron abundance */
	if ( i == cooltype_Compton )               weights_cooling[i] = 0.;  /* added analytically */ 
	if ( i == cooltype_Dust )                  weights_cooling[i] = 1.; 
    }
   
    /* set weights for heating rates */
    for (i = 0; i < colibre_cooling_N_heattypes-2; i++) {
  	if ( i <= colibre_cooling_N_elementtypes ) weights_heating[i] = abundance_ratio[i];
	if ( i == heattype_H2 )                    weights_heating[i] = 1.; 
	if ( i == heattype_COdiss )                weights_heating[i] = 1.; 
	if ( i == heattype_CosmicRay )             weights_heating[i] = 1.; 
	if ( i == heattype_UTA )                   weights_heating[i] = 1.; 
	if ( i == heattype_line )                  weights_heating[i] = 1.; 
	if ( i == heattype_Hlin )                  weights_heating[i] = 1.; 
	if ( i == heattype_ChaT )                  weights_heating[i] = 1.; 
	if ( i == heattype_HFF  )                  weights_heating[i] = 1.; 
	if ( i == heattype_Compton )               weights_heating[i] = 1.; 
	if ( i == heattype_Dust )                  weights_heating[i] = 1.; 
    }



    get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_u_cgs, &U_index, &d_U);
        
    /* n_e / n_H */
    electron_fraction = interpolation4d_plus_summation(cooling->table.Uelectron_fraction, abundance_ratio,
		   element_H, colibre_cooling_N_electrontypes - 4, 
                   red_index, U_index, met_index, n_H_index,
                   d_red, d_U, d_met, d_n_H,
                   colibre_cooling_N_redshifts, colibre_cooling_N_internalenergy,
                   colibre_cooling_N_metallicity, colibre_cooling_N_density, colibre_cooling_N_electrontypes);

    /* Lambda / n_H**2 */
    cooling_rate = interpolation4d_plus_summation(cooling->table.Ucooling, weights_cooling,
                   element_H, colibre_cooling_N_cooltypes - 3, 
                   red_index, U_index, met_index, n_H_index,
                   d_red, d_U, d_met, d_n_H, 
                   colibre_cooling_N_redshifts, colibre_cooling_N_internalenergy, 
                   colibre_cooling_N_metallicity, colibre_cooling_N_density, colibre_cooling_N_cooltypes);

    /* Gamma / n_H**2 */
    heating_rate = interpolation4d_plus_summation(cooling->table.Uheating, weights_heating,
                   element_H, colibre_cooling_N_heattypes - 3,
                   red_index, U_index, met_index, n_H_index,
                   d_red, d_U, d_met, d_n_H,
                   colibre_cooling_N_redshifts, colibre_cooling_N_internalenergy,
                   colibre_cooling_N_metallicity, colibre_cooling_N_density, colibre_cooling_N_heattypes);

    /* Temperature from internal energy */
    logtemp = interpolation4d(cooling->table.T_from_U, 
                   red_index, U_index, met_index, n_H_index,
                   d_red, d_U, d_met, d_n_H,
                   colibre_cooling_N_redshifts, colibre_cooling_N_internalenergy,
                   colibre_cooling_N_metallicity, colibre_cooling_N_density);  

    temp = pow(10., logtemp);

    /* Analytic Compton cooling rate: Lambda_Compton / n_H**2 */
    Compton_cooling_rate = cooling->Compton_constant_cgs * electron_fraction / n_H_cgs * 
                           pow(const_T_CMB_0_cgs * (1. + redshift), 4.) * (temp - const_T_CMB_0_cgs * (1. + redshift));

    printf("%.4f\t%.4f\t%.4e\t%.4e\t%.4e\t%.4e\n", log_u_cgs, logtemp, cooling_rate, heating_rate, Compton_cooling_rate, 
                                             electron_fraction);
    
    net_cooling_rate = cooling_rate + Compton_cooling_rate - heating_rate;

    return net_cooling_rate;  
}

