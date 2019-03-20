#include <math.h>
#include <stdio.h>

#include "interpolate.h"
#include "cooling_rates.h"

/**
 * @brief Uses the cooling tables to print the cooling function for 
 * a given redshift, hydrogen number density (n_H_cgs), and element abundance ratios
 * (set as particle metal masses, part->metal_mass)
 *
 * @param cooling Cooling data structure
 * @param part Particle data structure
 */

void cooling_function(struct cooling_function_data *cooling, struct particle_data *part) {

  /* ----- Cooling function is calculated for these parameter ------*/
  double redshift  = 0.0;
  double n_H_cgs   = 1.e-4;
  double log_u_cgs_min = 7.;
  double log_u_cgs_max = 18;
  double dlog_u_cgs    = 0.2;
  /* -----------------------------------------------------------*/


  float particle_mass = 2.8e4;
  int i, indx1d;

  /* set the masses for the individual elements */
  /* at the moment, this is done relative to the metal mass fractions in the tables */
  
  for (i = 0; i < colibre_cooling_N_elementtypes; i ++) {
    indx1d = row_major_index_2d(cooling->indxZsol, i, colibre_cooling_N_metallicity, colibre_cooling_N_abundances);
    part->metal_mass[i]  =  particle_mass * cooling->MassFractions[indx1d];
    if (i > 1) part->metal_mass[i] = particle_mass * cooling->MassFractions[indx1d] * 0.01;
  }

  float abundance_ratio[colibre_cooling_N_elementtypes];
  float logZZsol;
  logZZsol = abundance_ratio_to_solar(part, cooling, abundance_ratio);

  printf("\n----------------------------------------------------------------------------\n");
  printf("Abundance ratios (n_x / n_H) / (n_x / n_H)_table for logZZsol = %.2f \n", logZZsol);
  for (i = 0; i < colibre_cooling_N_elementtypes; i++){
        printf("%i\t%.8f\n", i, abundance_ratio[i]);
  }

  printf("----------------------------------------------------------------------------\n");

  double log_u_cgs;
  float  d_red, d_met, d_n_H;
  int    red_index, met_index, n_H_index;
  double cooling_rate;


  printf("log U\tlog T\tCooling rate\tHeating rate\tCompton cooling\tn_e/n_H \n");

  /* get indices for redshift, metallicity, and density */
  get_index_1d(cooling->Redshifts, colibre_cooling_N_redshifts, redshift, &red_index, &d_red);
  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol, &met_index, &d_met);
  get_index_1d(cooling->nH, colibre_cooling_N_density, log10(n_H_cgs), &n_H_index, &d_n_H);

  log_u_cgs = log_u_cgs_min;
  do {
    /* get total net cooling rate and print some individual components */
    cooling_rate = colibre_cooling_rate(log_u_cgs, redshift, n_H_cgs, pow(10., logZZsol), abundance_ratio,
                   n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);
    log_u_cgs += dlog_u_cgs;

  } while (log_u_cgs < log_u_cgs_max);


}
