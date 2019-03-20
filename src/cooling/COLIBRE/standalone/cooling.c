#include <stdio.h>
#include <string.h>

#include "physical_constants_cgs.h"
#include "cooling_tables.h"


/**
 * @brief Initialises properties stored in the cooling_function_data struct
 * and in the particle_data struct
 * @param cooling #cooling_function_data struct to initialize
 * @param part #particle_data struct to initialize
 */
void cooling_init_backend(struct cooling_function_data *cooling, struct particle_data *part) {


  /* Setting the path for the hdf5 file */
  strcpy(cooling->cooling_table_path,"/cosma7/data/dp004/dc-ploe1/CoolingTables/2019_04/UV_dust1_CR1_G1_shield1.hdf5");
  printf("\n----------------------------------------------------------------------------\n");
  printf("\n Reading in from file: \n");
  printf("     %s\n",cooling->cooling_table_path);

  /* Calculate constant for analytic Compton cooling */ 
  double radiation_constant =  4. * const_stefan_boltzmann_cgs / const_speed_light_c_cgs;
  cooling->Compton_constant_cgs = 4. * const_boltzmann_k_cgs / const_electron_mass_cgs / const_speed_light_c_cgs *
                                  const_thomson_cross_section_cgs * radiation_constant;

  /* Reading in the table bins and abundances */
  read_cooling_header(cooling);

  /* Allocating the memory and reading in the multi-dimensional datasets for the cooling rates et al. */
  read_cooling_tables(cooling);

  /* Allocate 'particle' data */
  allocate_part_data(part);
 
}

