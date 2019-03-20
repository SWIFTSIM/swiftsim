#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>

#include "cooling_tables.h"
#include "interpolate.h"

/**
 * @brief Reads in COLIBRE cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, metallicity,
 * abundance ratios, and elements used to index the cooling tables.
 *
 * @param cooling Cooling data structure
 */
void read_cooling_header(struct cooling_function_data *cooling) {

  hid_t dataset;
  herr_t status;

  /* read sizes of array dimensions */
  hid_t tempfile_id = H5Fopen(cooling->cooling_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0) printf("unable to open file %s\n", cooling->cooling_table_path);

  /* allocate arrays of bins */
  if (posix_memalign((void **)&cooling->Temp, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_temperature * sizeof(float)) != 0)
    printf("Failed to allocate temperature table\n");

  if (posix_memalign((void **)&cooling->Redshifts, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_redshifts * sizeof(float)) != 0)
    printf("Failed to allocate redshift table\n");

  if (posix_memalign((void **)&cooling->nH, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_density * sizeof(float)) != 0)
    printf("Failed to allocate density table\n");

  if (posix_memalign((void **)&cooling->Metallicity, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * sizeof(float)) != 0)
    printf("Failed to allocate metallicity table\n");

  if (posix_memalign((void **)&cooling->Therm, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_internalenergy * sizeof(float)) != 0)
    printf("Failed to allocate internal energy table\n");

  if (posix_memalign((void **)&cooling->LogAbundances,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * colibre_cooling_N_abundances
                     * sizeof(float)) != 0)
    printf("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&cooling->Abundances,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * colibre_cooling_N_abundances
                     * sizeof(float)) != 0)
    printf("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&cooling->Abundances_inv,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * colibre_cooling_N_abundances
                     * sizeof(float)) != 0)
    printf("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&cooling->atomicmass,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_elementtypes
                     * sizeof(float)) != 0)
    printf("Failed to allocate atomic masses array\n");

  if (posix_memalign((void **)&cooling->Zsol,
                     SWIFT_STRUCT_ALIGNMENT,
                     1 * sizeof(float)) != 0)
    printf("Failed to allocate solar metallicity array\n");

  if (posix_memalign((void **)&cooling->LogMassFractions,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * colibre_cooling_N_abundances
                     * sizeof(float)) != 0)
    printf("Failed to allocate log mass fraction array\n");

  if (posix_memalign((void **)&cooling->MassFractions,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * colibre_cooling_N_abundances
                     * sizeof(float)) != 0)
    printf("Failed to allocate mass fraction array\n");


  /* read in bins and misc information */
  dataset = H5Dopen(tempfile_id, "/TableBins/TemperatureBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->Temp);
  if (status < 0) printf("error reading temperature bins\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/RedshiftBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->Redshifts);
  if (status < 0) printf("error reading redshift bins\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/DensityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->nH);
  if (status < 0) printf("error reading density bins\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/MetallicityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->Metallicity);
  if (status < 0) printf("error reading metallicity bins\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/InternalEnergyBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->Therm);
  if (status < 0) printf("error reading internal energy bins\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TotalAbundances", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->LogAbundances);
  if (status < 0) printf("error reading total abundances\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/ElementMasses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->atomicmass);
  if (status < 0) printf("error reading element masses\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/SolarMetallicity", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->Zsol);
  if (status < 0) printf("error reading solar metallicity \n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TotalMassFractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->LogMassFractions);
  if (status < 0) printf("error reading total mass fractions\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* find the metallicity bin that refers to solar metallicity */
  int i,j, indx1d;
  float tol = 1.e-3;
  for (i = 0; i < colibre_cooling_N_metallicity; i++) {
	  if (fabs(cooling->Metallicity[i]) < tol) {
             cooling->indxZsol = i;
	  }
  }

  /* set some additional useful abundance arrays */
  for (i = 0; i < colibre_cooling_N_metallicity; i++){
      for (j = 0; j < colibre_cooling_N_abundances; j++){
          indx1d = row_major_index_2d(i, j, colibre_cooling_N_metallicity, colibre_cooling_N_abundances);
          cooling->Abundances[indx1d]     = pow(10., cooling->LogAbundances[indx1d]);
          cooling->Abundances_inv[indx1d] = 1. / cooling->Abundances[indx1d];
	  cooling->MassFractions[indx1d]  = pow(10., cooling->LogMassFractions[indx1d]);
      }
  }

}

/**
 *  @brief Allocate space for cooling tables and read them
 * 
 *  @param cooling #cooling_function_data structure
 **/
void read_cooling_tables(struct cooling_function_data *restrict cooling) {
  hid_t dataset;
  herr_t status;


  /* open hdf5 file */
  hid_t tempfile_id = H5Fopen(cooling->cooling_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0) printf("unable to open file %s\n", cooling->cooling_table_path);

  /* Allocate and read arrays to store cooling tables. */
  /* Cooling (temperature) */
  if (posix_memalign((void **)&cooling->table.Tcooling,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_temperature * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * colibre_cooling_N_cooltypes * 
                     sizeof(float)) != 0)
    printf("Failed to allocate Tcooling array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/Cooling", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.Tcooling);
  if (status < 0) printf("error reading Tcooling\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* Cooling (internal energy) */
  if (posix_memalign((void **)&cooling->table.Ucooling,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_internalenergy * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * colibre_cooling_N_cooltypes * 
                     sizeof(float)) != 0)
    printf("Failed to allocate Ucooling array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/Cooling", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.Ucooling);
  if (status < 0) printf("error reading Ucooling\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* Heating (temperature) */
  if (posix_memalign((void **)&cooling->table.Theating,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_temperature * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * colibre_cooling_N_heattypes * 
                     sizeof(float)) != 0)
    printf("Failed to allocate Theating array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/Heating", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.Theating);
  if (status < 0) printf("error reading Theating\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* Heating (internal energy) */
  if (posix_memalign((void **)&cooling->table.Uheating,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_internalenergy * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * colibre_cooling_N_heattypes * 
                     sizeof(float)) != 0)
    printf("Failed to allocate Uheating array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/Heating", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.Uheating);
  if (status < 0) printf("error reading Uheating\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* Electron fraction (temperature) */
  if (posix_memalign((void **)&cooling->table.Telectron_fraction,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_temperature * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * colibre_cooling_N_electrontypes * 
                     sizeof(float)) != 0)
    printf("Failed to allocate Telectron_fraction array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/ElectronFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.Telectron_fraction);
  if (status < 0) printf("error reading electron_fraction (temperature)\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* Electron fraction (internal energy) */
  if (posix_memalign((void **)&cooling->table.Uelectron_fraction,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_internalenergy * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * colibre_cooling_N_electrontypes * 
                     sizeof(float)) != 0)
    printf("Failed to allocate Uelectron_fraction array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/ElectronFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.Uelectron_fraction);
  if (status < 0) printf("error reading electron_fraction (internal energy)\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

  /* Internal energy from temperature */
  if (posix_memalign((void **)&cooling->table.U_from_T,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_temperature * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * 
                     sizeof(float)) != 0)
    printf("Failed to allocate U_from_T array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/U_from_T", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.U_from_T);
  if (status < 0) printf("error reading U_from_T array\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");


  /* Temperature from interal energy */
  if (posix_memalign((void **)&cooling->table.T_from_U,
                     SWIFT_STRUCT_ALIGNMENT, colibre_cooling_N_redshifts * 
                     colibre_cooling_N_internalenergy * colibre_cooling_N_metallicity * 
                     colibre_cooling_N_density * 
                     sizeof(float)) != 0)
    printf("Failed to allocate T_from_U array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/T_from_U", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cooling->table.T_from_U);
  if (status < 0) printf("error reading T_from_U array\n");
  status = H5Dclose(dataset);
  if (status < 0) printf("error closing cooling dataset");

}


/**
 *  @brief converts metal masses to abundance ratios relative to solar abundance ratios
 *  for their metallicity and returns the metallicity in log10 relative to solar
 *  @param part        #particle_data structure 
 *  @param cooling     #cooling_function_data structure 
 *  @param ratio_solar #(n_x / n_H) / (n_x / n_H)_sol  
 */


float abundance_ratio_to_solar(struct particle_data *part, struct cooling_function_data *cooling,
    float ratio_solar[colibre_cooling_N_elementtypes]) {

    int i, indx1d;
    float totmass = 0., metalmass = 0., logZZsol;
    int met_index;
    float d_met;
    float log_nx_nH_sol, log_nx_nH_min, log_nx_nH_max, log_nx_nH;

    /* from masses to abundances (nx/nH) */
    for (i = 0; i < colibre_cooling_N_elementtypes; i++){
	indx1d = row_major_index_2d(cooling->indxZsol, i, colibre_cooling_N_metallicity, colibre_cooling_N_abundances);
   	ratio_solar[i] = part->metal_mass[i] / part->metal_mass[element_H] * 
                         cooling->atomicmass[element_H] / cooling->atomicmass[i] * 
                         cooling->Abundances_inv[indx1d];
	totmass += part->metal_mass[i];
	if (i > element_He) metalmass += part->metal_mass[i];
    }

    /* at this point ratio_solar is (nx/nH) / (nx/nH)_sol */
    /* to multiply with the tables, we want the individual abundance ratio relative */
    /* to what is used in the tables for each metallicity */ 

    /* for example: for a metallicity of 1 per cent solar, the metallicity bin */
    /* for logZZsol = -2 has already the reduced cooling rates for each element; */
    /* it should therefore NOT be multiplied by 0.01 again */
    /* BUT: if e.g. Carbon is twice as abundant as the solar abundance ratio, */
    /* i.e. nC / nH = 0.02 * (nC/nH)_sol for the overall metallicity of 0.01, */
    /* the Carbon cooling rate is multiplied by 2 */

    logZZsol = log10(metalmass / totmass / cooling->Zsol[0]);

    for (i = 0; i < colibre_cooling_N_elementtypes; i++){
        get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol, &met_index, &d_met);
	log_nx_nH_min = cooling->LogAbundances[row_major_index_2d(met_index  , i, colibre_cooling_N_metallicity, colibre_cooling_N_abundances)];
	log_nx_nH_max = cooling->LogAbundances[row_major_index_2d(met_index+1, i, colibre_cooling_N_metallicity, colibre_cooling_N_abundances)];
	log_nx_nH_sol = cooling->LogAbundances[row_major_index_2d(cooling->indxZsol, i, colibre_cooling_N_metallicity, colibre_cooling_N_abundances)];
 	log_nx_nH     = (log_nx_nH_min * (1. - d_met) + log_nx_nH_max * d_met) - log_nx_nH_sol;
	ratio_solar[i] = ratio_solar[i] / pow(10., log_nx_nH);
    }

    /* at this point ratio_solar is (nx/nH) / (nx/nH)_table [Z], */
    /* the metallicity dependent abundance ratio for solar abundances */
    return logZZsol;

}


/**
 *  @brief Allocates particle data
 *  @param part        #particle_data structure 
 */
void allocate_part_data (struct particle_data *part) {

  /* Allocate particle data (metal mass fractions) */
  if (posix_memalign((void **)&part->metal_mass_fraction, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    printf("Failed to allocate particle data (metal_mass_fraction)\n");

  /* Allocate particle data (metal mass) */
  if (posix_memalign((void **)&part->metal_mass, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    printf("Failed to allocate particle data (metal_mass)\n");

}


