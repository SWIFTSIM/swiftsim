#ifndef SWIFT_DUST_M16_TABLES_H
#define SWIFT_DUST_M16_TABLES_H

/* dimensions for the COLIBRE cooling tables */
/* #define table_cooling_N_temperature 86 */
/* #define table_cooling_N_redshifts 46 */
/* #define table_cooling_N_density 7 */
/* #define table_cooling_N_metallicity 11 */
/* #define table_cooling_N_cooltypes 22 */
/* #define table_cooling_N_heattypes 24 */
/* #define table_cooling_N_electrontypes 14 */
/* #define table_cooling_N_elementtypes 12 */

#include <hdf5.h>
#include "inline.h"
#include "dust.h"

static INLINE void read_colibre_depletion(hid_t id, float **log_depletion_fractions,
					  const int table_cooling_N_redshifts,
					  const int table_cooling_N_temperature,
					  const int table_cooling_N_metallicity,
					  const int table_cooling_N_density,
					  const int table_cooling_N_elementtypes) {

  hid_t dataset;
  herr_t status;

  if (posix_memalign(
	  (void **)log_depletion_fractions, SWIFT_STRUCT_ALIGNMENT,
          table_cooling_N_redshifts * table_cooling_N_temperature *
	  table_cooling_N_metallicity * table_cooling_N_density *
	  (table_cooling_N_elementtypes - 1) * sizeof(float)) != 0)
    error("Failed to allocate depletion array\n");

  dataset = H5Dopen(id, "/Tdep/Depletion", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   *log_depletion_fractions);
  if (status < 0) error("error reading dust depletion (temperature)\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");
}

static INLINE void depletion_correct_rates(float *cooling_array_heating_rate,
					   float *cooling_array_cooling_rate,
					   float *log_depletion_fractions,
					   const int table_cooling_N_redshifts,
					   const int table_cooling_N_temperature,
					   const int table_cooling_N_metallicity,
					   const int table_cooling_N_density,
					   const int table_cooling_N_elementtypes,
					   const int table_cooling_N_cooltypes,
					   const int table_cooling_N_heattypes){
  /* index for depletion array */
  int idx = 0;

  /* indices for rate arrays */
  int idx_cool = 0;
  int idx_heat = 0;

  /* size difference of last array dimension relative to depletion array */
  int cool_len_diff = table_cooling_N_cooltypes - (table_cooling_N_elementtypes-1);
  int heat_len_diff = table_cooling_N_heattypes - (table_cooling_N_elementtypes-1);

  // iterate through cooling and heating table dimensions 
  for (int i = 0; i < table_cooling_N_redshifts; i++) {
    for (int j = 0; j < table_cooling_N_temperature; j++) {
      for (int k = 0; k < table_cooling_N_metallicity; k++) {
	for (int l = 0; l < table_cooling_N_density; l++) {
	  for (int m = 0; m < (table_cooling_N_elementtypes-1); m++) {
	    // For each cooling and heating rate we divide out by the
	    // fraction of element m in the gas phase to remove
	    // implicit depletion 
	    const float logfgas = log10f(1.-exp10f(log_depletion_fractions[idx]));

	    cooling_array_heating_rate[idx_heat] -= logfgas;
	    cooling_array_cooling_rate[idx_cool] -= logfgas;
	    
	    /* increment indices */
	    idx = idx + 1;
	    idx_cool = idx_cool + 1; 
	    idx_heat = idx_heat + 1;
	  }
	  // increment cooling and heating arrays to skip non named element channels
	  idx_cool = idx_cool + cool_len_diff;
	  idx_heat = idx_heat + heat_len_diff;
	}
      }
    }
  } 
  message("Scaled implicit dust depletion out of Colibre cooling tables.");
}

#endif /* SWIFT_DUST_M16_TABLES_H */
