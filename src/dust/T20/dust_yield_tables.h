#ifndef SWIFT_DUST_T20_TABLES_H
#define SWIFT_DUST_T20_TABLES_H

/* dimensions for the COLIBRE cooling tables */
#define table_cooling_N_temperature 86
#define table_cooling_N_redshifts 46
#define table_cooling_N_density 7
#define table_cooling_N_metallicity 11
#define table_cooling_N_cooltypes 22
#define table_cooling_N_heattypes 24
#define table_cooling_N_electrontypes 14
#define table_cooling_N_elementtypes 12

#include <hdf5.h>
#include "inline.h"
#include "dust.h"

static INLINE void depletion_correct_rates(struct cooling_function_data *cooling,
					   struct dustevo_props *dp){

  /* initialise variable to store (log) gas-phase element fraction */
  float logfgas;
  /* initialise index for 5D array */
  int idx = 0;

  /* iterate through cooling and heating table dimensions */
  for (int i = 0; i < table_cooling_N_redshifts; i++) {
    for (int j = 0; j < table_cooling_N_temperature; j++) {
      for (int k = 0; k < table_cooling_N_metallicity; k++) {
	for (int l = 0; l < table_cooling_N_density; l++) {
	  for (int m = 0; m < (table_cooling_N_elementtypes-1); m++) {

	    /* For each cooling and heating rate we divide out by the
	       fraction of element m in the gas phase to remove
	       implicit depletion */
	    logfgas = log10(1-pow(10, dp->logfD[idx]));

	    cooling->table.Theating[idx] =
	      cooling->table.Theating[idx] - logfgas;
	    cooling->table.Tcooling[idx] =
	      cooling->table.Tcooling[idx] - logfgas;

	    idx = idx + 1;
	  }
	}
      }
    }
  }
}


static INLINE void read_colibre_depletion(hid_t id, struct dustevo_props *dp) {
  
  /* HDF5 variables */

  hid_t dataset;
  herr_t status;

  if (posix_memalign(
          (void **)&dp->logfD, SWIFT_STRUCT_ALIGNMENT,
          table_cooling_N_redshifts * table_cooling_N_temperature *
	  table_cooling_N_metallicity * table_cooling_N_density *
	  (table_cooling_N_elementtypes - 1) * sizeof(float)) != 0)
    error("Failed to allocate depletion array\n");

  dataset = H5Dopen(id, "/Tdep/Depletion", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   dp->logfD);
  if (status < 0) error("error reading dust depletion (temperature)\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

}

#endif /* SWIFT_DUST_T20_TABLES_H */
