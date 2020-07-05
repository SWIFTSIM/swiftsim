#ifndef SWIFT_DUST_T20_TABLES_H
#define SWIFT_DUST_T20_TABLES_H

/* dimensions for the COLIBRE cooling tables */
#define colibre_cooling_N_temperature 86
#define colibre_cooling_N_redshifts 46
#define colibre_cooling_N_density 7
#define colibre_cooling_N_metallicity 11
#define colibre_cooling_N_cooltypes 22
#define colibre_cooling_N_heattypes 24
#define colibre_cooling_N_electrontypes 14
#define colibre_cooling_N_elementtypes 12

#include <hdf5.h>

#include "dust.h"

static INLINE void read_colibre_depletion(hid_t id, struct dustevo_props *dp) {
  
  /* HDF5 variables */

  hid_t dataset;
  herr_t status;

  if (posix_memalign(
          (void **)&dp->logfD, SWIFT_STRUCT_ALIGNMENT,
          colibre_cooling_N_redshifts * colibre_cooling_N_temperature *
	  colibre_cooling_N_metallicity * colibre_cooling_N_density *
	  (colibre_cooling_N_elementtypes - 1) * sizeof(float)) != 0)
    error("Failed to allocate depletion array\n");

  dataset = H5Dopen(id, "/Tdep/Depletion", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   dp->logfD);
  if (status < 0) error("error reading dust depletion (temperature)\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

}

#endif /* SWIFT_DUST_T20_TABLES_H */
