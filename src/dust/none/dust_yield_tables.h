#ifndef SWIFT_DUST_NONE_TABLES_H
#define SWIFT_DUST_NONE_TABLES_H

#include <hdf5.h>
#include "inline.h"
#include "dust.h"

static INLINE void depletion_correct_rates(float *cooling_array_heating_rate,
					   float *cooling_array_cooling_rate,
					   float *log_depletion_fraction,
					   const int table_cooling_N_redshifts,
					   const int table_cooling_N_temperature,
					   const int table_cooling_N_metallicity,
					   const int table_cooling_N_density,
					   const int table_cooling_N_elementtypes,
					   const int table_cooling_N_cooltypes,
					   const int table_cooling_N_heattypes) {
}

static INLINE void read_colibre_depletion(hid_t id,
					  float **log_depletion_fractions,
					  const int table_cooling_N_redshifts,
					  const int table_cooling_N_temperature,
					  const int table_cooling_N_metallicity,
					  const int table_cooling_N_density,
					  const int table_cooling_N_elementtypes) {
}

#endif /* SWIFT_DUST_NONE_TABLES_H */
