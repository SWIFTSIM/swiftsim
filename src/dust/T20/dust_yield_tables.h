#ifndef SWIFT_DUST_T20_TABLES_H
#define SWIFT_DUST_T20_TABLES_H

/* Config parameters. */
#include "config.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/* Avoid cyclic inclusions */
struct dustevo_props;
struct feedback_props;

/*! Number of elements considered for the SNII yields */
#define eagle_feedback_SNII_N_elements 11

/*! Number of mass bins considered for the SNII yields */
#define eagle_feedback_SNII_N_masses 11

/*! Number of metallicity bins considered for the SNII yields */
#define eagle_feedback_SNII_N_metals 5

/*! Number of elements considered for the AGB yields */
#define eagle_feedback_AGB_N_elements 11

/*! Number of mass bins considered for the AGB yields */
#define eagle_feedback_AGB_N_masses 23

/*! Number of metallicity bins considered for the AGB yields */
#define eagle_feedback_AGB_N_metals 3

/*! CC: Number of elements to be read from the yield tables */
#define enrichment_of_N_elements_from_yield_tables 9

/*! Number of bins used to define the IMF */
#define eagle_feedback_N_imf_bins 200

void read_colibre_depletion(hid_t id, float **log_depletion_fractions,
				  const int table_cooling_N_redshifts,
				  const int table_cooling_N_temperature,
				  const int table_cooling_N_metallicity,
				  const int table_cooling_N_density,
				  const int table_cooling_N_elementtypes);

void depletion_correct_rates(float *cooling_array_heating_rate,
				   float *cooling_array_cooling_rate,
				   float *log_depletion_fractions,
				   const int table_cooling_N_redshifts,
				   const int table_cooling_N_temperature,
				   const int table_cooling_N_metallicity,
				   const int table_cooling_N_density,
				   const int table_cooling_N_elementtypes,
				   const int table_cooling_N_cooltypes,
				   const int table_cooling_N_heattypes);

void initialise_dyield_tables(struct feedback_props *fp, 
			      struct dustevo_props *dp);

void print_dyield_tables(struct feedback_props *fp,
			 struct dustevo_props *dp);

void read_AGB_dyield_tables(struct dustevo_props *dp);

void resample_AGB_dyield(struct feedback_props *fp,
			 struct dustevo_props *dp);


void compute_SNII_dyield(struct feedback_props *fp, 
			 struct dustevo_props *dp);

void compute_AGB_dyield(struct feedback_props *fp, 
			struct dustevo_props *dp);


#endif /* SWIFT_DUST_T20_TABLES_H */
