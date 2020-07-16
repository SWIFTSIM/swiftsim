#ifndef SWIFT_DUST_M16_TABLES_H
#define SWIFT_DUST_M16_TABLES_H

#include <hdf5.h>
#include "inline.h"
#include "dust.h"
#include "feedback_properties.h"
#include "chemistry.h"
#include "interpolate.h"

/*! Number of elements considered for the SNII yields */
#define eagle_feedback_SNII_N_elements 11

/*! Number of mass bins considered for the SNII yields */
#define eagle_feedback_SNII_N_masses 11

/*! Number of metallicity bins considered for the SNII yields */
#define eagle_feedback_SNII_N_metals 5

/*! CC: Number of elements to be read from the yield tables */
#define enrichment_of_N_elements_from_yield_tables 9

/*! Number of bins used to define the IMF */
#define eagle_feedback_N_imf_bins 200

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

static INLINE void compute_SNII_yield(struct feedback_props *fp, 
				      struct dustevo_props *dp) {

  /*** NOTE: total metallicity yields untouched here, so Z represents the conserved dust + gas phase metals ***/
  message("Budgeting SNII dust yield from total metal yields using a Dwek et al (1998) prescription");

  /* variables to store temporary values in calculation */
  float rel_yield, base_yield, dust_yield, dust_yield_O;
  float solar_metalfrac, c_frac, mu_elem_inv; 
  int flat_index_3d, flat_index_2d, flat_index_3d_O;
  int eldx;
  bool isRefractory;
  const char* elname;

  
  /* iterate through each grain composition*/
  for (int grain = 0; grain < grain_species_count; grain++) {
    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      
      /* get constituent element chemistry index */
      eldx = dp->grain_element_indices[grain][elem];

      /* is this a refractory element? */
      isRefractory = (eldx != chemistry_element_C && eldx != chemistry_element_O);
      
      /* get element-specific constants here */
      solar_metalfrac = dp->abundance_pattern[eldx];
      mu_elem_inv = 1./dp->atomic_weight[eldx];
      c_frac = dp->condensation_frac[eldx];

      /* print constituent element */
      elname = chemistry_get_element_name((enum chemistry_element)eldx);
      message("\t Synthesise dust from element: %s", elname);

       for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
  	 for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {

	   /* indices for the resampled yield (3d) and ejecta (2d) arrays */
  	   flat_index_3d =
  	     row_major_index_3d(i,
  				eldx,
  				k,
  				eagle_feedback_SNII_N_metals,
  				enrichment_of_N_elements_from_yield_tables,
  				eagle_feedback_N_imf_bins);

  	   flat_index_2d =
  	     row_major_index_2d(i,
  				k,
  				eagle_feedback_SNII_N_metals,
  				eagle_feedback_N_imf_bins);
	   
	   /* throughput metal mass in the ejecta */
  	   base_yield = solar_metalfrac * \
	     		exp10f(fp->yield_SNII.metallicity[i]) * \
	     		fp->yield_SNII.ejecta_IMF_resampled[flat_index_2d];
	   
	   /* the additional ejected metal mass produced or destroyed by star */
	   rel_yield = fp->yield_SNII.yield_IMF_resampled[flat_index_3d];

	   /* ejected metal mass condensing in the dust phase */
	   dust_yield = (rel_yield+base_yield) * c_frac;

	   /* subtract metal mass in grains from gas yields */
	   fp->yield_SNII.yield_IMF_resampled[flat_index_3d] = 
	     fp->yield_SNII.yield_IMF_resampled[flat_index_3d] - dust_yield;

	   // <TD> add dust_yield to own array

	   if (isRefractory) {
	     /* oxygen depletion tied to depletion of refractory elements */

	     /* get 3d index for equivalent oxygen yield */
	     flat_index_3d_O =
	       row_major_index_3d(i,
				  chemistry_element_O,
				  k,
				  eagle_feedback_SNII_N_metals,
				  enrichment_of_N_elements_from_yield_tables,
				  eagle_feedback_N_imf_bins);

	     /* compute depleted oxygen mass for refractory element */
	     dust_yield_O = 10 * dust_yield * mu_elem_inv;

	     /* subtract depleted oxygen from gas-phase yields */
	     fp->yield_SNII.yield_IMF_resampled[flat_index_3d_O] = 
	       fp->yield_SNII.yield_IMF_resampled[flat_index_3d_O] - dust_yield_O;

	     // <TD> add oxygen in oxide dust_yield to own array	     
	     
	   }

	   if (0) {
	     message("   Element %s (%f) :: Index %d :: Mass %f :: Metallicity %f :: Rel Yield %f :: Base Yield %f :: Tot Yield %f",
		     elname, solar_metalfrac * 0.0127,
		     flat_index_3d, 
		     exp10f(fp->yield_mass_bins[k]), 
		     exp10f(fp->yield_SNII.metallicity[i]), 
		     rel_yield, 
		     base_yield,
		     dust_yield);
		     

	   }
	 }  
       }
     }
  }
  exit(0);
}
/* void compute_SNII_yield(struct feedback_props *fp,  */
/* 			struct dustevo_props *dp) { */
/* } */

#endif /* SWIFT_DUST_M16_TABLES_H */
