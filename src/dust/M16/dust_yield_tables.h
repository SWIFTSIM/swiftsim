#ifndef SWIFT_DUST_M16_TABLES_H
#define SWIFT_DUST_M16_TABLES_H

#include <hdf5.h>
#include <string.h>
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

static INLINE void initialise_dyield_tables(struct feedback_props *fp, 
					    struct dustevo_props *dp) {

  /* allocate memory */

  if (posix_memalign((void **)&dp->dyield_SNII.yield_IMF_resampled, 
		     SWIFT_STRUCT_ALIGNMENT,
		     eagle_feedback_SNII_N_metals * grain_species_count * 
		     eagle_feedback_N_imf_bins* sizeof(double)) != 0)
    error("Failed to allocate SNII dust yield array\n");

  if (posix_memalign((void **)&dp->dyield_AGB.yield_IMF_resampled, 
		     SWIFT_STRUCT_ALIGNMENT,
		     eagle_feedback_AGB_N_metals * grain_species_count * 
		     eagle_feedback_N_imf_bins* sizeof(double)) != 0)
    error("Failed to allocate AGB dust yield array\n");

  /* zero arrays */
  memset(dp->dyield_SNII.yield_IMF_resampled, 0.,
  	 eagle_feedback_SNII_N_metals * grain_species_count * eagle_feedback_N_imf_bins *sizeof(double));
  memset(dp->dyield_AGB.yield_IMF_resampled, 0.,
  	 eagle_feedback_AGB_N_metals * grain_species_count * eagle_feedback_N_imf_bins *sizeof(double));
}

static INLINE void print_dyield_tables(struct feedback_props *fp,
				       struct dustevo_props *dp) {
  
  message("downsampled printing of SNII grain yield table:");

  int dyield_index_3d;

  for (int grain = 0; grain < grain_species_count; grain++) {
    for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
      for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
  	dyield_index_3d =
  	  row_major_index_3d(i,
  			     grain,
  			     k,
  			     eagle_feedback_SNII_N_metals,
  			     grain_species_count,
  			     eagle_feedback_N_imf_bins);

  	  if (dyield_index_3d % 10 == -1) {
  	    message("\t\t Grain %d | Metallicity %f | Mass %f | Index %d | Yield %f",
  		    grain,
  		    exp10f(fp->yield_SNII.metallicity[i]),
  		    exp10f(fp->yield_mass_bins[k]),
  		    dyield_index_3d,
  		    dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d]);
  	  }
      }
    }
  }

  message("downsampled printing of AGB grain yield table:");

  for (int grain = 0; grain < grain_species_count; grain++) {
    for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
      for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
  	dyield_index_3d =
  	  row_major_index_3d(i,
  			     grain,
  			     k,
  			     eagle_feedback_AGB_N_metals,
  			     grain_species_count,
  			     eagle_feedback_N_imf_bins);

  	  if (dyield_index_3d % 20 == 0) {
  	    message("\t\t Grain %d | Metallicity %f | Mass %f | Index %d | Yield %f",
  		    grain,
  		    exp10f(fp->yield_AGB.metallicity[i]),
  		    exp10f(fp->yield_mass_bins[k]),
  		    dyield_index_3d,
  		    dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d]);
  	  }
      }
    }
  }
}

static INLINE void compute_SNII_dyield(struct feedback_props *fp, 
				       struct dustevo_props *dp) {

  message("Budgeting SNII dust yield from total metal yields using a Dwek et al (1998) prescription");

  /* variables to store temporary values in calculation */
  double rel_yield, base_yield, dust_yield, dust_yield_O;
  float solar_metalfrac, c_frac, mu_elem_inv; 
  int yield_index_3d, yield_index_2d, yield_index_3d_O;
  int dyield_index_3d, dyield_index_3d_O;
  int eldx;
  bool isRefractory;
  const char* elname;
  
  /* iterate through each grain composition*/
  for (int grain = 0; grain < grain_species_count; grain++) {

    /* grain-specific constants here */
    c_frac = dp->condensation_frac[grain];

    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      
      /* get constituent element chemistry index */
      eldx = dp->grain_element_indices[grain][elem];

      /* is this a refractory element? */
      isRefractory = (eldx != chemistry_element_C && eldx != chemistry_element_O);
      
      /* get element-specific constants here */
      solar_metalfrac = dp->abundance_pattern[eldx] / dp->solar_metallicity;
      mu_elem_inv = 1./dp->atomic_weight[eldx];

      /* print constituent element */
      elname = chemistry_get_element_name((enum chemistry_element)eldx);
      message("\t Synthesise dust from element: %s", elname);

      for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
	for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {

	  /* indices for the resampled gas yield (3d), ejecta (2d) and dust yield (3d) arrays */
	  yield_index_3d =
	    row_major_index_3d(i,eldx,k,
			       eagle_feedback_SNII_N_metals,
			       enrichment_of_N_elements_from_yield_tables,
			       eagle_feedback_N_imf_bins);

	  yield_index_2d =
	    row_major_index_2d(i,k,
			       eagle_feedback_SNII_N_metals,
			       eagle_feedback_N_imf_bins);
	  dyield_index_3d =
	    row_major_index_3d(i,grain,k,
	  		       eagle_feedback_SNII_N_metals,
	  		       grain_species_count,
	  		       eagle_feedback_N_imf_bins);	  

	   
	  /* throughput metal mass in the ejecta */
	  base_yield = solar_metalfrac * \
	    exp10f(fp->yield_SNII.metallicity[i]) * \
	    fp->yield_SNII.ejecta_IMF_resampled[yield_index_2d];
	   
	  /* the additional ejected metal mass produced or destroyed by star */
	  rel_yield = fp->yield_SNII.yield_IMF_resampled[yield_index_3d];

	  /* ejected metal mass condensing in the dust phase */
	  dust_yield = (rel_yield+base_yield) * c_frac;

	  /* subtract metal mass in grains from gas yields */
	  fp->yield_SNII.yield_IMF_resampled[yield_index_3d] = 
	    fp->yield_SNII.yield_IMF_resampled[yield_index_3d] - dust_yield;

	  /* add metal mass in grains to dust yields */
	  dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d] = 
	    dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d] + dust_yield;

	  if (isRefractory) {
	    /** oxygen depletion tied to depletion of refractory elements **/

	    /* get 3d indices for equivalent oxygen gas and dust yields */
	    yield_index_3d_O =
	      row_major_index_3d(i,
				 chemistry_element_O,
				 k,
				 eagle_feedback_SNII_N_metals,
				 enrichment_of_N_elements_from_yield_tables,
				 eagle_feedback_N_imf_bins);
	    dyield_index_3d_O =
	      row_major_index_3d(i,
	    			 grain_species_O,
	    			 k,
	    			 eagle_feedback_SNII_N_metals,
	    			 grain_species_count,
	    			 eagle_feedback_N_imf_bins);
	    
	    /* compute depleted oxygen mass for refractory element */
	    dust_yield_O = 10 * dust_yield * mu_elem_inv;

	    /* subtract depleted oxygen from gas-phase yields */
	    fp->yield_SNII.yield_IMF_resampled[yield_index_3d_O] = 
	      fp->yield_SNII.yield_IMF_resampled[yield_index_3d_O] - dust_yield_O;

	    /* add depleted oxygen to dust-phase yields */
	    dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d_O] = 
	      dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d_O] + dust_yield_O;

	  }

	  if (yield_index_3d % 50 == -1) {
	    message("\t Index %d :: Dindex %d :: :: Mass %f :: Metallicity %f :: Orig Yield %f :: Dust Yield %f :: Mod Yield %f",
		    //elname,
		    yield_index_3d, dyield_index_3d, 
		    exp10f(fp->yield_mass_bins[k]), 
		    exp10f(fp->yield_SNII.metallicity[i]), 
		    rel_yield, 
		    dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d],
		    fp->yield_SNII.yield_IMF_resampled[yield_index_3d]);
		     

	  }
	}  
      }
    }
  }
}

static INLINE void compute_AGB_dyield(struct feedback_props *fp, 
				       struct dustevo_props *dp) {

  message("Budgeting AGB dust yield from total metal yields using a Dwek et al (1998) prescription");

  /* variables to store temporary values in calculation */
  double rel_yield, base_yield, dust_yield, base_yield_O, rel_yield_O, dust_yield_O, base_yield_C, rel_yield_C;
  double dyield_O, dyield_C;
  float solar_metalfrac, c_frac, mu_elem_inv; 
  int yield_index_3d, yield_index_2d, yield_index_3d_O, yield_index_3d_C;
  int dyield_index_3d, dyield_index_3d_O, dyield_index_3d_C;
  int eldx;
  bool isRefractory, isCstar;
  const char* elname;
  
  /* iterate through each grain composition*/
  for (int grain = 0; grain < grain_species_count; grain++) {

    /* grain-specific constants here */
    c_frac = dp->condensation_frac[grain];

    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      
      /* get constituent element chemistry index */
      eldx = dp->grain_element_indices[grain][elem];

      /* is this a refractory element? */
      isRefractory = (eldx != chemistry_element_C && eldx != chemistry_element_O);
      
      /* get element-specific constants here */
      solar_metalfrac = dp->abundance_pattern[eldx] / dp->solar_metallicity;
      mu_elem_inv = 1./dp->atomic_weight[eldx];

      /* print constituent element */
      elname = chemistry_get_element_name((enum chemistry_element)eldx);
      message("\t Synthesise dust from element: %s", elname);

      if (eldx == chemistry_element_O) {
      	/* no direct yields for oxygen */
      	continue;
      }

      for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
	for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {
	  
	  /* indices for the resampled gas yield (3d), ejecta (2d) and dust yield (3d) arrays */
	  yield_index_3d =
	    row_major_index_3d(i, eldx, k,
			       eagle_feedback_AGB_N_metals,
			       enrichment_of_N_elements_from_yield_tables,
			       eagle_feedback_N_imf_bins);

	  yield_index_2d =
	    row_major_index_2d(i,k,
			       eagle_feedback_AGB_N_metals,
			       eagle_feedback_N_imf_bins);
	  dyield_index_3d =
	    row_major_index_3d(i, grain, k,
	  		       eagle_feedback_AGB_N_metals,
	  		       grain_species_count,
	  		       eagle_feedback_N_imf_bins);	  

	  /* specific yield indices for  C and O */
	  yield_index_3d_O =
	    row_major_index_3d(i, chemistry_element_O, k,
			       eagle_feedback_AGB_N_metals,
			       enrichment_of_N_elements_from_yield_tables,
			       eagle_feedback_N_imf_bins);

	  yield_index_3d_C =
	    row_major_index_3d(i, chemistry_element_C, k,
			       eagle_feedback_AGB_N_metals,
			       enrichment_of_N_elements_from_yield_tables,
			       eagle_feedback_N_imf_bins);

	  dyield_index_3d_O =
	    row_major_index_3d(i, grain_species_O, k,
	  		       eagle_feedback_AGB_N_metals,
	  		       grain_species_count,
	  		       eagle_feedback_N_imf_bins);	  

	  dyield_index_3d_C =
	    row_major_index_3d(i, grain_species_C, k,
	  		       eagle_feedback_AGB_N_metals,
	  		       grain_species_count,
	  		       eagle_feedback_N_imf_bins);	  
	   
	  /* throughput metal mass in the ejecta */
	  base_yield = solar_metalfrac * \
	    exp10f(fp->yield_AGB.metallicity[i]) * \
	    fp->yield_AGB.ejecta_IMF_resampled[yield_index_2d];

	  base_yield_O = (dp->abundance_pattern[chemistry_element_O] / \
			  dp->solar_metallicity) * \
	    exp10f(fp->yield_AGB.metallicity[i]) * \
	    fp->yield_AGB.ejecta_IMF_resampled[yield_index_2d];

	  base_yield_C = (dp->abundance_pattern[chemistry_element_C] /
			  dp->solar_metallicity) * \
	    exp10f(fp->yield_AGB.metallicity[i]) * \
	    fp->yield_AGB.ejecta_IMF_resampled[yield_index_2d];
	   
	  /* the additional ejected metal mass produced or destroyed by star */
	  rel_yield = fp->yield_AGB.yield_IMF_resampled[yield_index_3d];
	  rel_yield_O = fp->yield_AGB.yield_IMF_resampled[yield_index_3d_O];
	  rel_yield_C = fp->yield_AGB.yield_IMF_resampled[yield_index_3d_C];
	  
	  /* what mass of C or O is currently assigned to dust grains? */
	  dyield_O = dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d_O];
	  dyield_C = dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d_C];

	  /* is C/O ratio > 1 ? number density ratio of total C to O atoms in ejecta */
	  isCstar = ((rel_yield_C + base_yield_C + dyield_C) * 
		     dp->atomic_weight[chemistry_element_O]) > 
	    ((rel_yield_O + base_yield_O + dyield_O) *
	     dp->atomic_weight[chemistry_element_C]);

	  if (yield_index_3d % 1 == -1) {
	    message("%d %d %d %d M* %f | Z* %f \t\t C/O %e, C,O base %e %e, C,O rel. %e %e, abund %f Z %f yield %f" , 
		    i, 
		    k,
		    yield_index_3d_O, yield_index_3d_C,
		    exp10f(fp->yield_mass_bins[k]), 
		    exp10f(fp->yield_AGB.metallicity[i]), 
		    ((rel_yield_C + base_yield_C + dyield_C) * dp->atomic_weight[chemistry_element_O]) /
		    ((rel_yield_O + base_yield_O + dyield_O)  * dp->atomic_weight[chemistry_element_C]),
		    base_yield_C,base_yield_O, rel_yield_C,rel_yield_O, 
		    /* dp->atomic_weight[chemistry_element_C],dp->atomic_weight[chemistry_element_O] */
		    dp->abundance_pattern[chemistry_element_O],
		    exp10f(fp->yield_AGB.metallicity[chemistry_element_O]),
		    fp->yield_AGB.ejecta_IMF_resampled[yield_index_2d]);
	  }
	  /* Assume only carbon stars produce carbon dust yields */
	  if (isCstar) {
	    if (eldx == chemistry_element_C) {
	      /* dust yields depend on O abundance (via e.g. C locked in CO molecules), no c_frac here (assume 1) */
	      dust_yield = ((rel_yield+base_yield) - 0.75*(rel_yield_O+base_yield_O));

	      /* subtract metal mass in grains from gas yields */
	      fp->yield_AGB.yield_IMF_resampled[yield_index_3d] = 
		fp->yield_AGB.yield_IMF_resampled[yield_index_3d] - dust_yield;
	      
	      /* add metal mass in grains to dust yields */
	      dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d] = 
		dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d] + dust_yield;

	    }
	    
	  }
	  /* Assume only non-carbon stars produce silicate dust yields */
	  else {
	    if (isRefractory) {
	      /* compute refractory dust yields */
	      dust_yield = c_frac*(rel_yield+base_yield);

	      /* subtract metal mass in grains from gas yields */
	      fp->yield_AGB.yield_IMF_resampled[yield_index_3d] = 
		fp->yield_AGB.yield_IMF_resampled[yield_index_3d] - dust_yield;
   
	      /* add metal mass in grains to dust yields */
	      dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d] = 
		dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d] + dust_yield;

	      /** oxygen depletion tied to depletion of refractory elements **/
	      /* compute depleted oxygen mass for refractory element */
	      dust_yield_O = 10 * dust_yield * mu_elem_inv;

	      /* subtract depleted oxygen from gas-phase yields */
	      fp->yield_AGB.yield_IMF_resampled[yield_index_3d_O] =
		fp->yield_AGB.yield_IMF_resampled[yield_index_3d_O] - dust_yield_O;

	      /* add depleted oxygen to dust-phase yields */
	      dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d_O] =
		dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d_O] + dust_yield_O;
	    }
	  }

	  if (yield_index_3d % 10 == -1) {
	    message("\t Ind %d :: Dind %d :: DOind %d :: Oind %d :: Oind2 %d  :: Mass %f :: Metallicity %f :: Orig Yield %f :: Dust Yield %f :: Mod Yield %f",
		    yield_index_3d, dyield_index_3d,  
		    dyield_index_3d_O, yield_index_3d_O,
		    row_major_index_3d(i, chemistry_element_O, k,
				       eagle_feedback_AGB_N_metals,
				       enrichment_of_N_elements_from_yield_tables,
				       eagle_feedback_N_imf_bins),
		    exp10f(fp->yield_mass_bins[k]), 
		    exp10f(fp->yield_AGB.metallicity[i]), 
		    rel_yield, 
		    dp->dyield_AGB.yield_IMF_resampled[dyield_index_3d],
		    fp->yield_AGB.yield_IMF_resampled[yield_index_3d]);
	  }  
	}
      }
    }
  }
}

/* void compute_SNII_yield(struct feedback_props *fp,  */
/* 			struct dustevo_props *dp) { */
/* } */

#endif /* SWIFT_DUST_M16_TABLES_H */
