
/* This file's header */
#include "dust_yield_tables.h"

/* Standard headers */
#include <math.h>
#include <string.h>

/* Local headers */
#include "dust.h"
#include "feedback_properties.h"
#include "interpolate.h"
#include "chemistry.h"


void read_colibre_depletion(hid_t id, float **log_depletion_fractions,
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

void depletion_correct_rates(float *cooling_array_heating_rate,
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

void initialise_dyield_tables(struct feedback_props *fp, 
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

void print_dyield_tables(struct feedback_props *fp,
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



/**
 * @brief reads yield tables, flattens and stores them in stars_props data
 * struct
 *
 * @param feedback_props the #feedback_props data struct to read the table into.
 */
void read_AGB_dyield_tables(struct dustevo_props *dp) {

#ifdef HAVE_HDF5

  /* filenames to read HDF5 files */
  char fname[256], setname[100];

  hid_t file_id, dataset, dataset2, datatype, dataspace;
  herr_t status;

  /* Allocate array to store AGB dyield tables */
  if (swift_memalign("feedback-tables",
                     (void **)&dp->dyield_AGB.yield,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_AGB_N_metals * eagle_feedback_AGB_N_masses *
                         grain_species_count * sizeof(double)) != 0) {
    error("Failed to allocate AGB yield array");
  }

  /* Read AGB tables */
  sprintf(fname, "%s/AGB_dustyield.hdf5", dp->AGB_dyield_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
if (file_id < 0) error("unable to open file %s\n", fname);      

  double temp_yield_AGB[grain_species_count]
       [eagle_feedback_AGB_N_masses];
  char *metallicity_yield_table_name_AGB[eagle_feedback_AGB_N_metals];

  /* read metallicity names */
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset2 = H5Dopen(file_id, "Yield_names", H5P_DEFAULT);
  dataspace = H5Dget_space(dataset2);

  status = H5Dread(dataset2, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   metallicity_yield_table_name_AGB);
  if (status < 0) error("error reading yield table names");


  /* read AGB yield tables */
  for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
    /* read yields to temporary array */
    sprintf(setname, "/Yields/%s/Yield", metallicity_yield_table_name_AGB[i]);
    dataset = H5Dopen(file_id, setname, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temp_yield_AGB);
    if (status < 0) error("error reading AGB yield");
    status = H5Dclose(dataset);
    if (status < 0) error("error closing dataset");
    /* Flatten the temporary tables that were read, store in stars_props */
    for (int k = 0; k < eagle_feedback_AGB_N_masses; k++) {

      for (int j = 0; j < grain_species_count; j++) {
        const int flat_index_Z = row_major_index_3d(
            i, j, k, eagle_feedback_AGB_N_metals, grain_species_count,
            eagle_feedback_AGB_N_masses);
        dp->dyield_AGB.yield[flat_index_Z] = temp_yield_AGB[j][k];
      }
    }
  }
  
  /* Release the memory allocated by HDF5 for the strings */
  status = H5Dvlen_reclaim(datatype, dataspace, H5P_DEFAULT,
                           metallicity_yield_table_name_AGB);
  if (status < 0) error("error freeing string memory");
  status = H5Dclose(dataset2);
  if (status < 0) error("error closing dataset");
  status = H5Tclose(datatype);
  if (status < 0) error("error closing datatype");
  status = H5Sclose(dataspace);
  if (status < 0) error("error closing dataspace");

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");
#endif
}
/**
 * @brief resamples yields based on IMF mass bins
 *
 * @param feedback_props the #feedback_props data struct.
 */
INLINE  void resample_AGB_dyield(struct feedback_props *fp,
				       struct dustevo_props *dp) { //<HERE>

  int flat_index_3d  ;


  /* Declare temporary tables to accumulate yields */
  double AGB_yield[eagle_feedback_AGB_N_masses];
  float result;

  /* Resample yields for each element tracked in COLIBRE */
  for (enum grain_species grain = 0;
       grain < grain_species_count; grain++) {

    for (int i = 0; i < eagle_feedback_AGB_N_metals; i++) {
      for (int j = 0; j < eagle_feedback_AGB_N_masses; j++) {
	flat_index_3d = row_major_index_3d(
					   i, grain, j, eagle_feedback_AGB_N_metals,
					   eagle_feedback_AGB_N_elements, eagle_feedback_AGB_N_masses);
	AGB_yield[j] = fp->yield_AGB.yield[flat_index_3d] *
	  exp(M_LN10 * (-fp->yield_AGB.mass[j]));
      }

      for (int j = 0; j < eagle_feedback_N_imf_bins; j++) {
	if (fp->yield_mass_bins[j] <
	    fp->yield_AGB.mass[0])
	  result = AGB_yield[0];
	else if (fp->yield_mass_bins[j] >
		 fp->yield_AGB
		 .mass[eagle_feedback_AGB_N_masses - 1])
	  result = AGB_yield[eagle_feedback_AGB_N_masses - 1];
	else
	  result = interpolate_1D_non_uniform(
					      fp->yield_AGB.mass, AGB_yield,
					      eagle_feedback_AGB_N_masses,
					      fp->yield_mass_bins[j]);

	flat_index_3d =
	  row_major_index_3d(i, grain, j, eagle_feedback_AGB_N_metals,
			     enrichment_of_N_elements_from_yield_tables,
			     eagle_feedback_N_imf_bins);
	fp->yield_AGB.yield_IMF_resampled[flat_index_3d] =
	  exp(M_LN10 * fp->yield_mass_bins[j]) * result;
      }
    }
  }  
}


void compute_SNII_dyield(struct feedback_props *fp, 
				       struct dustevo_props *dp) {

  message("Budgeting SNII dust yield from total metal yields using a Zhukovska, Gail & Trieloff (2008) prescription");

  /* variables to store temporary values in calculation */
  double rel_yield, base_yield, dust_yield/*, dust_yield_O*/;
  float solar_metalfrac, c_frac, elfrac;/* , mu_elem_inv */; 
  int yield_index_3d, yield_index_2d/* , yield_index_3d_O */;
  int dyield_index_3d/* , dyield_index_3d_O */;
  int eldx;
  const char* elname;
  
  /* first pass through each grain composition to find bottlenecked yield for each grain type */
  for (int grain = 0; grain < grain_species_count; grain++) {
    
    /* grain-specific constants here */
    c_frac = dp->condensation_frac[grain];

    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {  
      /* get constituent element chemistry index */
      eldx = dp->grain_element_indices[grain][elem];

      /* only Mg and Si are considered key elements for silicate grains */
      if ((grain==1) && (eldx != chemistry_element_Mg) && (eldx != chemistry_element_Si)){
	continue;
      }

      /* fraction of grain molecular mass constituted by this element */
      elfrac = dp->grain_element_mfrac[grain][elem];

      /* get element-specific constants here */
      solar_metalfrac = dp->abundance_pattern[eldx] / dp->solar_metallicity;

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
	  dust_yield = (rel_yield+base_yield) * c_frac / elfrac;
       
	  /* if dust yield is unset or estimated larger than predicted for this element, set */
	    if (dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d] != 0){
	      dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d] = 
		min(dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d], dust_yield);
	    }
	    else {
	      dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d] = dust_yield;
	    }
	}  
      }
    }
    
    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
     
      /* get constituent element chemistry index */
      eldx = dp->grain_element_indices[grain][elem];

      /* fraction of grain molecular mass constituted by this element */
      elfrac = dp->grain_element_mfrac[grain][elem];

      /* print constituent element */
      elname = chemistry_get_element_name((enum chemistry_element)eldx);
      message("\t Synthesise dust from element: %s", elname);

      for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
	for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {

	  /* indices for the resampled gas yield (3d) and dust yield (3d) arrays */
	  yield_index_3d =
	    row_major_index_3d(i,eldx,k,
			       eagle_feedback_SNII_N_metals,
			       enrichment_of_N_elements_from_yield_tables,
			       eagle_feedback_N_imf_bins);
	  dyield_index_3d =
	    row_major_index_3d(i,grain,k,
	  		       eagle_feedback_SNII_N_metals,
	  		       grain_species_count,
	  		       eagle_feedback_N_imf_bins);	  

	   
	    /* subtract metal mass in grains from gas yields */
	    fp->yield_SNII.yield_IMF_resampled[yield_index_3d] =
	      fp->yield_SNII.yield_IMF_resampled[yield_index_3d] - 
	      (dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d] * elfrac);
	}
      }
    }
  }
}

void compute_AGB_dyield(struct feedback_props *fp, 
				       struct dustevo_props *dp) {
  message("Reading AGB dust yield from %s", dp->AGB_dyield_path);
  read_AGB_dyield_tables(dp);
  resample_AGB_dyield(fp, dp); 

  message("Budgeting AGB dust yield from tabulated metal yields");
  /* variables to store temporary values in calculation */
  double dust_contr, elfrac;
  int yield_index_3d;
  int dyield_index_3d;
  int eldx;
  const char* elname;
  
  /* iterate through each grain composition*/
  for (int grain = 0; grain < grain_species_count; grain++) {
    for (int elem = 0; elem < dp->grain_element_count[grain]; elem++) {
      
      /* get constituent element chemistry index */
      eldx = dp->grain_element_indices[grain][elem];

      /* fraction of grain molecular mass constituted by this element */
      elfrac = dp->grain_element_mfrac[grain][elem];

      /* print constituent element */
      elname = chemistry_get_element_name((enum chemistry_element)eldx);
      message("\t Synthesise dust from element: %s", elname);

      for (int i = 0; i < eagle_feedback_SNII_N_metals; i++) {
	for (int k = 0; k < eagle_feedback_N_imf_bins; k++) {

	  /* indices for the resampled gas yield (3d) and dust yield (3d) arrays */
	  yield_index_3d =
	    row_major_index_3d(i,eldx,k,
			       eagle_feedback_SNII_N_metals,
			       enrichment_of_N_elements_from_yield_tables,
			       eagle_feedback_N_imf_bins);

	  dyield_index_3d =
	    row_major_index_3d(i,grain,k,
	  		       eagle_feedback_SNII_N_metals,
	  		       grain_species_count,
	  		       eagle_feedback_N_imf_bins);	  

	  /* element mass in yield contributing to dust yield */
	  dust_contr = elfrac*dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d];
	   
	  /* subtract metal mass in grains from gas yields */
	  fp->yield_SNII.yield_IMF_resampled[yield_index_3d] = 
	    fp->yield_SNII.yield_IMF_resampled[yield_index_3d] - dust_contr;

	  if (yield_index_3d % 50 == -1) {
	    message("\t Index %d :: Dindex %d :: :: Mass %f :: Metallicity %f :: Dust Yield %f :: Mod Yield %f",
		    //elname,
		    yield_index_3d, dyield_index_3d, 
		    exp10f(fp->yield_mass_bins[k]), 
		    exp10f(fp->yield_SNII.metallicity[i]), 
		    dp->dyield_SNII.yield_IMF_resampled[dyield_index_3d],
		    fp->yield_SNII.yield_IMF_resampled[yield_index_3d]);
		     

	  }
	}  
      }
    }
  }
}
