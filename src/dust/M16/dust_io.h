#ifndef SWIFT_DUST_IO_M16_H
#define SWIFT_DUST_IO_M16_H

#include "dust.h"
#include "io_properties.h"

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology flag indicating if running with cosmology.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int dust_write_particles(const struct part* parts,
                                            struct io_props* list,
                                            const int with_cosmology) {

  /* List what we want to write */
  list[0] = io_make_output_field(
      "DustMassFractions", FLOAT, grain_species_count,
      UNIT_CONV_NO_UNITS, 0.f, parts, dust_data.grain_mass_fraction,
      "Fractions of the particles' masses that are in a given species of dust grain");

  return 1;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of dust evolution to the file
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The #engine.
 */
INLINE static void dust_write_flavour(hid_t h_grp, hid_t h_grp_columns,
				      const struct engine* e) {

  io_write_attribute_s(h_grp, "Dust Model", "Following McKinnon et. al (2016)");

  /* Creatoe an array of grain names */
  const int grain_name_length = 32;
  char grain_names[grain_species_count][grain_name_length];
  for (int grain = 0; grain < grain_species_count; ++grain) {
    sprintf(grain_names[grain], "%s",
            dust_get_grain_name((enum grain_species)grain));
  }

  /* Add to the named columns */
  hsize_t dims[1] = {grain_species_count};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, grain_name_length);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "DustMassFractions", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, grain_names[0]);
  H5Dclose(dset);
  /* dset = H5Dcreate(h_grp_columns, "MetalDiffusionRates", type, space, */
  /*                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); */
  /* H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names[0]); */
  /* H5Dclose(dset); */

  H5Tclose(type);
  H5Sclose(space);
}

/**
 * @brief Create and write array for mapping elements to dust grains
 * @param h_grp The HDF5 group in which to write
 * @param e The #engine.
 */
INLINE static void dust_write_composition(hid_t h_grp,
					  const struct engine* e) {

  const struct dustevo_props* dp = e->dustevo;

  /* Create an array mapping grains to elements */
  int eldx;
  float grain_mapping[grain_species_count][chemistry_element_count] = {{0.}};
  for (int grain = 0; grain < grain_species_count; ++grain) {
  int elgrain = 0;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
      eldx = dp->grain_element_indices[grain][elgrain];
      if (eldx == elem) {
	grain_mapping[grain][elem] = dp->grain_element_mfrac[grain][elgrain];
        if (elgrain < dp->grain_element_count[grain]-1){
          elgrain += 1;
        }
      }
    }
  }


  /* Write out array */
  hsize_t dims[2] = {grain_species_count, chemistry_element_count};
  hid_t space = H5Screate_simple(2, dims, NULL);
  hid_t dset = H5Dcreate(h_grp, "GrainToElementMapping", H5T_NATIVE_FLOAT, 
			   space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) error("Error while creating grain to element map array");

  H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &grain_mapping[0][0]);
  io_write_attribute_s(dset, "Description", "Mass fraction of respective grain type (row)"
			                    "constituted by given element (column)");
  H5Dclose(dset);
  H5Sclose(space); 
}


#endif


#endif /* SWIFT_DUST_M16_PROPERTIES_H */
