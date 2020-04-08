/****************************************************************************
 * This file is part of CHIMES.
 * Copyright (c) 2020 Alexander Richings (alexander.j.richings@durham.ac.uk)
 *
 * CHIMES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***************************************************************************/

#ifdef CHIMES_ENABLE_GNU_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif

#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include "chimes_proto.h"
#include "chimes_vars.h"

// CHIMES data tables.
struct chimes_table_bins_struct
    chimes_table_bins; /*!< Structure containing table bins for all of the
                          rate tables. */
struct chimes_T_dependent_struct chimes_table_T_dependent; /*!< Structure
                                                              containing the
                                                              rates for the
                                                              T-dependent
                                                              reaction group. */
struct chimes_constant_struct
    chimes_table_constant; /*!< Structure containing the rates for the
                              constant reaction group. */
struct chimes_recombination_AB_struct
    chimes_table_recombination_AB; /*!< Structure containing the rates for the
                                      reacombination_AB group. */
struct chimes_grain_recombination_struct
    chimes_table_grain_recombination; /*!< Structure containing the rates from
                                         the grain_recombination group. */
struct chimes_cosmic_ray_struct
    chimes_table_cosmic_ray; /*!< Structure containing the rates for the
                                cosmic_ray group. */
struct chimes_CO_cosmic_ray_struct
    chimes_table_CO_cosmic_ray; /*!< Structure containing the rates for the
                                   CO_cosmic_ray group. */
struct chimes_H2_dust_formation_struct
    chimes_table_H2_dust_formation; /*!< Structure containing the rates for the
                                       H2_dust_formation group. */
struct chimes_H2_collis_dissoc_struct
    chimes_table_H2_collis_dissoc; /*!< Structure containing the rates for the
                                      H2_collis_dissoc group. */
struct chimes_photoion_fuv_struct chimes_table_photoion_fuv; /*!< Structure
                                                                containing the
                                                                rates for the
                                                                photoion_fuv
                                                                group. */
struct chimes_photoion_euv_struct chimes_table_photoion_euv; /*!< Structure
                                                                containing the
                                                                rates for the
                                                                photoion_euv
                                                                group. */
struct chimes_photoion_auger_fuv_struct
    chimes_table_photoion_auger_fuv; /*!< Structure containing the rates for the
                                        photoion_auger_fuv group. */
struct chimes_photoion_auger_euv_struct
    chimes_table_photoion_auger_euv; /*!< Structure containing the rates for the
                                        photoion_auger_euv group. */
struct chimes_photodissoc_group1_struct
    chimes_table_photodissoc_group1; /*!< Structure containing the rates for the
                                        photodissoc_group1 group. */
struct chimes_photodissoc_group2_struct
    chimes_table_photodissoc_group2; /*!< Structure containing the rates for the
                                        photodissoc_group2 group. */
struct chimes_H2_photodissoc_struct
    chimes_table_H2_photodissoc; /*!< Structure containing the rates for the
                                    H2_photodissoc group. */
struct chimes_CO_photodissoc_struct
    chimes_table_CO_photodissoc; /*!< Structure containing the rates for the
                                    CO_photodissoc group. */
struct chimes_spectra_struct
    chimes_table_spectra; /*!< Structure containing the spectra information. */
struct chimes_redshift_dependent_UVB_struct
    chimes_table_redshift_dependent_UVB; /*!< Structure containing information
                                            needed to interpolate the
                                            redshift-dependent UVB. */
struct chimes_cooling_struct
    chimes_table_cooling; /*!< Structure containing the cooling and heating
                             rates. */
struct chimes_eqm_abundances_struct
    chimes_table_eqm_abundances; /*!< Structure containing the equilibrium
                                    abundance tables. */

/**
 * @brief Allocate memory to eqm tables.
 *
 * Reads in the size of each dimension of the equilibrium
 * abundance tables and then allocates memory for the
 * arrays that will store these tables.
 *
 * @param filename String containing the path to the eqm abundance table file.
 * @param my_eqm_abundances The #chimes_eqm_abundances_struct struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void allocate_eqm_table_memory(
    char *filename, struct chimes_eqm_abundances_struct *my_eqm_abundances,
    struct globalVariables *myGlobalVars) {
  hid_t file_id, dataset;
  int i;
  float *T, *nH, *Z;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) {
    printf("CHIMES ERROR: unable to open equilibrium abundance data file: %s\n",
           filename);
    exit(EXIT_FAILURE);
  }

  /* Check that the number of species
   * in the equilibrium table matches
   * the size of the network. */
  int N_species;
  dataset = H5Dopen1(file_id, "TableBins/N_species");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_species);
  H5Dclose(dataset);

  if (N_species != myGlobalVars->totalNumberOfSpecies) {
    printf(
        "CHIMES ERROR: equilibrium table %s contains %d species, but we are "
        "using a network with %d species.\n",
        filename, N_species, myGlobalVars->totalNumberOfSpecies);
    exit(EXIT_FAILURE);
  }

  /* Read in size of each dimension */
  dataset = H5Dopen1(file_id, "TableBins/N_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(my_eqm_abundances->N_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_Densities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(my_eqm_abundances->N_Densities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_Metallicities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(my_eqm_abundances->N_Metallicities));
  H5Dclose(dataset);

  // Allocate memory to buffers
  T = (float *)malloc(my_eqm_abundances->N_Temperatures * sizeof(float));
  nH = (float *)malloc(my_eqm_abundances->N_Densities * sizeof(float));
  Z = (float *)malloc(my_eqm_abundances->N_Metallicities * sizeof(float));

  // Allocate memory to tables.
  my_eqm_abundances->Temperatures = (ChimesFloat *)malloc(
      my_eqm_abundances->N_Temperatures * sizeof(ChimesFloat));
  my_eqm_abundances->Densities = (ChimesFloat *)malloc(
      my_eqm_abundances->N_Densities * sizeof(ChimesFloat));
  my_eqm_abundances->Metallicities = (ChimesFloat *)malloc(
      my_eqm_abundances->N_Metallicities * sizeof(ChimesFloat));

  my_eqm_abundances->Abundances = (ChimesFloat *)malloc(
      myGlobalVars->totalNumberOfSpecies * my_eqm_abundances->N_Temperatures *
      my_eqm_abundances->N_Densities * my_eqm_abundances->N_Metallicities *
      sizeof(ChimesFloat));

  // Read in table bins
  dataset = H5Dopen1(file_id, "TableBins/Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, T);
  H5Dclose(dataset);

  for (i = 0; i < my_eqm_abundances->N_Temperatures; i++)
    my_eqm_abundances->Temperatures[i] = (ChimesFloat)T[i];

  dataset = H5Dopen1(file_id, "TableBins/Densities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nH);
  H5Dclose(dataset);

  for (i = 0; i < my_eqm_abundances->N_Densities; i++)
    my_eqm_abundances->Densities[i] = (ChimesFloat)nH[i];

  dataset = H5Dopen1(file_id, "TableBins/Metallicities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Z);
  H5Dclose(dataset);

  for (i = 0; i < my_eqm_abundances->N_Metallicities; i++)
    my_eqm_abundances->Metallicities[i] = (ChimesFloat)Z[i];

  H5Fclose(file_id);

  // Free buffer memory
  free(T);
  free(nH);
  free(Z);
}

/**
 * @brief Loads the eqm table.
 *
 * Reads in the equilibrium abundance table. Note that
 * the memory for arrays that will store these tables
 * needs to be allocated before this routine is called.
 *
 * @param filename String containing the path to the eqm abundance table file.
 * @param my_eqm_abundances The #chimes_eqm_abundances_struct struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void load_eqm_table(char *filename,
                    struct chimes_eqm_abundances_struct *my_eqm_abundances,
                    struct globalVariables *myGlobalVars) {
  hid_t file_id, dataset, dataspace_id, memspace_id;
  hsize_t dims[2];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int rank, i, j, k, l;
  float *array1D;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) {
    printf("CHIMES ERROR: unable to open equilibrium abundance data file: %s\n",
           filename);
    exit(EXIT_FAILURE);
  }

  // Allocate memory to buffers
  array1D = (float *)malloc(my_eqm_abundances->N_Temperatures * sizeof(float));

  // Read in the abundance array.
  dataset = H5Dopen1(file_id, "Abundances");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);

  dims[0] = my_eqm_abundances->N_Temperatures;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (i = 0; i < my_eqm_abundances->N_Densities; i++) {
    for (j = 0; j < my_eqm_abundances->N_Metallicities; j++) {
      for (k = 0; k < myGlobalVars->totalNumberOfSpecies; k++) {
        /* Note that if we switch off individual
         * elements, you will need to create a
         * new EqAbundances table that does NOT
         * include these excluded elements. */
        offset4D[0] = 0;
        offset4D[1] = i;
        offset4D[2] = j;
        offset4D[3] = k;
        count4D[0] = my_eqm_abundances->N_Temperatures;
        count4D[1] = 1;
        count4D[2] = 1;
        count4D[3] = 1;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL,
                            count4D, NULL);

        H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                H5P_DEFAULT, array1D);

        for (l = 0; l < my_eqm_abundances->N_Temperatures; l++)
          my_eqm_abundances->Abundances[chimes_flatten_index_4d(
              k, l, i, j, my_eqm_abundances->N_Temperatures,
              my_eqm_abundances->N_Densities,
              my_eqm_abundances->N_Metallicities)] = chimes_max((ChimesFloat)array1D[l], chimes_log10(CHIMES_FLT_MIN)); 
      }
    }
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  H5Fclose(file_id);

  // Free buffer memory
  free(array1D);
}

/**
 * @brief Determines whether to include a reaction.
 *
 * Compares the elements that are needed for a given
 * reaction to those that are to be included in the
 * network, as specified by the user. If all required
 * elements are included, this routine returns 1,
 * otherwise it returns 0.
 *
 * @param reaction_array Array of 9 integer flags specifying the elements needed
 * for the reaction.
 * @param network_array Array of 9 integer flags specifying the elements
 * included in the network.
 */
int compare_element_incl_arrays(int *reaction_array, int *network_array) {
  int i, include_reaction;

  include_reaction = 1;
  for (i = 0; i < 9; i++) {
    if ((reaction_array[i] == 1) && (network_array[i] == 0)) {
      include_reaction = 0;
      break;
    }
  }

  return include_reaction;
}

int compare_element_incl_arrays_flatten(int *reaction_array, int reaction_idx,
                                        int *network_array) {
  int i, include_reaction;

  include_reaction = 1;
  for (i = 0; i < 9; i++) {
    if ((reaction_array[i + (reaction_idx * 9)] == 1) &&
        (network_array[i] == 0)) {
      include_reaction = 0;
      break;
    }
  }

  return include_reaction;
}

/**
 * @brief Loads the main reaction data.
 *
 * Reads in the main reaction data tables. These define
 * all of the reactions in the network, and tabulate the
 * reaction rate coefficients.
 *
 * @param myGlobalVars The #globalVariables struct.
 */
void initialise_main_data(struct globalVariables *myGlobalVars) {
  char fname[500];
  int i, j, k, l, m, rank, N_reactions_all, incl_index;
  int base_index, base_incl_index, found_base;
  int *array_buffer_int;
  float *array_buffer_float;
  hid_t file_id, dataset, dataspace_id, memspace_id;
  hsize_t dims[2];
  hsize_t dims2D[2], count2D[2], offset2D[2];
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  hsize_t dims5D[5], count5D[5], offset5D[5];

  /* When we read the main data file, we will first
   * read each reaction group in to a 'master'
   * table structure. We will then copy them over
   * to global table structures, but including
   * only those reactions that are to be included
   * in the network, based on the element_incl flags. */
  struct chimes_T_dependent_struct chimes_master_T_dependent;
  struct chimes_constant_struct chimes_master_constant;
  struct chimes_recombination_AB_struct chimes_master_recombination_AB;
  struct chimes_grain_recombination_struct chimes_master_grain_recombination;
  struct chimes_cosmic_ray_struct chimes_master_cosmic_ray;
  struct chimes_CO_cosmic_ray_struct chimes_master_CO_cosmic_ray;
  struct chimes_photoion_fuv_struct chimes_master_photoion_fuv;
  struct chimes_photoion_euv_struct chimes_master_photoion_euv;
  struct chimes_photoion_auger_fuv_struct chimes_master_photoion_auger_fuv;
  struct chimes_photoion_auger_euv_struct chimes_master_photoion_auger_euv;
  struct chimes_photodissoc_group1_struct chimes_master_photodissoc_group1;
  struct chimes_photodissoc_group2_struct chimes_master_photodissoc_group2;
  struct chimes_CO_photodissoc_struct chimes_master_CO_photodissoc;
  struct chimes_cooling_struct chimes_master_cooling;

  ChimesFloat log_chimes_flt_min = chimes_log10(CHIMES_FLT_MIN);

  // Open main data file
  sprintf(fname, "%s", myGlobalVars->MainDataTablePath);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) {
    printf("CHIMES ERROR: unable to open main data file: %s\n", fname);
    exit(EXIT_FAILURE);
  }

  /****************
   ** Table bins **
   ****************/

  // Read size of each table bin
  dataset = H5Dopen1(file_id, "TableBins/N_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_Dust_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_Dust_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_Psi");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_Psi));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_secondary_cosmic_ray_xHII");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_secondary_cosmic_ray_xHII));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_H2self_column_densities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_H2self_column_densities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_b_turbulence");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_b_turbulence));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_COself_column_densities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_COself_column_densities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_H2CO_column_densities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_H2CO_column_densities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_mol_cool_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_mol_cool_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_2d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_2d_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_hiT_2d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_hiT_2d_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_2d_ElectronDensities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_2d_ElectronDensities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_4d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_4d_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_hiT_4d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_hiT_4d_Temperatures));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_4d_HIDensities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_4d_HIDensities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_4d_ElectronDensities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_4d_ElectronDensities));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "TableBins/N_cool_4d_HIIDensities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_bins.N_cool_4d_HIIDensities));
  H5Dclose(dataset);

  if ((myGlobalVars->element_included[0] == 1) &&
      (myGlobalVars->element_included[2] == 1)) {
    dataset = H5Dopen1(file_id, "TableBins/N_CO_cool_rot_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_bins.N_CO_cool_rot_ColumnDensities));
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "TableBins/N_CO_cool_vib_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_bins.N_CO_cool_vib_ColumnDensities));
    H5Dclose(dataset);
  }

  if (myGlobalVars->element_included[2] == 1) {
    dataset = H5Dopen1(file_id, "TableBins/N_H2O_cool_hiT_Temperatures");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_bins.N_H2O_cool_hiT_Temperatures));
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "TableBins/N_H2O_cool_lowT_Temperatures");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_bins.N_H2O_cool_lowT_Temperatures));
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "TableBins/N_H2O_cool_rot_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_bins.N_H2O_cool_rot_ColumnDensities));
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "TableBins/N_H2O_cool_vib_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_bins.N_H2O_cool_vib_ColumnDensities));
    H5Dclose(dataset);
  }

  // Allocate memory
  chimes_table_bins.Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.Dust_Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_Dust_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.Psi =
      (ChimesFloat *)malloc(chimes_table_bins.N_Psi * sizeof(ChimesFloat));
  chimes_table_bins.secondary_cosmic_ray_xHII = (ChimesFloat *)malloc(
      chimes_table_bins.N_secondary_cosmic_ray_xHII * sizeof(ChimesFloat));
  chimes_table_bins.H2self_column_densities = (ChimesFloat *)malloc(
      chimes_table_bins.N_H2self_column_densities * sizeof(ChimesFloat));
  chimes_table_bins.b_turbulence = (ChimesFloat *)malloc(
      chimes_table_bins.N_b_turbulence * sizeof(ChimesFloat));
  chimes_table_bins.COself_column_densities = (ChimesFloat *)malloc(
      chimes_table_bins.N_COself_column_densities * sizeof(ChimesFloat));
  chimes_table_bins.H2CO_column_densities = (ChimesFloat *)malloc(
      chimes_table_bins.N_H2CO_column_densities * sizeof(ChimesFloat));
  chimes_table_bins.mol_cool_Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.cool_2d_Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_2d_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.cool_hiT_2d_Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_hiT_2d_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.cool_2d_ElectronDensities = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_2d_ElectronDensities * sizeof(ChimesFloat));
  chimes_table_bins.cool_4d_Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_4d_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.cool_hiT_4d_Temperatures = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_hiT_4d_Temperatures * sizeof(ChimesFloat));
  chimes_table_bins.cool_4d_HIDensities = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_4d_HIDensities * sizeof(ChimesFloat));
  chimes_table_bins.cool_4d_ElectronDensities = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_4d_ElectronDensities * sizeof(ChimesFloat));
  chimes_table_bins.cool_4d_HIIDensities = (ChimesFloat *)malloc(
      chimes_table_bins.N_cool_4d_HIIDensities * sizeof(ChimesFloat));

  if ((myGlobalVars->element_included[0] == 1) &&
      (myGlobalVars->element_included[2] == 1)) {
    chimes_table_bins.CO_cool_rot_ColumnDensities = (ChimesFloat *)malloc(
        chimes_table_bins.N_CO_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_bins.CO_cool_vib_ColumnDensities = (ChimesFloat *)malloc(
        chimes_table_bins.N_CO_cool_vib_ColumnDensities * sizeof(ChimesFloat));
  }

  if (myGlobalVars->element_included[2] == 1) {
    chimes_table_bins.H2O_cool_hiT_Temperatures = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures * sizeof(ChimesFloat));
    chimes_table_bins.H2O_cool_lowT_Temperatures = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures * sizeof(ChimesFloat));
    chimes_table_bins.H2O_cool_rot_ColumnDensities = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_bins.H2O_cool_vib_ColumnDensities = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_vib_ColumnDensities * sizeof(ChimesFloat));
  }

  // Read data arrays
  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
    chimes_table_bins.Temperatures[i] = (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_Dust_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/Dust_Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Dust_Temperatures; i++)
    chimes_table_bins.Dust_Temperatures[i] = (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(chimes_table_bins.N_Psi * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/Psi");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Psi; i++)
    chimes_table_bins.Psi[i] = (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_secondary_cosmic_ray_xHII * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/secondary_cosmic_ray_xHII");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_secondary_cosmic_ray_xHII; i++)
    chimes_table_bins.secondary_cosmic_ray_xHII[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_H2self_column_densities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/H2self_column_densities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_H2self_column_densities; i++)
    chimes_table_bins.H2self_column_densities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_b_turbulence * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/b_turbulence");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_b_turbulence; i++)
    chimes_table_bins.b_turbulence[i] = (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_COself_column_densities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/COself_column_densities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_COself_column_densities; i++)
    chimes_table_bins.COself_column_densities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_H2CO_column_densities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/H2CO_column_densities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_H2CO_column_densities; i++)
    chimes_table_bins.H2CO_column_densities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_cool_2d_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_2d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_2d_Temperatures; i++)
    chimes_table_bins.cool_2d_Temperatures[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_cool_2d_ElectronDensities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_2d_ElectronDensities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_2d_ElectronDensities; i++)
    chimes_table_bins.cool_2d_ElectronDensities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_cool_hiT_2d_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_hiT_2d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_hiT_2d_Temperatures; i++)
    chimes_table_bins.cool_hiT_2d_Temperatures[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_cool_4d_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_4d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_4d_Temperatures; i++)
    chimes_table_bins.cool_4d_Temperatures[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_cool_hiT_4d_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_hiT_4d_Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_hiT_4d_Temperatures; i++)
    chimes_table_bins.cool_hiT_4d_Temperatures[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_cool_4d_HIDensities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_4d_HIDensities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_4d_HIDensities; i++)
    chimes_table_bins.cool_4d_HIDensities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_cool_4d_ElectronDensities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_4d_ElectronDensities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_4d_ElectronDensities; i++)
    chimes_table_bins.cool_4d_ElectronDensities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_cool_4d_HIIDensities * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/cool_4d_HIIDensities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_cool_4d_HIIDensities; i++)
    chimes_table_bins.cool_4d_HIIDensities[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "TableBins/mol_cool_Temperatures");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_table_bins.mol_cool_Temperatures[i] =
        (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  if ((myGlobalVars->element_included[0] == 1) &&
      (myGlobalVars->element_included[2] == 1)) {
    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_CO_cool_rot_ColumnDensities * sizeof(float));
    dataset = H5Dopen1(file_id, "TableBins/CO_cool_rot_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_CO_cool_rot_ColumnDensities; i++)
      chimes_table_bins.CO_cool_rot_ColumnDensities[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_CO_cool_vib_ColumnDensities * sizeof(float));
    dataset = H5Dopen1(file_id, "TableBins/CO_cool_vib_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_CO_cool_vib_ColumnDensities; i++)
      chimes_table_bins.CO_cool_vib_ColumnDensities[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);
  }

  if (myGlobalVars->element_included[2] == 1) {
    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures * sizeof(float));
    dataset = H5Dopen1(file_id, "TableBins/H2O_cool_hiT_Temperatures");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_hiT_Temperatures; i++)
      chimes_table_bins.H2O_cool_hiT_Temperatures[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures * sizeof(float));
    dataset = H5Dopen1(file_id, "TableBins/H2O_cool_lowT_Temperatures");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
      chimes_table_bins.H2O_cool_lowT_Temperatures[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(float));
    dataset = H5Dopen1(file_id, "TableBins/H2O_cool_rot_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; i++)
      chimes_table_bins.H2O_cool_rot_ColumnDensities[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_H2O_cool_vib_ColumnDensities * sizeof(float));
    dataset = H5Dopen1(file_id, "TableBins/H2O_cool_vib_ColumnDensities");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_vib_ColumnDensities; i++)
      chimes_table_bins.H2O_cool_vib_ColumnDensities[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);
  }

  /***************************
   ** T_dependent reactions **
   ***************************/

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "T_dependent/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_T_dependent.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_master_T_dependent.N_reactions[1];

  // Allocate memory in table structure
  chimes_master_T_dependent.reactants =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_master_T_dependent.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_master_T_dependent.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_master_T_dependent.rates = (ChimesFloat *)malloc(
      N_reactions_all * chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_master_T_dependent.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));

  array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "T_dependent/reactants");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 3; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_T_dependent.reactants[chimes_flatten_index_2d(i, j, 3)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "T_dependent/products");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 3; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_T_dependent.products[chimes_flatten_index_2d(i, j, 3)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "T_dependent/rates");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_T_dependent.rates[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "T_dependent/element_incl");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 9; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_T_dependent.element_incl[chimes_flatten_index_2d(i, j, 9)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "T_dependent/molecular_flag");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_T_dependent.molecular_flag);
  H5Dclose(dataset);

  dataset =
      H5Dopen1(file_id, "T_dependent/H2_collis_dissoc_heating_reaction_index");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_master_T_dependent.H2_collis_dissoc_heating_reaction_index));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "T_dependent/H2_form_heating_reaction_index");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_master_T_dependent.H2_form_heating_reaction_index));
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_float);
  free(array_buffer_int);

  /************************
   ** constant reactions **
   ************************/

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "constant/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_constant.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_master_constant.N_reactions[1];

  // Allocate memory in table structure
  chimes_master_constant.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_master_constant.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_master_constant.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_master_constant.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_constant.rates =
      (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));

  array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "constant/reactants");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 2; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_constant.reactants[chimes_flatten_index_2d(i, j, 2)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "constant/products");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 3; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_constant.products[chimes_flatten_index_2d(i, j, 3)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "constant/element_incl");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 9; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_constant.element_incl[chimes_flatten_index_2d(i, j, 9)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "constant/molecular_flag");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_constant.molecular_flag);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "constant/rates");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < N_reactions_all; i++)
    chimes_master_constant.rates[i] = (ChimesFloat)array_buffer_float[i];

  dataset = H5Dopen1(file_id, "constant/H2_form_heating_reaction_index");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_master_constant.H2_form_heating_reaction_index));
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_float);
  free(array_buffer_int);

  /********************************
   ** recombination_AB reactions **
   ********************************/

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "recombination_AB/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_recombination_AB.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_master_recombination_AB.N_reactions[1];

  // Allocate memory in table structure
  chimes_master_recombination_AB.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_master_recombination_AB.products =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_recombination_AB.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_master_recombination_AB.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_recombination_AB.rates = (ChimesFloat *)malloc(
      N_reactions_all * 2 * chimes_table_bins.N_Temperatures *
      sizeof(ChimesFloat));

  array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "recombination_AB/reactants");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 2; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_recombination_AB
          .reactants[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "recombination_AB/rates_caseA");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_recombination_AB.rates[chimes_flatten_index_3d(
          i, 0, j, 2, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "recombination_AB/rates_caseB");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_recombination_AB.rates[chimes_flatten_index_3d(
          i, 1, j, 2, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "recombination_AB/element_incl");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 9; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_recombination_AB
          .element_incl[chimes_flatten_index_2d(i, j, 9)] = array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "recombination_AB/products");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_recombination_AB.products);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "recombination_AB/molecular_flag");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_recombination_AB.molecular_flag);
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_float);
  free(array_buffer_int);

  /********************************
   ** grain_recombination reactions **
   ********************************/

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "grain_recombination/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_grain_recombination.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_master_grain_recombination.N_reactions[1];

  // Allocate memory in table structure
  chimes_master_grain_recombination.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_master_grain_recombination.products =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_grain_recombination.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_master_grain_recombination.rates =
      (ChimesFloat *)malloc(N_reactions_all * chimes_table_bins.N_Temperatures *
                            chimes_table_bins.N_Psi * sizeof(ChimesFloat));

  array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "grain_recombination/reactants");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 2; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_grain_recombination
          .reactants[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "grain_recombination/rates");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
  dims3D[0] = N_reactions_all;
  dims3D[1] = 1;
  dims3D[2] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims3D, NULL);

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++) {
    for (j = 0; j < chimes_table_bins.N_Psi; j++) {
      offset3D[0] = 0;
      offset3D[1] = i;
      offset3D[2] = j;
      count3D[0] = N_reactions_all;
      count3D[1] = 1;
      count3D[2] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (k = 0; k < N_reactions_all; k++)
        chimes_master_grain_recombination.rates[chimes_flatten_index_3d(
            k, i, j, chimes_table_bins.N_Temperatures,
            chimes_table_bins.N_Psi)] =
            chimes_max((ChimesFloat)array_buffer_float[k], log_chimes_flt_min);
    }
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "grain_recombination/element_incl");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 9; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_grain_recombination
          .element_incl[chimes_flatten_index_2d(i, j, 9)] = array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "grain_recombination/products");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_grain_recombination.products);
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_float);
  free(array_buffer_int);

  /************************
   ** cosmic_ray reactions **
   ************************/

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "cosmic_ray/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cosmic_ray.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_master_cosmic_ray.N_reactions[1];

  // Allocate memory in table structure
  chimes_master_cosmic_ray.reactants =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_cosmic_ray.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_master_cosmic_ray.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_master_cosmic_ray.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_cosmic_ray.rates =
      (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
  chimes_master_cosmic_ray.secondary_base_reaction =
      (int *)malloc(2 * sizeof(int));
  chimes_master_cosmic_ray.secondary_ratio = (ChimesFloat *)malloc(
      2 * chimes_table_bins.N_secondary_cosmic_ray_xHII * sizeof(ChimesFloat));

  array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "cosmic_ray/reactants");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cosmic_ray.reactants);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cosmic_ray/products");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 3; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_cosmic_ray.products[chimes_flatten_index_2d(i, j, 3)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cosmic_ray/element_incl");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 9; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_cosmic_ray.element_incl[chimes_flatten_index_2d(i, j, 9)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cosmic_ray/molecular_flag");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cosmic_ray.molecular_flag);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cosmic_ray/rates");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < N_reactions_all; i++)
    chimes_master_cosmic_ray.rates[i] = (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  dataset = H5Dopen1(file_id, "cosmic_ray/secondary_base_reaction");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cosmic_ray.secondary_base_reaction);
  H5Dclose(dataset);

  array_buffer_float = (float *)malloc(2 * sizeof(float));

  dataset = H5Dopen1(file_id, "cosmic_ray/secondary_ratio");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = 2;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_secondary_cosmic_ray_xHII; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = 2;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < 2; i++)
      chimes_master_cosmic_ray.secondary_ratio[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_secondary_cosmic_ray_xHII)] =
          (ChimesFloat)array_buffer_float[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_int);
  free(array_buffer_float);

  /*****************************
   ** CO_cosmic_ray reactions **
   *****************************/

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "CO_cosmic_ray/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_CO_cosmic_ray.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_master_CO_cosmic_ray.N_reactions[1];

  // Allocate memory in table structure
  chimes_master_CO_cosmic_ray.reactants =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_master_CO_cosmic_ray.products =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_master_CO_cosmic_ray.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_master_CO_cosmic_ray.rates = (ChimesFloat *)malloc(
      N_reactions_all * chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));

  array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "CO_cosmic_ray/reactants");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_CO_cosmic_ray.reactants);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "CO_cosmic_ray/products");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 2; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_CO_cosmic_ray.products[chimes_flatten_index_2d(i, j, 2)] =
          array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "CO_cosmic_ray/element_incl");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 9; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_CO_cosmic_ray
          .element_incl[chimes_flatten_index_2d(i, j, 9)] = array_buffer_int[i];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "CO_cosmic_ray/rates");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_CO_cosmic_ray.rates[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_int);
  free(array_buffer_float);

  if (myGlobalVars->N_spectra > 0) {
    /***************************
     ** photoion_fuv reactions **
     ***************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "photoion_fuv/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_fuv.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_photoion_fuv.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_photoion_fuv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_photoion_fuv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_photoion_fuv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_photoion_fuv.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "photoion_fuv/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_fuv.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_fuv/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_fuv.products[chimes_flatten_index_2d(i, j, 2)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_fuv/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_fuv
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_fuv/gamma");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_photoion_fuv.gamma[i] = (ChimesFloat)array_buffer_float[i];

    // Free buffers
    free(array_buffer_float);
    free(array_buffer_int);

    /***************************
     ** photoion_euv reactions **
     ***************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "photoion_euv/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_euv.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_photoion_euv.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_photoion_euv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_photoion_euv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_photoion_euv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_photoion_euv.E_thresh =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_master_photoion_euv.molecular_flag =
        (int *)malloc(N_reactions_all * sizeof(int));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "photoion_euv/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_euv.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_euv/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_euv.products[chimes_flatten_index_2d(i, j, 2)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_euv/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_euv
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_euv/molecular_flag");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_euv.molecular_flag);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_euv/E_thresh");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_photoion_euv.E_thresh[i] =
          (ChimesFloat)array_buffer_float[i];

    // Free buffers
    free(array_buffer_int);
    free(array_buffer_float);

    /**********************************
     ** photoion_auger_fuv reactions **
     **********************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "photoion_auger_fuv/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_auger_fuv.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_photoion_auger_fuv.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_photoion_auger_fuv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_photoion_auger_fuv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_photoion_auger_fuv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_photoion_auger_fuv.base_reaction =
        (int *)malloc(N_reactions_all * sizeof(int));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "photoion_auger_fuv/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_auger_fuv.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_auger_fuv/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_auger_fuv
            .products[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_auger_fuv/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_auger_fuv
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_auger_fuv/base_reaction");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_auger_fuv.base_reaction);
    H5Dclose(dataset);

    /**********************************
     ** photoion_auger_euv reactions **
     **********************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "photoion_auger_euv/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_auger_euv.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_photoion_auger_euv.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_photoion_auger_euv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_photoion_auger_euv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_photoion_auger_euv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_photoion_auger_euv.base_reaction =
        (int *)malloc(N_reactions_all * sizeof(int));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "photoion_auger_euv/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_auger_euv.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_auger_euv/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_auger_euv
            .products[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_auger_euv/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photoion_auger_euv
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_auger_euv/base_reaction");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photoion_auger_euv.base_reaction);
    H5Dclose(dataset);

    // Free buffers
    free(array_buffer_float);
    free(array_buffer_int);

    /**********************************
     ** photodissoc_group1 reactions **
     **********************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "photodissoc_group1/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photodissoc_group1.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_photodissoc_group1.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_photodissoc_group1.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_photodissoc_group1.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_photodissoc_group1.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_photodissoc_group1.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_master_photodissoc_group1.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_master_photodissoc_group1.molecular_flag =
        (int *)malloc(N_reactions_all * sizeof(int));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "photodissoc_group1/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photodissoc_group1.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photodissoc_group1/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photodissoc_group1
            .products[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photodissoc_group1/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photodissoc_group1
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photodissoc_group1/gamma");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_photodissoc_group1.gamma[i] =
          (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "photodissoc_group1/rates");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_photodissoc_group1.rates[i] =
          (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "photodissoc_group1/molecular_flag");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photodissoc_group1.molecular_flag);
    H5Dclose(dataset);

    // Free buffers
    free(array_buffer_float);
    free(array_buffer_int);

    /**********************************
     ** photodissoc_group2 reactions **
     **********************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "photodissoc_group2/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photodissoc_group2.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_photodissoc_group2.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_photodissoc_group2.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_photodissoc_group2.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_photodissoc_group2.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_photodissoc_group2.gamma_coeff =
        (ChimesFloat *)malloc(3 * sizeof(ChimesFloat));
    chimes_master_photodissoc_group2.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "photodissoc_group2/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_photodissoc_group2.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photodissoc_group2/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photodissoc_group2
            .products[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photodissoc_group2/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_photodissoc_group2
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photodissoc_group2/rates");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_photodissoc_group2.rates[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(3 * sizeof(float));
    dataset = H5Dopen1(file_id, "photodissoc_group2/gamma_coeff");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < 3; i++)
      chimes_master_photodissoc_group2.gamma_coeff[i] =
          (ChimesFloat)array_buffer_float[i];

    // Free buffers
    free(array_buffer_float);
    free(array_buffer_int);

    /*************************************
     ** H2_photodissoc reactions.       **
     ** All reactions involve H only,   **
     ** so we can read this group       **
     ** straight into the global table. **
     *************************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "H2_photodissoc/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_table_H2_photodissoc.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_table_H2_photodissoc.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_H2_photodissoc.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_H2_photodissoc.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_H2_photodissoc.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_H2_photodissoc.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_H2_photodissoc.self_shielding = (ChimesFloat *)malloc(
        N_reactions_all * chimes_table_bins.N_Temperatures *
        chimes_table_bins.N_H2self_column_densities *
        chimes_table_bins.N_b_turbulence * sizeof(ChimesFloat));
    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "H2_photodissoc/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_int);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_photodissoc.reactants[i] =
          myGlobalVars->speciesIndices[array_buffer_int[i]];

    dataset = H5Dopen1(file_id, "H2_photodissoc/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_table_H2_photodissoc.products[chimes_flatten_index_2d(i, j, 2)] =
            myGlobalVars->speciesIndices[array_buffer_int[i]];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "H2_photodissoc/gamma");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_photodissoc.gamma[i] = (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "H2_photodissoc/rates");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_photodissoc.rates[i] = (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "H2_photodissoc/self_shielding");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
    dims4D[0] = N_reactions_all;
    dims4D[1] = 1;
    dims4D[2] = 1;
    dims4D[3] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims4D, NULL);

    for (i = 0; i < chimes_table_bins.N_Temperatures; i++) {
      for (j = 0; j < chimes_table_bins.N_H2self_column_densities; j++) {
        for (k = 0; k < chimes_table_bins.N_b_turbulence; k++) {
          offset4D[0] = 0;
          offset4D[1] = i;
          offset4D[2] = j;
          offset4D[3] = k;
          count4D[0] = N_reactions_all;
          count4D[1] = 1;
          count4D[2] = 1;
          count4D[3] = 1;
          H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL,
                              count4D, NULL);

          H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                  H5P_DEFAULT, array_buffer_float);

          for (l = 0; l < N_reactions_all; l++)
            chimes_table_H2_photodissoc.self_shielding[chimes_flatten_index_4d(
                l, i, j, k, chimes_table_bins.N_Temperatures,
                chimes_table_bins.N_H2self_column_densities,
                chimes_table_bins.N_b_turbulence)] =
                (ChimesFloat)array_buffer_float[l];
        }
      }
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    // Free buffers
    free(array_buffer_float);
    free(array_buffer_int);

    /*******************************
     ** CO_photodissoc reactions. **
     *******************************/

    // Read number of reactions in group
    dataset = H5Dopen1(file_id, "CO_photodissoc/N_reactions");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_CO_photodissoc.N_reactions);
    H5Dclose(dataset);

    N_reactions_all = chimes_master_CO_photodissoc.N_reactions[1];

    // Allocate memory in table structure
    chimes_master_CO_photodissoc.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_master_CO_photodissoc.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_master_CO_photodissoc.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_master_CO_photodissoc.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_master_CO_photodissoc.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_master_CO_photodissoc.self_shielding = (ChimesFloat *)malloc(
        N_reactions_all * chimes_table_bins.N_COself_column_densities *
        chimes_table_bins.N_H2CO_column_densities * sizeof(ChimesFloat));

    array_buffer_int = (int *)malloc(N_reactions_all * sizeof(int));
    array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

    // Read data arrays
    dataset = H5Dopen1(file_id, "CO_photodissoc/reactants");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            chimes_master_CO_photodissoc.reactants);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "CO_photodissoc/products");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 2; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_CO_photodissoc
            .products[chimes_flatten_index_2d(i, j, 2)] = array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "CO_photodissoc/element_incl");

    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = N_reactions_all;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < 9; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = N_reactions_all;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_int);

      for (i = 0; i < N_reactions_all; i++)
        chimes_master_CO_photodissoc
            .element_incl[chimes_flatten_index_2d(i, j, 9)] =
            array_buffer_int[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "CO_photodissoc/gamma");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_CO_photodissoc.gamma[i] =
          (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "CO_photodissoc/rates");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < N_reactions_all; i++)
      chimes_master_CO_photodissoc.rates[i] =
          (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "CO_photodissoc/self_shielding");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
    dims3D[0] = N_reactions_all;
    dims3D[1] = 1;
    dims3D[2] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims3D, NULL);

    for (i = 0; i < chimes_table_bins.N_COself_column_densities; i++) {
      for (j = 0; j < chimes_table_bins.N_H2CO_column_densities; j++) {
        offset3D[0] = 0;
        offset3D[1] = i;
        offset3D[2] = j;
        count3D[0] = N_reactions_all;
        count3D[1] = 1;
        count3D[2] = 1;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL,
                            count3D, NULL);

        H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                H5P_DEFAULT, array_buffer_float);

        for (l = 0; l < N_reactions_all; l++)
          chimes_master_CO_photodissoc.self_shielding[chimes_flatten_index_3d(
              l, i, j, chimes_table_bins.N_COself_column_densities,
              chimes_table_bins.N_H2CO_column_densities)] =
              (ChimesFloat)array_buffer_float[l];
      }
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    // Free buffers
    free(array_buffer_float);
    free(array_buffer_int);
  }

  /*********************************
   ** H2_dust_formation reactions **
   *********************************/

  /* This group has only one reactions, and
   * it only contains hydrogen, so we can
   * read this straight into the global
   * table, no need to go through a master
   * table first. */

  // Allocate memory in table structure
  chimes_table_H2_dust_formation.reactants = (int *)malloc(2 * sizeof(int));
  chimes_table_H2_dust_formation.products = (int *)malloc(sizeof(int));
  chimes_table_H2_dust_formation.rates = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * chimes_table_bins.N_Dust_Temperatures *
      sizeof(ChimesFloat));

  array_buffer_int = (int *)malloc(3 * sizeof(int));
  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_Temperatures * sizeof(float));

  // Read data arrays
  dataset = H5Dopen1(file_id, "H2_dust_formation/reactants");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_int);
  H5Dclose(dataset);

  for (i = 0; i < 2; i++)
    chimes_table_H2_dust_formation.reactants[i] =
        myGlobalVars->speciesIndices[array_buffer_int[i]];

  dataset = H5Dopen1(file_id, "H2_dust_formation/products");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_int);
  H5Dclose(dataset);

  chimes_table_H2_dust_formation.products[0] =
      myGlobalVars->speciesIndices[array_buffer_int[0]];

  dataset = H5Dopen1(file_id, "H2_dust_formation/rates");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = chimes_table_bins.N_Temperatures;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Dust_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = chimes_table_bins.N_Temperatures;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
      chimes_table_H2_dust_formation.rates[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Dust_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  // Free buffers
  free(array_buffer_float);
  free(array_buffer_int);

  /*********************************
   ** H2_collis_dissoc reactions **
   *********************************/

  /* The reactions in this group only involve
   * H and He, so no need to first read them
   * into a 'master' table - we can just read
   * straight into the global table. Also,
   * they all involve molecules (i.e. molecular
   * flag of 1). */

  // Read number of reactions in group
  dataset = H5Dopen1(file_id, "H2_collis_dissoc/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_table_H2_collis_dissoc.N_reactions);
  H5Dclose(dataset);

  N_reactions_all = chimes_table_H2_collis_dissoc.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_H2_collis_dissoc.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_table_H2_collis_dissoc.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_table_H2_collis_dissoc.k0 = (ChimesFloat *)malloc(
      N_reactions_all * chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_table_H2_collis_dissoc.kLTE = (ChimesFloat *)malloc(
      N_reactions_all * chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_table_H2_collis_dissoc.critical_density_H = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_table_H2_collis_dissoc.critical_density_H2 = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_table_H2_collis_dissoc.critical_density_He = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));

  // Read data arrays
  array_buffer_int = (int *)malloc(2 * sizeof(int));
  dataset = H5Dopen1(file_id, "H2_collis_dissoc/reactants");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 2; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_collis_dissoc
          .reactants[chimes_flatten_index_2d(i, j, 2)] =
          myGlobalVars->speciesIndices[array_buffer_int[i]];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);
  free(array_buffer_int);

  array_buffer_int = (int *)malloc(3 * sizeof(int));
  dataset = H5Dopen1(file_id, "H2_collis_dissoc/products");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < 3; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_int);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_collis_dissoc.products[chimes_flatten_index_2d(i, j, 3)] =
          myGlobalVars->speciesIndices[array_buffer_int[i]];
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);
  free(array_buffer_int);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_Temperatures * sizeof(float));
  dataset = H5Dopen1(file_id, "H2_collis_dissoc/critical_density_H");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
    chimes_table_H2_collis_dissoc.critical_density_H[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "H2_collis_dissoc/critical_density_H2");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
    chimes_table_H2_collis_dissoc.critical_density_H2[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "H2_collis_dissoc/critical_density_He");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
    chimes_table_H2_collis_dissoc.critical_density_He[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));
  dataset = H5Dopen1(file_id, "H2_collis_dissoc/k0");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_collis_dissoc
          .k0[chimes_flatten_index_2d(i, j, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "H2_collis_dissoc/kLTE");

  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = N_reactions_all;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = N_reactions_all;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < N_reactions_all; i++)
      chimes_table_H2_collis_dissoc.kLTE[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  free(array_buffer_float);

  dataset = H5Dopen1(file_id, "H2_collis_dissoc/Heating_reaction_index");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_table_H2_collis_dissoc.Heating_reaction_index));
  H5Dclose(dataset);

  /*************
   ** cooling **
   *************/

  // Read number of coolants
  dataset = H5Dopen1(file_id, "cooling/N_coolants");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_master_cooling.N_coolants));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/N_coolants_2d");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_master_cooling.N_coolants_2d));
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/N_coolants_4d");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(chimes_master_cooling.N_coolants_4d));
  H5Dclose(dataset);

  // Allocate memory in table structure
  chimes_master_cooling.coolants =
      (int *)malloc(chimes_master_cooling.N_coolants * sizeof(int));
  chimes_master_cooling.coolants_2d =
      (int *)malloc(chimes_master_cooling.N_coolants_2d * sizeof(int));
  chimes_master_cooling.coolants_4d =
      (int *)malloc(chimes_master_cooling.N_coolants_4d * sizeof(int));

  chimes_master_cooling.rates = (ChimesFloat *)malloc(
      chimes_master_cooling.N_coolants * chimes_table_bins.N_Temperatures *
      sizeof(ChimesFloat));

  chimes_master_cooling.rates_2d = (ChimesFloat *)malloc(
      chimes_master_cooling.N_coolants_2d *
      chimes_table_bins.N_cool_2d_Temperatures *
      chimes_table_bins.N_cool_2d_ElectronDensities * sizeof(ChimesFloat));
  chimes_master_cooling.rates_hiT_2d = (ChimesFloat *)malloc(
      chimes_master_cooling.N_coolants_2d *
      chimes_table_bins.N_cool_hiT_2d_Temperatures * sizeof(ChimesFloat));

  chimes_master_cooling.rates_4d = (ChimesFloat *)malloc(
      chimes_master_cooling.N_coolants_4d *
      chimes_table_bins.N_cool_4d_Temperatures *
      chimes_table_bins.N_cool_4d_HIDensities *
      chimes_table_bins.N_cool_4d_ElectronDensities *
      chimes_table_bins.N_cool_4d_HIIDensities * sizeof(ChimesFloat));
  chimes_master_cooling.rates_hiT_4d = (ChimesFloat *)malloc(
      chimes_master_cooling.N_coolants_4d *
      chimes_table_bins.N_cool_hiT_4d_Temperatures * sizeof(ChimesFloat));

  chimes_master_cooling.photoelectric_heating =
      (ChimesFloat *)malloc(chimes_table_bins.N_Temperatures *
                            chimes_table_bins.N_Psi * sizeof(ChimesFloat));
  chimes_master_cooling.grain_recombination =
      (ChimesFloat *)malloc(chimes_table_bins.N_Temperatures *
                            chimes_table_bins.N_Psi * sizeof(ChimesFloat));
  chimes_master_cooling.gas_grain_transfer = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));

  chimes_master_cooling.H2_cool_lowDens_H2 = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_master_cooling.H2_cool_lowDens_HI = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_master_cooling.H2_cool_lowDens_HII = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_master_cooling.H2_cool_lowDens_HeI = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_master_cooling.H2_cool_lowDens_elec = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_master_cooling.H2_cool_LTE = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));

  /* CO_cool and H2O_cool arrays can be read
   * directly into chimes_table_cooling, if
   * the required elements are present. */
  if ((myGlobalVars->element_included[0] == 1) &&
      (myGlobalVars->element_included[2] == 1)) {
    chimes_table_cooling.CO_cool_rot_L0 = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
    chimes_table_cooling.CO_cool_vib_L0 = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
    chimes_table_cooling.CO_cool_rot_Llte = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures *
        chimes_table_bins.N_CO_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.CO_cool_rot_nhalf = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures *
        chimes_table_bins.N_CO_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.CO_cool_rot_a = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures *
        chimes_table_bins.N_CO_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.CO_cool_vib_Llte = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures *
        chimes_table_bins.N_CO_cool_vib_ColumnDensities * sizeof(ChimesFloat));
  }

  if (myGlobalVars->element_included[2] == 1) {
    chimes_table_cooling.H2O_cool_rot_hiT_L0 = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures * sizeof(ChimesFloat));
    chimes_table_cooling.H2O_cool_rot_hiT_Llte = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2O_cool_rot_hiT_nhalf = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2O_cool_rot_hiT_a = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));

    chimes_table_cooling.H2Oortho_cool_rot_lowT_L0 = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures * sizeof(ChimesFloat));
    chimes_table_cooling.H2Opara_cool_rot_lowT_L0 = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures * sizeof(ChimesFloat));
    chimes_table_cooling.H2Oortho_cool_rot_lowT_Llte = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2Oortho_cool_rot_lowT_nhalf = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2Oortho_cool_rot_lowT_a = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2Opara_cool_rot_lowT_Llte = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2Opara_cool_rot_lowT_nhalf = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));
    chimes_table_cooling.H2Opara_cool_rot_lowT_a = (ChimesFloat *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures *
        chimes_table_bins.N_H2O_cool_rot_ColumnDensities * sizeof(ChimesFloat));

    chimes_table_cooling.H2O_cool_vib_L0 = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
    chimes_table_cooling.H2O_cool_vib_Llte = (ChimesFloat *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures *
        chimes_table_bins.N_H2O_cool_vib_ColumnDensities * sizeof(ChimesFloat));
  }

  // Read data arrays
  dataset = H5Dopen1(file_id, "cooling/coolants");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cooling.coolants);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/coolants_2d");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cooling.coolants_2d);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/coolants_4d");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          chimes_master_cooling.coolants_4d);
  H5Dclose(dataset);

  array_buffer_float =
      (float *)malloc(chimes_master_cooling.N_coolants * sizeof(float));

  dataset = H5Dopen1(file_id, "cooling/rates");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = chimes_master_cooling.N_coolants;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = chimes_master_cooling.N_coolants;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < chimes_master_cooling.N_coolants; i++)
      chimes_master_cooling.rates[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_master_cooling.N_coolants_2d * sizeof(float));

  dataset = H5Dopen1(file_id, "cooling/rates_2d");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
  dims3D[0] = chimes_master_cooling.N_coolants_2d;
  dims3D[1] = 1;
  dims3D[2] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims3D, NULL);

  for (j = 0; j < chimes_table_bins.N_cool_2d_Temperatures; j++) {
    for (k = 0; k < chimes_table_bins.N_cool_2d_ElectronDensities; k++) {
      offset3D[0] = 0;
      offset3D[1] = j;
      offset3D[2] = k;
      count3D[0] = chimes_master_cooling.N_coolants_2d;
      count3D[1] = 1;
      count3D[2] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_master_cooling.N_coolants_2d; i++)
        chimes_master_cooling.rates_2d[chimes_flatten_index_3d(
            i, j, k, chimes_table_bins.N_cool_2d_Temperatures,
            chimes_table_bins.N_cool_2d_ElectronDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/rates_hiT_2d");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = chimes_master_cooling.N_coolants_2d;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_cool_hiT_2d_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = chimes_master_cooling.N_coolants_2d;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < chimes_master_cooling.N_coolants_2d; i++)
      chimes_master_cooling.rates_hiT_2d[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_cool_hiT_2d_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_master_cooling.N_coolants_4d * sizeof(float));

  dataset = H5Dopen1(file_id, "cooling/rates_4d");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims5D, NULL);
  dims5D[0] = chimes_master_cooling.N_coolants_4d;
  dims5D[1] = 1;
  dims5D[2] = 1;
  dims5D[3] = 1;
  dims5D[4] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims5D, NULL);

  for (j = 0; j < chimes_table_bins.N_cool_4d_Temperatures; j++) {
    for (k = 0; k < chimes_table_bins.N_cool_4d_HIDensities; k++) {
      for (l = 0; l < chimes_table_bins.N_cool_4d_ElectronDensities; l++) {
        for (m = 0; m < chimes_table_bins.N_cool_4d_HIIDensities; m++) {
          offset5D[0] = 0;
          offset5D[1] = j;
          offset5D[2] = k;
          offset5D[3] = l;
          offset5D[4] = l;
          count5D[0] = chimes_master_cooling.N_coolants_4d;
          count5D[1] = 1;
          count5D[2] = 1;
          count5D[3] = 1;
          count5D[4] = 1;
          H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset5D, NULL,
                              count5D, NULL);

          H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                  H5P_DEFAULT, array_buffer_float);

          for (i = 0; i < chimes_master_cooling.N_coolants_4d; i++)
            chimes_master_cooling.rates_4d[chimes_flatten_index_5d(
                i, j, k, l, m, chimes_table_bins.N_cool_4d_Temperatures,
                chimes_table_bins.N_cool_4d_HIDensities,
                chimes_table_bins.N_cool_4d_ElectronDensities,
                chimes_table_bins.N_cool_4d_HIIDensities)] =
                chimes_max((ChimesFloat)array_buffer_float[i],
                           log_chimes_flt_min);
        }
      }
    }
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/rates_hiT_4d");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = chimes_master_cooling.N_coolants_4d;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_cool_hiT_4d_Temperatures; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = chimes_master_cooling.N_coolants_4d;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < chimes_master_cooling.N_coolants_4d; i++)
      chimes_master_cooling.rates_hiT_4d[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_cool_hiT_4d_Temperatures)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_Temperatures * sizeof(float));

  dataset = H5Dopen1(file_id, "cooling/photoelectric_heating");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = chimes_table_bins.N_Temperatures;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Psi; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = chimes_table_bins.N_Temperatures;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
      chimes_master_cooling.photoelectric_heating[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Psi)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "cooling/grain_recombination");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
  dims[0] = chimes_table_bins.N_Temperatures;
  dims[1] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for (j = 0; j < chimes_table_bins.N_Psi; j++) {
    offset2D[0] = 0;
    offset2D[1] = j;
    count2D[0] = chimes_table_bins.N_Temperatures;
    count2D[1] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                        NULL);

    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
            array_buffer_float);

    for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
      chimes_master_cooling.grain_recombination[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Psi)] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  free(array_buffer_float);

  array_buffer_float =
      (float *)malloc(chimes_table_bins.N_Temperatures * sizeof(float));

  dataset = H5Dopen1(file_id, "cooling/gas_grain_transfer");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++)
    chimes_master_cooling.gas_grain_transfer[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  free(array_buffer_float);

  array_buffer_float = (float *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(float));

  dataset = H5Dopen1(file_id, "cooling/H2_cool_lowDens_H2");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_master_cooling.H2_cool_lowDens_H2[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "cooling/H2_cool_lowDens_HI");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_master_cooling.H2_cool_lowDens_HI[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "cooling/H2_cool_lowDens_HII");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_master_cooling.H2_cool_lowDens_HII[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "cooling/H2_cool_lowDens_HeI");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_master_cooling.H2_cool_lowDens_HeI[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "cooling/H2_cool_lowDens_elec");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_master_cooling.H2_cool_lowDens_elec[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  dataset = H5Dopen1(file_id, "cooling/H2_cool_LTE");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
    chimes_master_cooling.H2_cool_LTE[i] =
        chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

  if ((myGlobalVars->element_included[0] == 1) &&
      (myGlobalVars->element_included[2] == 1)) {
    dataset = H5Dopen1(file_id, "cooling/CO_cool_rot_L0");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
      chimes_table_cooling.CO_cool_rot_L0[i] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

    dataset = H5Dopen1(file_id, "cooling/CO_cool_vib_L0");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
      chimes_table_cooling.CO_cool_vib_L0[i] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

    dataset = H5Dopen1(file_id, "cooling/CO_cool_rot_Llte");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_mol_cool_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_CO_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_mol_cool_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
        chimes_table_cooling.CO_cool_rot_Llte[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_CO_cool_rot_ColumnDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/CO_cool_rot_nhalf");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_mol_cool_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_CO_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_mol_cool_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
        chimes_table_cooling.CO_cool_rot_nhalf[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_CO_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/CO_cool_rot_a");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_mol_cool_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_CO_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_mol_cool_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
        chimes_table_cooling.CO_cool_rot_a[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_CO_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/CO_cool_vib_Llte");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_mol_cool_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_CO_cool_vib_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_mol_cool_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
        chimes_table_cooling.CO_cool_vib_Llte[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_CO_cool_vib_ColumnDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);
  }

  free(array_buffer_float);

  if (myGlobalVars->element_included[2] == 1) {
    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_H2O_cool_hiT_Temperatures * sizeof(float));

    dataset = H5Dopen1(file_id, "cooling/H2O_cool_rot_hiT_L0");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_hiT_Temperatures; i++)
      chimes_table_cooling.H2O_cool_rot_hiT_L0[i] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

    dataset = H5Dopen1(file_id, "cooling/H2O_cool_rot_hiT_Llte");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_hiT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_hiT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_hiT_Temperatures; i++)
        chimes_table_cooling.H2O_cool_rot_hiT_Llte[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2O_cool_rot_hiT_nhalf");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_hiT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_hiT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_hiT_Temperatures; i++)
        chimes_table_cooling.H2O_cool_rot_hiT_nhalf[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2O_cool_rot_hiT_a");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_hiT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_hiT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_hiT_Temperatures; i++)
        chimes_table_cooling.H2O_cool_rot_hiT_a[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_H2O_cool_lowT_Temperatures * sizeof(float));

    dataset = H5Dopen1(file_id, "cooling/H2Oortho_cool_rot_lowT_L0");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
      chimes_table_cooling.H2Oortho_cool_rot_lowT_L0[i] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

    dataset = H5Dopen1(file_id, "cooling/H2Oortho_cool_rot_lowT_Llte");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
        chimes_table_cooling
            .H2Oortho_cool_rot_lowT_Llte[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2Oortho_cool_rot_lowT_nhalf");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
        chimes_table_cooling
            .H2Oortho_cool_rot_lowT_nhalf[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2Oortho_cool_rot_lowT_a");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
        chimes_table_cooling.H2Oortho_cool_rot_lowT_a[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2Opara_cool_rot_lowT_L0");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
      chimes_table_cooling.H2Opara_cool_rot_lowT_L0[i] =
          chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);

    dataset = H5Dopen1(file_id, "cooling/H2Opara_cool_rot_lowT_Llte");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
        chimes_table_cooling.H2Opara_cool_rot_lowT_Llte[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2Opara_cool_rot_lowT_nhalf");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
        chimes_table_cooling
            .H2Opara_cool_rot_lowT_nhalf[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "cooling/H2Opara_cool_rot_lowT_a");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_rot_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_H2O_cool_lowT_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_H2O_cool_lowT_Temperatures; i++)
        chimes_table_cooling.H2Opara_cool_rot_lowT_a[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_rot_ColumnDensities)] =
            (ChimesFloat)array_buffer_float[i];
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    free(array_buffer_float);

    array_buffer_float = (float *)malloc(
        chimes_table_bins.N_mol_cool_Temperatures * sizeof(float));

    dataset = H5Dopen1(file_id, "cooling/H2O_cool_vib_L0");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
      chimes_table_cooling.H2O_cool_vib_L0[i] =
          (ChimesFloat)array_buffer_float[i];

    dataset = H5Dopen1(file_id, "cooling/H2O_cool_vib_Llte");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
    dims[0] = chimes_table_bins.N_mol_cool_Temperatures;
    dims[1] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims, NULL);

    for (j = 0; j < chimes_table_bins.N_H2O_cool_vib_ColumnDensities; j++) {
      offset2D[0] = 0;
      offset2D[1] = j;
      count2D[0] = chimes_table_bins.N_mol_cool_Temperatures;
      count2D[1] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++)
        chimes_table_cooling.H2O_cool_vib_Llte[chimes_flatten_index_2d(
            i, j, chimes_table_bins.N_H2O_cool_vib_ColumnDensities)] =
            chimes_max((ChimesFloat)array_buffer_float[i], log_chimes_flt_min);
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    free(array_buffer_float);
  }

  // finished with chimes_main_data.hdf5
  H5Fclose(file_id);

  /* We now need to read the cross-sections
   * tables for each spectrum into the
   * various chimes_master_photoion_XYZ
   * structures. */
  if (myGlobalVars->N_spectra > 0)
    read_cross_sections_tables(
        &chimes_table_bins, &chimes_master_photoion_fuv,
        &chimes_master_photoion_euv, &chimes_master_photoion_auger_fuv,
        &chimes_master_photoion_auger_euv, &chimes_table_spectra, myGlobalVars);

  /********************************************
   ** Now that we have read the data into    **
   ** the master tables, we need to          **
   ** determine which reactions are included **
   ** in the network, and copy them over to  **
   ** to the global tables.                  **
   ********************************************/

  /***************************
   ** T_dependent reactions **
   ***************************/
  chimes_table_T_dependent.N_reactions[0] = 0;
  chimes_table_T_dependent.N_reactions[1] = 0;
  for (i = 0; i < chimes_master_T_dependent.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_T_dependent.element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_T_dependent.N_reactions[1] += 1;
      if (chimes_master_T_dependent.molecular_flag[i] == 0)
        chimes_table_T_dependent.N_reactions[0] += 1;
    }
  }

  N_reactions_all = chimes_table_T_dependent.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_T_dependent.reactants =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_table_T_dependent.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_table_T_dependent.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_table_T_dependent.rates = (ChimesFloat *)malloc(
      N_reactions_all * chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));
  chimes_table_T_dependent.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));

  // Copy table data
  incl_index = 0;
  for (i = 0; i < chimes_master_T_dependent.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_T_dependent.element_incl, i,
            myGlobalVars->element_included)) {
      for (j = 0; j < 3; j++) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. However, an index of -1 should be left
         * as -1. */
        if (chimes_master_T_dependent
                .reactants[chimes_flatten_index_2d(i, j, 3)] < 0)
          chimes_table_T_dependent
              .reactants[chimes_flatten_index_2d(incl_index, j, 3)] =
              chimes_master_T_dependent
                  .reactants[chimes_flatten_index_2d(i, j, 3)];
        else
          chimes_table_T_dependent
              .reactants[chimes_flatten_index_2d(incl_index, j, 3)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_T_dependent.reactants
                                       [chimes_flatten_index_2d(i, j, 3)]];

        if (chimes_master_T_dependent
                .products[chimes_flatten_index_2d(i, j, 3)] < 0)
          chimes_table_T_dependent
              .products[chimes_flatten_index_2d(incl_index, j, 3)] =
              chimes_master_T_dependent
                  .products[chimes_flatten_index_2d(i, j, 3)];
        else
          chimes_table_T_dependent
              .products[chimes_flatten_index_2d(incl_index, j, 3)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_T_dependent.products
                                       [chimes_flatten_index_2d(i, j, 3)]];
      }

      for (j = 0; j < 9; j++)
        chimes_table_T_dependent
            .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
            chimes_master_T_dependent
                .element_incl[chimes_flatten_index_2d(i, j, 9)];

      for (j = 0; j < chimes_table_bins.N_Temperatures; j++)
        chimes_table_T_dependent.rates[chimes_flatten_index_2d(
            incl_index, j, chimes_table_bins.N_Temperatures)] =
            chimes_master_T_dependent.rates[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_Temperatures)];

      chimes_table_T_dependent.molecular_flag[incl_index] =
          chimes_master_T_dependent.molecular_flag[i];

      if (i ==
          chimes_master_T_dependent.H2_collis_dissoc_heating_reaction_index)
        chimes_table_T_dependent.H2_collis_dissoc_heating_reaction_index =
            incl_index;

      if (i == chimes_master_T_dependent.H2_form_heating_reaction_index)
        chimes_table_T_dependent.H2_form_heating_reaction_index = incl_index;

      incl_index++;
    }
  }

  /************************
   ** constant reactions **
   ************************/
  chimes_table_constant.N_reactions[0] = 0;
  chimes_table_constant.N_reactions[1] = 0;
  for (i = 0; i < chimes_master_constant.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(chimes_master_constant.element_incl,
                                            i,
                                            myGlobalVars->element_included)) {
      chimes_table_constant.N_reactions[1] += 1;
      if (chimes_master_constant.molecular_flag[i] == 0)
        chimes_table_constant.N_reactions[0] += 1;
    }
  }

  N_reactions_all = chimes_table_constant.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_constant.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_table_constant.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_table_constant.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_table_constant.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_constant.rates =
      (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));

  // Copy table data
  incl_index = 0;
  for (i = 0; i < chimes_master_constant.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(chimes_master_constant.element_incl,
                                            i,
                                            myGlobalVars->element_included)) {
      /* The reactants and products arrays need to be translated into
       * the indices in the reduced network, as given by
       * myGlobalVars->speciesIndices. However, an index of -1 should be left as
       * -1. */
      for (j = 0; j < 2; j++) {
        if (chimes_master_constant.reactants[chimes_flatten_index_2d(i, j, 2)] <
            0)
          chimes_table_constant
              .reactants[chimes_flatten_index_2d(incl_index, j, 2)] =
              chimes_master_constant
                  .reactants[chimes_flatten_index_2d(i, j, 2)];
        else
          chimes_table_constant
              .reactants[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_constant.reactants
                                       [chimes_flatten_index_2d(i, j, 2)]];
      }

      for (j = 0; j < 3; j++) {
        if (chimes_master_constant
                .products[chimes_flatten_index_2d(incl_index, j, 3)] < 0)
          chimes_table_constant
              .products[chimes_flatten_index_2d(incl_index, j, 3)] =
              chimes_master_constant.products[chimes_flatten_index_2d(i, j, 3)];
        else
          chimes_table_constant
              .products[chimes_flatten_index_2d(incl_index, j, 3)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_constant.products
                                       [chimes_flatten_index_2d(i, j, 3)]];
      }

      for (j = 0; j < 9; j++)
        chimes_table_constant
            .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
            chimes_master_constant
                .element_incl[chimes_flatten_index_2d(i, j, 9)];

      chimes_table_constant.rates[incl_index] = chimes_master_constant.rates[i];
      chimes_table_constant.molecular_flag[incl_index] =
          chimes_master_constant.molecular_flag[i];

      if (i == chimes_master_constant.H2_form_heating_reaction_index)
        chimes_table_constant.H2_form_heating_reaction_index = incl_index;

      incl_index++;
    }
  }

  /********************************
   ** recombination_AB reactions **
   ********************************/
  chimes_table_recombination_AB.N_reactions[0] = 0;
  chimes_table_recombination_AB.N_reactions[1] = 0;
  for (i = 0; i < chimes_master_recombination_AB.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_recombination_AB.element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_recombination_AB.N_reactions[1] += 1;
      if (chimes_master_recombination_AB.molecular_flag[i] == 0)
        chimes_table_recombination_AB.N_reactions[0] += 1;
    }
  }

  N_reactions_all = chimes_table_recombination_AB.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_recombination_AB.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_table_recombination_AB.products =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_recombination_AB.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_table_recombination_AB.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_recombination_AB.rates = (ChimesFloat *)malloc(
      N_reactions_all * 2 * chimes_table_bins.N_Temperatures *
      sizeof(ChimesFloat));

  // Copy table data
  incl_index = 0;
  for (i = 0; i < chimes_master_recombination_AB.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_recombination_AB.element_incl, i,
            myGlobalVars->element_included)) {
      for (j = 0; j < 2; j++) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. However, an index of -1 should be left
         * as -1. */
        if (chimes_master_recombination_AB
                .reactants[chimes_flatten_index_2d(i, j, 2)] < 0)
          chimes_table_recombination_AB
              .reactants[chimes_flatten_index_2d(incl_index, j, 2)] =
              chimes_master_recombination_AB
                  .reactants[chimes_flatten_index_2d(i, j, 2)];
        else
          chimes_table_recombination_AB
              .reactants[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_recombination_AB.reactants
                                       [chimes_flatten_index_2d(i, j, 2)]];
      }

      if (chimes_master_recombination_AB.products[i] < 0)
        chimes_table_recombination_AB.products[incl_index] =
            chimes_master_recombination_AB.products[i];
      else
        chimes_table_recombination_AB.products[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_recombination_AB.products[i]];

      for (j = 0; j < 9; j++)
        chimes_table_recombination_AB
            .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
            chimes_master_recombination_AB
                .element_incl[chimes_flatten_index_2d(i, j, 9)];

      for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
        chimes_table_recombination_AB.rates[chimes_flatten_index_3d(
            incl_index, 0, j, 2, chimes_table_bins.N_Temperatures)] =
            chimes_master_recombination_AB.rates[chimes_flatten_index_3d(
                i, 0, j, 2, chimes_table_bins.N_Temperatures)];
        chimes_table_recombination_AB.rates[chimes_flatten_index_3d(
            incl_index, 1, j, 2, chimes_table_bins.N_Temperatures)] =
            chimes_master_recombination_AB.rates[chimes_flatten_index_3d(
                i, 1, j, 2, chimes_table_bins.N_Temperatures)];
      }

      chimes_table_recombination_AB.molecular_flag[incl_index] =
          chimes_master_recombination_AB.molecular_flag[i];
      incl_index++;
    }
  }

  /********************************
   ** grain_recombination reactions **
   ********************************/
  chimes_table_grain_recombination.N_reactions[0] = 0;
  chimes_table_grain_recombination.N_reactions[1] = 0;
  for (i = 0; i < chimes_master_grain_recombination.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_grain_recombination.element_incl, i,
            myGlobalVars->element_included)) {
      // All reactions are non-molecular
      chimes_table_grain_recombination.N_reactions[0] += 1;
      chimes_table_grain_recombination.N_reactions[1] += 1;
    }
  }

  N_reactions_all = chimes_table_grain_recombination.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_grain_recombination.reactants =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_table_grain_recombination.products =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_grain_recombination.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_table_grain_recombination.rates =
      (ChimesFloat *)malloc(N_reactions_all * chimes_table_bins.N_Temperatures *
                            chimes_table_bins.N_Psi * sizeof(ChimesFloat));

  // Copy table data
  incl_index = 0;
  for (i = 0; i < chimes_master_grain_recombination.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_grain_recombination.element_incl, i,
            myGlobalVars->element_included)) {
      for (j = 0; j < 2; j++)
        chimes_table_grain_recombination
            .reactants[chimes_flatten_index_2d(incl_index, j, 2)] =
            myGlobalVars->speciesIndices[chimes_master_grain_recombination
                                             .reactants[chimes_flatten_index_2d(
                                                 i, j, 2)]];

      chimes_table_grain_recombination.products[incl_index] =
          myGlobalVars
              ->speciesIndices[chimes_master_grain_recombination.products[i]];

      for (j = 0; j < 9; j++)
        chimes_table_grain_recombination
            .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
            chimes_master_grain_recombination
                .element_incl[chimes_flatten_index_2d(i, j, 9)];

      for (j = 0; j < chimes_table_bins.N_Temperatures; j++) {
        for (k = 0; k < chimes_table_bins.N_Psi; k++)
          chimes_table_grain_recombination.rates[chimes_flatten_index_3d(
              incl_index, j, k, chimes_table_bins.N_Temperatures,
              chimes_table_bins.N_Psi)] =
              chimes_master_grain_recombination.rates[chimes_flatten_index_3d(
                  i, j, k, chimes_table_bins.N_Temperatures,
                  chimes_table_bins.N_Psi)];
      }

      incl_index++;
    }
  }

  /**************************
   ** cosmic_ray reactions **
   **************************/
  chimes_table_cosmic_ray.N_reactions[0] = 0;
  chimes_table_cosmic_ray.N_reactions[1] = 0;
  for (i = 0; i < chimes_master_cosmic_ray.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_cosmic_ray.element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_cosmic_ray.N_reactions[1] += 1;
      if (chimes_master_cosmic_ray.molecular_flag[i] == 0)
        chimes_table_cosmic_ray.N_reactions[0] += 1;
    }
  }

  N_reactions_all = chimes_table_cosmic_ray.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_cosmic_ray.reactants =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_cosmic_ray.products =
      (int *)malloc(N_reactions_all * 3 * sizeof(int));
  chimes_table_cosmic_ray.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_table_cosmic_ray.molecular_flag =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_cosmic_ray.rates =
      (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
  chimes_table_cosmic_ray.secondary_base_reaction =
      (int *)malloc(2 * sizeof(int));
  chimes_table_cosmic_ray.secondary_ratio = (ChimesFloat *)malloc(
      2 * chimes_table_bins.N_secondary_cosmic_ray_xHII * sizeof(ChimesFloat));

  // Copy table data
  incl_index = 0;
  for (i = 0; i < chimes_master_cosmic_ray.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_cosmic_ray.element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_cosmic_ray.reactants[incl_index] =
          myGlobalVars->speciesIndices[chimes_master_cosmic_ray.reactants[i]];

      for (j = 0; j < 3; j++) {
        if (chimes_master_cosmic_ray
                .products[chimes_flatten_index_2d(i, j, 3)] < 0)
          chimes_table_cosmic_ray
              .products[chimes_flatten_index_2d(incl_index, j, 3)] =
              chimes_master_cosmic_ray
                  .products[chimes_flatten_index_2d(i, j, 3)];
        else
          chimes_table_cosmic_ray
              .products[chimes_flatten_index_2d(incl_index, j, 3)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_cosmic_ray.products
                                       [chimes_flatten_index_2d(i, j, 3)]];
      }

      for (j = 0; j < 9; j++)
        chimes_table_cosmic_ray
            .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
            chimes_master_cosmic_ray
                .element_incl[chimes_flatten_index_2d(i, j, 9)];

      chimes_table_cosmic_ray.rates[incl_index] =
          chimes_master_cosmic_ray.rates[i];
      chimes_table_cosmic_ray.molecular_flag[incl_index] =
          chimes_master_cosmic_ray.molecular_flag[i];

      for (j = 0; j < 2; j++) {
        if (i == chimes_master_cosmic_ray.secondary_base_reaction[j]) {
          chimes_table_cosmic_ray.secondary_base_reaction[j] = incl_index;

          for (k = 0; k < chimes_table_bins.N_secondary_cosmic_ray_xHII; k++)
            chimes_table_cosmic_ray.secondary_ratio[chimes_flatten_index_2d(
                j, k, chimes_table_bins.N_secondary_cosmic_ray_xHII)] =
                chimes_master_cosmic_ray
                    .secondary_ratio[chimes_flatten_index_2d(
                        j, k, chimes_table_bins.N_secondary_cosmic_ray_xHII)];
        }
      }

      incl_index++;
    }
  }

  /*****************************
   ** CO_cosmic_ray reactions **
   *****************************/
  chimes_table_CO_cosmic_ray.N_reactions[0] = 0;
  chimes_table_CO_cosmic_ray.N_reactions[1] = 0;
  for (i = 0; i < chimes_master_CO_cosmic_ray.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_CO_cosmic_ray.element_incl, i,
            myGlobalVars->element_included)) {
      // Only molecular reactions in this group.
      chimes_table_CO_cosmic_ray.N_reactions[1] += 1;
    }
  }

  N_reactions_all = chimes_table_CO_cosmic_ray.N_reactions[1];

  // Allocate memory in table structure
  chimes_table_CO_cosmic_ray.reactants =
      (int *)malloc(N_reactions_all * sizeof(int));
  chimes_table_CO_cosmic_ray.products =
      (int *)malloc(N_reactions_all * 2 * sizeof(int));
  chimes_table_CO_cosmic_ray.element_incl =
      (int *)malloc(N_reactions_all * 9 * sizeof(int));
  chimes_table_CO_cosmic_ray.rates = (ChimesFloat *)malloc(
      N_reactions_all * chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));

  // Copy table data
  incl_index = 0;
  for (i = 0; i < chimes_master_CO_cosmic_ray.N_reactions[1]; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_master_CO_cosmic_ray.element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_CO_cosmic_ray.reactants[incl_index] =
          myGlobalVars
              ->speciesIndices[chimes_master_CO_cosmic_ray.reactants[i]];

      for (j = 0; j < 2; j++)
        chimes_table_CO_cosmic_ray
            .products[chimes_flatten_index_2d(incl_index, j, 2)] =
            myGlobalVars->speciesIndices[chimes_master_CO_cosmic_ray
                                             .products[chimes_flatten_index_2d(
                                                 i, j, 2)]];

      for (j = 0; j < 9; j++)
        chimes_table_CO_cosmic_ray
            .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
            chimes_master_CO_cosmic_ray
                .element_incl[chimes_flatten_index_2d(i, j, 9)];

      for (j = 0; j < chimes_table_bins.N_Temperatures; j++)
        chimes_table_CO_cosmic_ray.rates[chimes_flatten_index_2d(
            incl_index, j, chimes_table_bins.N_Temperatures)] =
            chimes_master_CO_cosmic_ray.rates[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_Temperatures)];

      incl_index++;
    }
  }

  if (myGlobalVars->N_spectra > 0) {
    /***************************
     ** photoion_fuv reactions **
     ***************************/
    chimes_table_photoion_fuv.N_reactions[0] = 0;
    chimes_table_photoion_fuv.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_photoion_fuv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_fuv.element_incl, i,
              myGlobalVars->element_included)) {
        // No molecular reactions in this group
        chimes_table_photoion_fuv.N_reactions[0] += 1;
        chimes_table_photoion_fuv.N_reactions[1] += 1;
      }
    }

    N_reactions_all = chimes_table_photoion_fuv.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_photoion_fuv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_fuv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_photoion_fuv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_photoion_fuv.sigmaPhot = (ChimesFloat *)malloc(
        myGlobalVars->N_spectra * N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photoion_fuv.epsilonPhot = (ChimesFloat *)malloc(
        myGlobalVars->N_spectra * N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photoion_fuv.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_photoion_fuv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_fuv.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_photoion_fuv.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_photoion_fuv.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_photoion_fuv
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_photoion_fuv.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_photoion_fuv
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_photoion_fuv
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        for (j = 0; j < myGlobalVars->N_spectra; j++) {
          chimes_table_photoion_fuv.sigmaPhot[chimes_flatten_index_2d(
              j, incl_index, N_reactions_all)] =
              chimes_master_photoion_fuv.sigmaPhot[chimes_flatten_index_2d(
                  j, i, chimes_master_photoion_fuv.N_reactions[1])];
          chimes_table_photoion_fuv.epsilonPhot[chimes_flatten_index_2d(
              j, incl_index, N_reactions_all)] =
              chimes_master_photoion_fuv.epsilonPhot[chimes_flatten_index_2d(
                  j, i, chimes_master_photoion_fuv.N_reactions[1])];
        }

        chimes_table_photoion_fuv.gamma[incl_index] =
            chimes_master_photoion_fuv.gamma[i];
        incl_index++;
      }
    }

    /***************************
     ** photoion_euv reactions **
     ***************************/
    chimes_table_photoion_euv.N_reactions[0] = 0;
    chimes_table_photoion_euv.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_photoion_euv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_euv.element_incl, i,
              myGlobalVars->element_included)) {
        chimes_table_photoion_euv.N_reactions[1] += 1;
        if (chimes_master_photoion_euv.molecular_flag[i] == 0)
          chimes_table_photoion_euv.N_reactions[0] += 1;
      }
    }

    N_reactions_all = chimes_table_photoion_euv.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_photoion_euv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_euv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_photoion_euv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_photoion_euv.E_thresh =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photoion_euv.molecular_flag =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_euv.sigmaPhot = (ChimesFloat *)malloc(
        myGlobalVars->N_spectra * N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photoion_euv.shieldFactor_1D = (ChimesFloat *)malloc(
        N_reactions_all * myGlobalVars->N_spectra * 3 *
        chimes_table_bins.N_Column_densities * sizeof(ChimesFloat));
    chimes_table_photoion_euv.shieldFactor_2D = (ChimesFloat *)malloc(
        N_reactions_all * myGlobalVars->N_spectra * 6 *
        chimes_table_bins.N_Column_densities *
        chimes_table_bins.N_Column_densities * sizeof(ChimesFloat));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_photoion_euv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_euv.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_photoion_euv.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_photoion_euv.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_photoion_euv
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_photoion_euv.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_photoion_euv
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_photoion_euv
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        for (j = 0; j < myGlobalVars->N_spectra; j++) {
          chimes_table_photoion_euv.sigmaPhot[chimes_flatten_index_2d(
              j, incl_index, N_reactions_all)] =
              chimes_master_photoion_euv.sigmaPhot[chimes_flatten_index_2d(
                  j, i, chimes_master_photoion_euv.N_reactions[1])];

          for (k = 0; k < 3; k++) {
            for (l = 0; l < chimes_table_bins.N_Column_densities; l++)
              chimes_table_photoion_euv.shieldFactor_1D[chimes_flatten_index_4d(
                  incl_index, j, k, l, myGlobalVars->N_spectra, 3,
                  chimes_table_bins.N_Column_densities)] =
                  chimes_master_photoion_euv
                      .shieldFactor_1D[chimes_flatten_index_4d(
                          i, j, k, l, myGlobalVars->N_spectra, 3,
                          chimes_table_bins.N_Column_densities)];
          }

          for (k = 0; k < 6; k++) {
            for (l = 0; l < chimes_table_bins.N_Column_densities; l++) {
              for (m = 0; m < chimes_table_bins.N_Column_densities; m++)
                chimes_table_photoion_euv
                    .shieldFactor_2D[chimes_flatten_index_5d(
                        incl_index, j, k, l, m, myGlobalVars->N_spectra, 6,
                        chimes_table_bins.N_Column_densities,
                        chimes_table_bins.N_Column_densities)] =
                    chimes_master_photoion_euv
                        .shieldFactor_2D[chimes_flatten_index_5d(
                            i, j, k, l, m, myGlobalVars->N_spectra, 6,
                            chimes_table_bins.N_Column_densities,
                            chimes_table_bins.N_Column_densities)];
            }
          }
        }

        chimes_table_photoion_euv.molecular_flag[incl_index] =
            chimes_master_photoion_euv.molecular_flag[i];
        chimes_table_photoion_euv.E_thresh[incl_index] =
            chimes_master_photoion_euv.E_thresh[i];
        incl_index++;
      }
    }

    /**********************************
     ** photoion_auger_fuv reactions **
     **********************************/
    chimes_table_photoion_auger_fuv.N_reactions[0] = 0;
    chimes_table_photoion_auger_fuv.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_photoion_auger_fuv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_auger_fuv.element_incl, i,
              myGlobalVars->element_included)) {
        // No molecular reactions in this group
        chimes_table_photoion_auger_fuv.N_reactions[0] += 1;
        chimes_table_photoion_auger_fuv.N_reactions[1] += 1;
      }
    }

    N_reactions_all = chimes_table_photoion_auger_fuv.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_photoion_auger_fuv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_auger_fuv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_photoion_auger_fuv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_photoion_auger_fuv.sigmaPhot = (ChimesFloat *)malloc(
        myGlobalVars->N_spectra * N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photoion_auger_fuv.base_reaction =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_auger_fuv.number_of_electrons =
        (int *)malloc(N_reactions_all * sizeof(int));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_photoion_auger_fuv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_auger_fuv.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_photoion_auger_fuv.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_photoion_auger_fuv.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_photoion_auger_fuv
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_photoion_auger_fuv.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_photoion_auger_fuv
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_photoion_auger_fuv
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        for (j = 0; j < myGlobalVars->N_spectra; j++)
          chimes_table_photoion_auger_fuv.sigmaPhot[chimes_flatten_index_2d(
              j, incl_index, N_reactions_all)] =
              chimes_master_photoion_auger_fuv
                  .sigmaPhot[chimes_flatten_index_2d(
                      j, i, chimes_master_photoion_auger_fuv.N_reactions[1])];

        base_incl_index = 0;
        found_base = 0;
        for (base_index = 0;
             base_index < chimes_master_photoion_fuv.N_reactions[1];
             base_index++) {
          if (compare_element_incl_arrays_flatten(
                  chimes_master_photoion_fuv.element_incl, base_index,
                  myGlobalVars->element_included)) {
            if (base_index ==
                chimes_master_photoion_auger_fuv.base_reaction[i]) {
              chimes_table_photoion_auger_fuv.base_reaction[incl_index] =
                  base_incl_index;
              found_base = 1;
              break;
            }
            base_incl_index++;
          }
        }

        if (found_base == 0) {
          printf(
              "CHIMES ERROR: In photoion_auger_fuv group, base reaction not "
              "found in reduced network.\n");
          exit(EXIT_FAILURE);
        }

        /* Auger ionisations involve a single photon releasing
         * multiple electrons. Here we record how many electrons
         * are released in the reaction, which is given by
         * the difference between the final and initial
         * ionisation state. */
        chimes_table_photoion_auger_fuv.number_of_electrons[incl_index] =
            chimes_master_photoion_auger_fuv
                .products[chimes_flatten_index_2d(i, 0, 2)] -
            chimes_master_photoion_auger_fuv.reactants[i];

        incl_index++;
      }
    }

    /**********************************
     ** photoion_auger_euv reactions **
     **********************************/
    chimes_table_photoion_auger_euv.N_reactions[0] = 0;
    chimes_table_photoion_auger_euv.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_photoion_auger_euv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_auger_euv.element_incl, i,
              myGlobalVars->element_included)) {
        // No molecular reactions in this group
        chimes_table_photoion_auger_euv.N_reactions[0] += 1;
        chimes_table_photoion_auger_euv.N_reactions[1] += 1;
      }
    }

    N_reactions_all = chimes_table_photoion_auger_euv.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_photoion_auger_euv.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_auger_euv.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_photoion_auger_euv.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_photoion_auger_euv.sigmaPhot = (ChimesFloat *)malloc(
        myGlobalVars->N_spectra * N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photoion_auger_euv.base_reaction =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photoion_auger_euv.number_of_electrons =
        (int *)malloc(N_reactions_all * sizeof(int));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_photoion_auger_euv.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photoion_auger_euv.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_photoion_auger_euv.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_photoion_auger_euv.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_photoion_auger_euv
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_photoion_auger_euv.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_photoion_auger_euv
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_photoion_auger_euv
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        for (j = 0; j < myGlobalVars->N_spectra; j++)
          chimes_table_photoion_auger_euv.sigmaPhot[chimes_flatten_index_2d(
              j, incl_index, N_reactions_all)] =
              chimes_master_photoion_auger_euv
                  .sigmaPhot[chimes_flatten_index_2d(
                      j, i, chimes_master_photoion_auger_euv.N_reactions[1])];

        base_incl_index = 0;
        found_base = 0;
        for (base_index = 0;
             base_index < chimes_master_photoion_euv.N_reactions[1];
             base_index++) {
          if (compare_element_incl_arrays_flatten(
                  chimes_master_photoion_euv.element_incl, base_index,
                  myGlobalVars->element_included)) {
            if (base_index ==
                chimes_master_photoion_auger_euv.base_reaction[i]) {
              chimes_table_photoion_auger_euv.base_reaction[incl_index] =
                  base_incl_index;
              found_base = 1;
              break;
            }
            base_incl_index++;
          }
        }

        if (found_base == 0) {
          printf(
              "CHIMES ERROR: In photoion_auger_euv group, base reaction not "
              "found in reduced network.\n");
          exit(EXIT_FAILURE);
        }

        /* Auger ionisations involve a single photon releasing
         * multiple electrons. Here we record how many electrons
         * are released in the reaction, which is given by
         * the difference between the final and initial
         * ionisation state. */
        chimes_table_photoion_auger_euv.number_of_electrons[incl_index] =
            chimes_master_photoion_auger_euv
                .products[chimes_flatten_index_2d(i, 0, 2)] -
            chimes_master_photoion_auger_euv.reactants[i];

        incl_index++;
      }
    }

    /**********************************
     ** photodissoc_group1 reactions **
     **********************************/
    chimes_table_photodissoc_group1.N_reactions[0] = 0;
    chimes_table_photodissoc_group1.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_photodissoc_group1.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photodissoc_group1.element_incl, i,
              myGlobalVars->element_included)) {
        chimes_table_photodissoc_group1.N_reactions[1] += 1;
        if (chimes_master_photodissoc_group1.molecular_flag[i] == 0)
          chimes_table_photodissoc_group1.N_reactions[0] += 1;
      }
    }

    N_reactions_all = chimes_table_photodissoc_group1.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_photodissoc_group1.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photodissoc_group1.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_photodissoc_group1.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_photodissoc_group1.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photodissoc_group1.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photodissoc_group1.molecular_flag =
        (int *)malloc(N_reactions_all * sizeof(int));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_photodissoc_group1.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photodissoc_group1.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_photodissoc_group1.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_photodissoc_group1.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_photodissoc_group1
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_photodissoc_group1.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_photodissoc_group1
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_photodissoc_group1
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        chimes_table_photodissoc_group1.gamma[incl_index] =
            chimes_master_photodissoc_group1.gamma[i];
        chimes_table_photodissoc_group1.rates[incl_index] =
            chimes_master_photodissoc_group1.rates[i];
        chimes_table_photodissoc_group1.molecular_flag[incl_index] =
            chimes_master_photodissoc_group1.molecular_flag[i];
        incl_index++;
      }
    }

    /**********************************
     ** photodissoc_group2 reactions **
     **********************************/
    chimes_table_photodissoc_group2.N_reactions[0] = 0;
    chimes_table_photodissoc_group2.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_photodissoc_group2.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photodissoc_group2.element_incl, i,
              myGlobalVars->element_included)) {
        // Only molecular reactions
        chimes_table_photodissoc_group2.N_reactions[1] += 1;
      }
    }

    N_reactions_all = chimes_table_photodissoc_group2.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_photodissoc_group2.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_photodissoc_group2.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_photodissoc_group2.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_photodissoc_group2.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_photodissoc_group2.gamma_coeff =
        (ChimesFloat *)malloc(3 * sizeof(ChimesFloat));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_photodissoc_group2.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_photodissoc_group2.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_photodissoc_group2.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_photodissoc_group2.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_photodissoc_group2
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_photodissoc_group2.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_photodissoc_group2
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_photodissoc_group2
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        chimes_table_photodissoc_group2.rates[incl_index] =
            chimes_master_photodissoc_group2.rates[i];
        incl_index++;
      }
    }

    for (i = 0; i < 3; i++)
      chimes_table_photodissoc_group2.gamma_coeff[i] =
          chimes_master_photodissoc_group2.gamma_coeff[i];

    /**********************************
     ** CO_photodissoc reactions **
     **********************************/
    chimes_table_CO_photodissoc.N_reactions[0] = 0;
    chimes_table_CO_photodissoc.N_reactions[1] = 0;
    for (i = 0; i < chimes_master_CO_photodissoc.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_CO_photodissoc.element_incl, i,
              myGlobalVars->element_included)) {
        // Only molecular reactions.
        chimes_table_CO_photodissoc.N_reactions[1] += 1;
      }
    }

    N_reactions_all = chimes_table_CO_photodissoc.N_reactions[1];

    // Allocate memory in table structure
    chimes_table_CO_photodissoc.reactants =
        (int *)malloc(N_reactions_all * sizeof(int));
    chimes_table_CO_photodissoc.products =
        (int *)malloc(N_reactions_all * 2 * sizeof(int));
    chimes_table_CO_photodissoc.element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));
    chimes_table_CO_photodissoc.gamma =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_CO_photodissoc.rates =
        (ChimesFloat *)malloc(N_reactions_all * sizeof(ChimesFloat));
    chimes_table_CO_photodissoc.self_shielding = (ChimesFloat *)malloc(
        N_reactions_all * chimes_table_bins.N_COself_column_densities *
        chimes_table_bins.N_H2CO_column_densities * sizeof(ChimesFloat));

    // Copy table data
    incl_index = 0;
    for (i = 0; i < chimes_master_CO_photodissoc.N_reactions[1]; i++) {
      if (compare_element_incl_arrays_flatten(
              chimes_master_CO_photodissoc.element_incl, i,
              myGlobalVars->element_included)) {
        /* The reactants and products arrays need to be translated into
         * the indices in the reduced network, as given by
         * myGlobalVars->speciesIndices. */
        chimes_table_CO_photodissoc.reactants[incl_index] =
            myGlobalVars
                ->speciesIndices[chimes_master_CO_photodissoc.reactants[i]];

        for (j = 0; j < 2; j++)
          chimes_table_CO_photodissoc
              .products[chimes_flatten_index_2d(incl_index, j, 2)] =
              myGlobalVars
                  ->speciesIndices[chimes_master_CO_photodissoc.products
                                       [chimes_flatten_index_2d(i, j, 2)]];

        for (j = 0; j < 9; j++)
          chimes_table_CO_photodissoc
              .element_incl[chimes_flatten_index_2d(incl_index, j, 9)] =
              chimes_master_CO_photodissoc
                  .element_incl[chimes_flatten_index_2d(i, j, 9)];

        chimes_table_CO_photodissoc.gamma[incl_index] =
            chimes_master_CO_photodissoc.gamma[i];
        chimes_table_CO_photodissoc.rates[incl_index] =
            chimes_master_CO_photodissoc.rates[i];

        for (j = 0; j < chimes_table_bins.N_COself_column_densities; j++) {
          for (k = 0; k < chimes_table_bins.N_H2CO_column_densities; k++)
            chimes_table_CO_photodissoc.self_shielding[chimes_flatten_index_3d(
                incl_index, j, k, chimes_table_bins.N_COself_column_densities,
                chimes_table_bins.N_H2CO_column_densities)] =
                chimes_master_CO_photodissoc
                    .self_shielding[chimes_flatten_index_3d(
                        i, j, k, chimes_table_bins.N_COself_column_densities,
                        chimes_table_bins.N_H2CO_column_densities)];
        }

        incl_index++;
      }
    }
  }

  /*************
   ** cooling **
   *************/
  chimes_table_cooling.N_coolants = 0;
  chimes_table_cooling.N_coolants_2d = 0;
  chimes_table_cooling.N_coolants_4d = 0;

  /* Use the speciesIndices array to determine
   * which coolants are included in the network.
   * It will be -1 for excluded species. */
  for (i = 0; i < chimes_master_cooling.N_coolants; i++) {
    if (myGlobalVars->speciesIndices[chimes_master_cooling.coolants[i]] >= 0)
      chimes_table_cooling.N_coolants++;
  }

  for (i = 0; i < chimes_master_cooling.N_coolants_2d; i++) {
    if (myGlobalVars->speciesIndices[chimes_master_cooling.coolants_2d[i]] >= 0)
      chimes_table_cooling.N_coolants_2d++;
  }

  for (i = 0; i < chimes_master_cooling.N_coolants_4d; i++) {
    if (myGlobalVars->speciesIndices[chimes_master_cooling.coolants_4d[i]] >= 0)
      chimes_table_cooling.N_coolants_4d++;
  }

  // Allocate memory in table structure
  chimes_table_cooling.coolants =
      (int *)malloc(chimes_table_cooling.N_coolants * sizeof(int));
  chimes_table_cooling.coolants_2d =
      (int *)malloc(chimes_table_cooling.N_coolants_2d * sizeof(int));
  chimes_table_cooling.coolants_4d =
      (int *)malloc(chimes_table_cooling.N_coolants_4d * sizeof(int));

  chimes_table_cooling.rates = (ChimesFloat *)malloc(
      chimes_table_cooling.N_coolants * chimes_table_bins.N_Temperatures *
      sizeof(ChimesFloat));

  chimes_table_cooling.rates_2d = (ChimesFloat *)malloc(
      chimes_table_cooling.N_coolants_2d *
      chimes_table_bins.N_cool_2d_Temperatures *
      chimes_table_bins.N_cool_2d_ElectronDensities * sizeof(ChimesFloat));
  chimes_table_cooling.rates_hiT_2d = (ChimesFloat *)malloc(
      chimes_table_cooling.N_coolants_2d *
      chimes_table_bins.N_cool_hiT_2d_Temperatures * sizeof(ChimesFloat));

  chimes_table_cooling.rates_4d = (ChimesFloat *)malloc(
      chimes_table_cooling.N_coolants_4d *
      chimes_table_bins.N_cool_4d_Temperatures *
      chimes_table_bins.N_cool_4d_HIDensities *
      chimes_table_bins.N_cool_4d_ElectronDensities *
      chimes_table_bins.N_cool_4d_HIIDensities * sizeof(ChimesFloat));
  chimes_table_cooling.rates_hiT_4d = (ChimesFloat *)malloc(
      chimes_table_cooling.N_coolants_4d *
      chimes_table_bins.N_cool_hiT_4d_Temperatures * sizeof(ChimesFloat));

  chimes_table_cooling.photoelectric_heating =
      (ChimesFloat *)malloc(chimes_table_bins.N_Temperatures *
                            chimes_table_bins.N_Psi * sizeof(ChimesFloat));
  chimes_table_cooling.grain_recombination =
      (ChimesFloat *)malloc(chimes_table_bins.N_Temperatures *
                            chimes_table_bins.N_Psi * sizeof(ChimesFloat));
  chimes_table_cooling.gas_grain_transfer = (ChimesFloat *)malloc(
      chimes_table_bins.N_Temperatures * sizeof(ChimesFloat));

  chimes_table_cooling.H2_cool_lowDens_H2 = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_table_cooling.H2_cool_lowDens_HI = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_table_cooling.H2_cool_lowDens_HII = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_table_cooling.H2_cool_lowDens_HeI = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_table_cooling.H2_cool_lowDens_elec = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));
  chimes_table_cooling.H2_cool_LTE = (ChimesFloat *)malloc(
      chimes_table_bins.N_mol_cool_Temperatures * sizeof(ChimesFloat));

  // Copy over data from master tables
  incl_index = 0;
  for (i = 0; i < chimes_master_cooling.N_coolants; i++) {
    if (myGlobalVars->speciesIndices[chimes_master_cooling.coolants[i]] >= 0) {
      chimes_table_cooling.coolants[incl_index] =
          myGlobalVars->speciesIndices[chimes_master_cooling.coolants[i]];
      for (j = 0; j < chimes_table_bins.N_Temperatures; j++)
        chimes_table_cooling.rates[chimes_flatten_index_2d(
            incl_index, j, chimes_table_bins.N_Temperatures)] =
            chimes_master_cooling.rates[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_Temperatures)];

      incl_index++;
    }
  }

  incl_index = 0;
  for (i = 0; i < chimes_master_cooling.N_coolants_2d; i++) {
    if (myGlobalVars->speciesIndices[chimes_master_cooling.coolants_2d[i]] >=
        0) {
      chimes_table_cooling.coolants_2d[incl_index] =
          myGlobalVars->speciesIndices[chimes_master_cooling.coolants_2d[i]];

      for (j = 0; j < chimes_table_bins.N_cool_2d_Temperatures; j++) {
        for (k = 0; k < chimes_table_bins.N_cool_2d_ElectronDensities; k++)
          chimes_table_cooling.rates_2d[chimes_flatten_index_3d(
              incl_index, j, k, chimes_table_bins.N_cool_2d_Temperatures,
              chimes_table_bins.N_cool_2d_ElectronDensities)] =
              chimes_master_cooling.rates_2d[chimes_flatten_index_3d(
                  i, j, k, chimes_table_bins.N_cool_2d_Temperatures,
                  chimes_table_bins.N_cool_2d_ElectronDensities)];
      }

      for (j = 0; j < chimes_table_bins.N_cool_hiT_2d_Temperatures; j++)
        chimes_table_cooling.rates_hiT_2d[chimes_flatten_index_2d(
            incl_index, j, chimes_table_bins.N_cool_hiT_2d_Temperatures)] =
            chimes_master_cooling.rates_hiT_2d[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_cool_hiT_2d_Temperatures)];

      incl_index++;
    }
  }

  incl_index = 0;
  for (i = 0; i < chimes_master_cooling.N_coolants_4d; i++) {
    if (myGlobalVars->speciesIndices[chimes_master_cooling.coolants_4d[i]] >=
        0) {
      chimes_table_cooling.coolants_4d[incl_index] =
          myGlobalVars->speciesIndices[chimes_master_cooling.coolants_4d[i]];

      for (j = 0; j < chimes_table_bins.N_cool_4d_Temperatures; j++) {
        for (k = 0; k < chimes_table_bins.N_cool_4d_HIDensities; k++) {
          for (l = 0; l < chimes_table_bins.N_cool_4d_ElectronDensities; l++) {
            for (m = 0; m < chimes_table_bins.N_cool_4d_HIIDensities; m++)
              chimes_table_cooling.rates_4d[chimes_flatten_index_5d(
                  incl_index, j, k, l, m,
                  chimes_table_bins.N_cool_4d_Temperatures,
                  chimes_table_bins.N_cool_4d_HIDensities,
                  chimes_table_bins.N_cool_4d_ElectronDensities,
                  chimes_table_bins.N_cool_4d_HIIDensities)] =
                  chimes_master_cooling.rates_4d[chimes_flatten_index_5d(
                      i, j, k, l, m, chimes_table_bins.N_cool_4d_Temperatures,
                      chimes_table_bins.N_cool_4d_HIDensities,
                      chimes_table_bins.N_cool_4d_ElectronDensities,
                      chimes_table_bins.N_cool_4d_HIIDensities)];
          }
        }
      }

      for (j = 0; j < chimes_table_bins.N_cool_hiT_4d_Temperatures; j++)
        chimes_table_cooling.rates_hiT_4d[chimes_flatten_index_2d(
            incl_index, j, chimes_table_bins.N_cool_hiT_4d_Temperatures)] =
            chimes_master_cooling.rates_hiT_4d[chimes_flatten_index_2d(
                i, j, chimes_table_bins.N_cool_hiT_4d_Temperatures)];

      incl_index++;
    }
  }

  for (i = 0; i < chimes_table_bins.N_Temperatures; i++) {
    for (j = 0; j < chimes_table_bins.N_Psi; j++) {
      chimes_table_cooling.photoelectric_heating[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Psi)] =
          chimes_master_cooling.photoelectric_heating[chimes_flatten_index_2d(
              i, j, chimes_table_bins.N_Psi)];
      chimes_table_cooling.grain_recombination[chimes_flatten_index_2d(
          i, j, chimes_table_bins.N_Psi)] =
          chimes_master_cooling.grain_recombination[chimes_flatten_index_2d(
              i, j, chimes_table_bins.N_Psi)];
    }

    chimes_table_cooling.gas_grain_transfer[i] =
        chimes_master_cooling.gas_grain_transfer[i];
  }

  for (i = 0; i < chimes_table_bins.N_mol_cool_Temperatures; i++) {
    chimes_table_cooling.H2_cool_lowDens_H2[i] =
        chimes_master_cooling.H2_cool_lowDens_H2[i];
    chimes_table_cooling.H2_cool_lowDens_HI[i] =
        chimes_master_cooling.H2_cool_lowDens_HI[i];
    chimes_table_cooling.H2_cool_lowDens_HII[i] =
        chimes_master_cooling.H2_cool_lowDens_HII[i];
    chimes_table_cooling.H2_cool_lowDens_HeI[i] =
        chimes_master_cooling.H2_cool_lowDens_HeI[i];
    chimes_table_cooling.H2_cool_lowDens_elec[i] =
        chimes_master_cooling.H2_cool_lowDens_elec[i];
    chimes_table_cooling.H2_cool_LTE[i] = chimes_master_cooling.H2_cool_LTE[i];
  }

  if (myGlobalVars->redshift_dependent_UVB_index >= 0) {
    if (myGlobalVars->redshift_dependent_UVB_index >= myGlobalVars->N_spectra) {
      printf(
          "CHIMES ERROR: redshift_dependent_UVB_index = %d, N_spectra = %d.\n",
          myGlobalVars->redshift_dependent_UVB_index, myGlobalVars->N_spectra);
      exit(EXIT_FAILURE);
    }

    allocate_redshift_dependent_UVB_memory(myGlobalVars);

    /* We need to record the element_incl flags
     * from the full photoion reaction lists.
     * These will be used to re-construct the
     * reduced reaction lists every time we
     * load a new UVB spectrum. */
    N_reactions_all = chimes_master_photoion_fuv.N_reactions[1];
    chimes_table_redshift_dependent_UVB.photoion_fuv_element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));

    for (i = 0; i < N_reactions_all; i++) {
      for (j = 0; j < 9; j++)
        chimes_table_redshift_dependent_UVB
            .photoion_fuv_element_incl[chimes_flatten_index_2d(i, j, 9)] =
            chimes_master_photoion_fuv
                .element_incl[chimes_flatten_index_2d(i, j, 9)];
    }

    N_reactions_all = chimes_master_photoion_euv.N_reactions[1];
    chimes_table_redshift_dependent_UVB.photoion_euv_element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));

    for (i = 0; i < N_reactions_all; i++) {
      for (j = 0; j < 9; j++)
        chimes_table_redshift_dependent_UVB
            .photoion_euv_element_incl[chimes_flatten_index_2d(i, j, 9)] =
            chimes_master_photoion_euv
                .element_incl[chimes_flatten_index_2d(i, j, 9)];
    }

    N_reactions_all = chimes_master_photoion_auger_fuv.N_reactions[1];
    chimes_table_redshift_dependent_UVB.photoion_auger_fuv_element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));

    for (i = 0; i < N_reactions_all; i++) {
      for (j = 0; j < 9; j++)
        chimes_table_redshift_dependent_UVB
            .photoion_auger_fuv_element_incl[chimes_flatten_index_2d(i, j, 9)] =
            chimes_master_photoion_auger_fuv
                .element_incl[chimes_flatten_index_2d(i, j, 9)];
    }

    N_reactions_all = chimes_master_photoion_auger_euv.N_reactions[1];
    chimes_table_redshift_dependent_UVB.photoion_auger_euv_element_incl =
        (int *)malloc(N_reactions_all * 9 * sizeof(int));

    for (i = 0; i < N_reactions_all; i++) {
      for (j = 0; j < 9; j++)
        chimes_table_redshift_dependent_UVB
            .photoion_auger_euv_element_incl[chimes_flatten_index_2d(i, j, 9)] =
            chimes_master_photoion_auger_euv
                .element_incl[chimes_flatten_index_2d(i, j, 9)];
    }
  }

  /**********************************
   ** Free memory in master tables **
   **********************************/
  free(chimes_master_T_dependent.reactants);
  free(chimes_master_T_dependent.products);
  free(chimes_master_T_dependent.element_incl);
  free(chimes_master_T_dependent.rates);
  free(chimes_master_T_dependent.molecular_flag);

  free(chimes_master_constant.reactants);
  free(chimes_master_constant.products);
  free(chimes_master_constant.element_incl);
  free(chimes_master_constant.molecular_flag);
  free(chimes_master_constant.rates);

  free(chimes_master_recombination_AB.reactants);
  free(chimes_master_recombination_AB.products);
  free(chimes_master_recombination_AB.element_incl);
  free(chimes_master_recombination_AB.molecular_flag);
  free(chimes_master_recombination_AB.rates);

  free(chimes_master_grain_recombination.reactants);
  free(chimes_master_grain_recombination.products);
  free(chimes_master_grain_recombination.element_incl);
  free(chimes_master_grain_recombination.rates);

  free(chimes_master_cosmic_ray.secondary_ratio);
  free(chimes_master_cosmic_ray.secondary_base_reaction);
  free(chimes_master_cosmic_ray.reactants);
  free(chimes_master_cosmic_ray.products);
  free(chimes_master_cosmic_ray.element_incl);
  free(chimes_master_cosmic_ray.rates);
  free(chimes_master_cosmic_ray.molecular_flag);

  free(chimes_master_CO_cosmic_ray.reactants);
  free(chimes_master_CO_cosmic_ray.products);
  free(chimes_master_CO_cosmic_ray.element_incl);
  free(chimes_master_CO_cosmic_ray.rates);

  if (myGlobalVars->N_spectra > 0) {
    free(chimes_master_photoion_fuv.reactants);
    free(chimes_master_photoion_fuv.products);
    free(chimes_master_photoion_fuv.element_incl);
    free(chimes_master_photoion_fuv.gamma);
    free(chimes_master_photoion_fuv.sigmaPhot);
    free(chimes_master_photoion_fuv.epsilonPhot);

    free(chimes_master_photoion_euv.reactants);
    free(chimes_master_photoion_euv.products);
    free(chimes_master_photoion_euv.element_incl);
    free(chimes_master_photoion_euv.molecular_flag);
    free(chimes_master_photoion_euv.E_thresh);
    free(chimes_master_photoion_euv.sigmaPhot);
    free(chimes_master_photoion_euv.shieldFactor_1D);
    free(chimes_master_photoion_euv.shieldFactor_2D);

    free(chimes_master_photoion_auger_fuv.reactants);
    free(chimes_master_photoion_auger_fuv.products);
    free(chimes_master_photoion_auger_fuv.element_incl);
    free(chimes_master_photoion_auger_fuv.base_reaction);
    free(chimes_master_photoion_auger_fuv.sigmaPhot);

    free(chimes_master_photoion_auger_euv.reactants);
    free(chimes_master_photoion_auger_euv.products);
    free(chimes_master_photoion_auger_euv.element_incl);
    free(chimes_master_photoion_auger_euv.base_reaction);
    free(chimes_master_photoion_auger_euv.sigmaPhot);

    free(chimes_master_photodissoc_group1.reactants);
    free(chimes_master_photodissoc_group1.products);
    free(chimes_master_photodissoc_group1.element_incl);
    free(chimes_master_photodissoc_group1.gamma);
    free(chimes_master_photodissoc_group1.rates);
    free(chimes_master_photodissoc_group1.molecular_flag);

    free(chimes_master_photodissoc_group2.reactants);
    free(chimes_master_photodissoc_group2.products);
    free(chimes_master_photodissoc_group2.element_incl);
    free(chimes_master_photodissoc_group2.gamma_coeff);
    free(chimes_master_photodissoc_group2.rates);

    free(chimes_master_CO_photodissoc.reactants);
    free(chimes_master_CO_photodissoc.products);
    free(chimes_master_CO_photodissoc.element_incl);
    free(chimes_master_CO_photodissoc.gamma);
    free(chimes_master_CO_photodissoc.rates);
    free(chimes_master_CO_photodissoc.self_shielding);
  }

  free(chimes_master_cooling.rates);
  free(chimes_master_cooling.coolants);

  free(chimes_master_cooling.rates_2d);
  free(chimes_master_cooling.rates_hiT_2d);
  free(chimes_master_cooling.coolants_2d);

  free(chimes_master_cooling.rates_4d);
  free(chimes_master_cooling.rates_hiT_4d);
  free(chimes_master_cooling.coolants_4d);

  free(chimes_master_cooling.photoelectric_heating);
  free(chimes_master_cooling.grain_recombination);
  free(chimes_master_cooling.gas_grain_transfer);
  free(chimes_master_cooling.H2_cool_lowDens_H2);
  free(chimes_master_cooling.H2_cool_lowDens_HI);
  free(chimes_master_cooling.H2_cool_lowDens_HII);
  free(chimes_master_cooling.H2_cool_lowDens_HeI);
  free(chimes_master_cooling.H2_cool_lowDens_elec);
  free(chimes_master_cooling.H2_cool_LTE);

  return;
}

/**
 * @brief Read cross sections tables.
 *
 * Reads in all of the photoionisation cross sections
 * tables, one for each spectrum.
 *
 * @param my_table_bins The #chimes_table_bins_struct struct.
 * @param my_photoion_fuv The #chimes_photoion_fuv_struct struct.
 * @param my_photoion_euv The #chimes_photoion_euv_struct struct.
 * @param my_photoion_auger_fuv The #chimes_photoion_auger_fuv_struct struct.
 * @param my_photoion_auger_euv The #chimes_photoion_auger_euv_struct struct.
 * @param my_spectra The #chimes_spectra_struct struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void read_cross_sections_tables(
    struct chimes_table_bins_struct *my_table_bins,
    struct chimes_photoion_fuv_struct *my_photoion_fuv,
    struct chimes_photoion_euv_struct *my_photoion_euv,
    struct chimes_photoion_auger_fuv_struct *my_photoion_auger_fuv,
    struct chimes_photoion_auger_euv_struct *my_photoion_auger_euv,
    struct chimes_spectra_struct *my_spectra,
    struct globalVariables *myGlobalVars) {
  char fname[1000];
  int i, j, k, l, m, rank;
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  hid_t file_id, dataset, dataspace_id, memspace_id;
  float *array_buffer_float;

  /* If using a redshift-dependent UVB,
   * read in the redshift bins. */
  if (myGlobalVars->redshift_dependent_UVB_index >= 0) {
    sprintf(
        fname, "%s/redshifts.hdf5",
        myGlobalVars
            ->PhotoIonTablePath[myGlobalVars->redshift_dependent_UVB_index]);
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    dataset = H5Dopen1(file_id, "N_redshifts");
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            &(chimes_table_redshift_dependent_UVB.N_redshifts));
    H5Dclose(dataset);

    array_buffer_float = (float *)malloc(
        chimes_table_redshift_dependent_UVB.N_redshifts * sizeof(float));
    chimes_table_redshift_dependent_UVB.redshift_bins = (ChimesFloat *)malloc(
        chimes_table_redshift_dependent_UVB.N_redshifts * sizeof(ChimesFloat));

    dataset = H5Dopen1(file_id, "redshift_bins");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (i = 0; i < chimes_table_redshift_dependent_UVB.N_redshifts; i++)
      chimes_table_redshift_dependent_UVB.redshift_bins[i] =
          (ChimesFloat)array_buffer_float[i];

    free(array_buffer_float);

    H5Fclose(file_id);
  }

  /* Read in the column density bins.
   * These are the same for all cross
   * sections files, so take from the
   * first cross sections file. */
  if (myGlobalVars->redshift_dependent_UVB_index == 0)
    sprintf(fname, "%s/z%.3f_cross_sections.hdf5",
            myGlobalVars->PhotoIonTablePath[0],
            chimes_table_redshift_dependent_UVB.redshift_bins[0]);
  else
    sprintf(fname, "%s", myGlobalVars->PhotoIonTablePath[0]);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) {
    printf("CHIMES ERROR: unable to open cross sections file: %s\n", fname);
    exit(EXIT_FAILURE);
  }

  dataset = H5Dopen1(file_id, "TableBins/N_Column_densities");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &(my_table_bins->N_Column_densities));
  H5Dclose(dataset);

  my_table_bins->Column_densities = (ChimesFloat *)malloc(
      my_table_bins->N_Column_densities * sizeof(ChimesFloat));
  array_buffer_float =
      (float *)malloc(my_table_bins->N_Column_densities * sizeof(float));

  dataset = H5Dopen1(file_id, "TableBins/Column_densities");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  for (i = 0; i < my_table_bins->N_Column_densities; i++)
    my_table_bins->Column_densities[i] = (ChimesFloat)array_buffer_float[i];

  free(array_buffer_float);

  H5Fclose(file_id);

  /* Allocate memory to the cross sections
   * tables in the given structures. */
  my_photoion_fuv->sigmaPhot = (ChimesFloat *)malloc(
      myGlobalVars->N_spectra * my_photoion_fuv->N_reactions[1] *
      sizeof(ChimesFloat));
  my_photoion_fuv->epsilonPhot = (ChimesFloat *)malloc(
      myGlobalVars->N_spectra * my_photoion_fuv->N_reactions[1] *
      sizeof(ChimesFloat));

  my_photoion_euv->sigmaPhot = (ChimesFloat *)malloc(
      myGlobalVars->N_spectra * my_photoion_euv->N_reactions[1] *
      sizeof(ChimesFloat));
  my_photoion_euv->shieldFactor_1D = (ChimesFloat *)malloc(
      my_photoion_euv->N_reactions[1] * myGlobalVars->N_spectra * 3 *
      chimes_table_bins.N_Column_densities * sizeof(ChimesFloat));
  my_photoion_euv->shieldFactor_2D = (ChimesFloat *)malloc(
      my_photoion_euv->N_reactions[1] * myGlobalVars->N_spectra * 6 *
      chimes_table_bins.N_Column_densities *
      chimes_table_bins.N_Column_densities * sizeof(ChimesFloat));

  my_photoion_auger_fuv->sigmaPhot = (ChimesFloat *)malloc(
      myGlobalVars->N_spectra * my_photoion_auger_fuv->N_reactions[1] *
      sizeof(ChimesFloat));

  my_photoion_auger_euv->sigmaPhot = (ChimesFloat *)malloc(
      myGlobalVars->N_spectra * my_photoion_auger_euv->N_reactions[1] *
      sizeof(ChimesFloat));

  my_spectra->isotropic_photon_density =
      (ChimesFloat *)malloc(myGlobalVars->N_spectra * sizeof(ChimesFloat));
  my_spectra->G0_parameter =
      (ChimesFloat *)malloc(myGlobalVars->N_spectra * sizeof(ChimesFloat));
  my_spectra->H2_dissocJ =
      (ChimesFloat *)malloc(myGlobalVars->N_spectra * sizeof(ChimesFloat));

  /* Loop through the spectra, reading in
   * the cross sections tables for each one. */
  for (i = 0; i < myGlobalVars->N_spectra; i++) {
    /* If using a redshift-dependent UVB,
     * skip that spectrum index. This will
     * be loaded separately. */
    if (i == myGlobalVars->redshift_dependent_UVB_index) continue;

    sprintf(fname, "%s", myGlobalVars->PhotoIonTablePath[i]);
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0) {
      printf("CHIMES ERROR: unable to open cross sections file: %s\n", fname);
      exit(EXIT_FAILURE);
    }

    // photoion_fuv
    array_buffer_float =
        (float *)malloc(my_photoion_fuv->N_reactions[1] * sizeof(float));

    dataset = H5Dopen1(file_id, "photoion_fuv/sigmaPhot");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (j = 0; j < my_photoion_fuv->N_reactions[1]; j++)
      my_photoion_fuv->sigmaPhot[chimes_flatten_index_2d(
          i, j, my_photoion_fuv->N_reactions[1])] =
          (ChimesFloat)array_buffer_float[j];

    dataset = H5Dopen1(file_id, "photoion_fuv/epsilonPhot");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (j = 0; j < my_photoion_fuv->N_reactions[1]; j++)
      my_photoion_fuv->epsilonPhot[chimes_flatten_index_2d(
          i, j, my_photoion_fuv->N_reactions[1])] =
          (ChimesFloat)array_buffer_float[j];

    free(array_buffer_float);

    // photoion_euv
    array_buffer_float =
        (float *)malloc(my_photoion_euv->N_reactions[1] * sizeof(float));

    dataset = H5Dopen1(file_id, "photoion_euv/sigmaPhot");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (j = 0; j < my_photoion_euv->N_reactions[1]; j++)
      my_photoion_euv->sigmaPhot[chimes_flatten_index_2d(
          i, j, my_photoion_euv->N_reactions[1])] =
          (ChimesFloat)array_buffer_float[j];

    dataset = H5Dopen1(file_id, "photoion_euv/shieldFactor_1D");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
    dims3D[0] = my_photoion_euv->N_reactions[1];
    dims3D[1] = 1;
    dims3D[2] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims3D, NULL);

    for (j = 0; j < 3; j++) {
      for (k = 0; k < my_table_bins->N_Column_densities; k++) {
        offset3D[0] = 0;
        offset3D[1] = j;
        offset3D[2] = k;
        count3D[0] = my_photoion_euv->N_reactions[1];
        count3D[1] = 1;
        count3D[2] = 1;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL,
                            count3D, NULL);

        H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                H5P_DEFAULT, array_buffer_float);

        for (l = 0; l < my_photoion_euv->N_reactions[1]; l++)
          my_photoion_euv->shieldFactor_1D[chimes_flatten_index_4d(
              l, i, j, k, myGlobalVars->N_spectra, 3,
              chimes_table_bins.N_Column_densities)] =
              (ChimesFloat)array_buffer_float[l];
      }
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    dataset = H5Dopen1(file_id, "photoion_euv/shieldFactor_2D");
    dataspace_id = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
    dims4D[0] = my_photoion_euv->N_reactions[1];
    dims4D[1] = 1;
    dims4D[2] = 1;
    dims4D[3] = 1;
    rank = 1;
    memspace_id = H5Screate_simple(rank, dims4D, NULL);

    for (j = 0; j < 6; j++) {
      for (k = 0; k < my_table_bins->N_Column_densities; k++) {
        for (l = 0; l < my_table_bins->N_Column_densities; l++) {
          offset4D[0] = 0;
          offset4D[1] = j;
          offset4D[2] = k;
          offset4D[3] = l;
          count4D[0] = my_photoion_euv->N_reactions[1];
          count4D[1] = 1;
          count4D[2] = 1;
          count4D[3] = 1;
          H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL,
                              count4D, NULL);

          H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                  H5P_DEFAULT, array_buffer_float);

          for (m = 0; m < my_photoion_euv->N_reactions[1]; m++)
            my_photoion_euv->shieldFactor_2D[chimes_flatten_index_5d(
                m, i, j, k, l, myGlobalVars->N_spectra, 6,
                chimes_table_bins.N_Column_densities,
                chimes_table_bins.N_Column_densities)] =
                (ChimesFloat)array_buffer_float[m];
        }
      }
    }
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset);

    free(array_buffer_float);

    // photoion_auger_fuv
    array_buffer_float =
        (float *)malloc(my_photoion_auger_fuv->N_reactions[1] * sizeof(float));
    dataset = H5Dopen1(file_id, "photoion_auger_fuv/sigmaPhot");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (j = 0; j < my_photoion_auger_fuv->N_reactions[1]; j++)
      my_photoion_auger_fuv->sigmaPhot[chimes_flatten_index_2d(
          i, j, my_photoion_auger_fuv->N_reactions[1])] =
          (ChimesFloat)array_buffer_float[j];

    free(array_buffer_float);

    // photoion_auger_euv
    array_buffer_float =
        (float *)malloc(my_photoion_auger_euv->N_reactions[1] * sizeof(float));
    dataset = H5Dopen1(file_id, "photoion_auger_euv/sigmaPhot");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    for (j = 0; j < my_photoion_auger_euv->N_reactions[1]; j++)
      my_photoion_auger_euv->sigmaPhot[chimes_flatten_index_2d(
          i, j, my_photoion_auger_euv->N_reactions[1])] =
          (ChimesFloat)array_buffer_float[j];

    free(array_buffer_float);

    /* We read in the isotropic_photon_density,
     * G0_parameter and H2_dissocJ and store
     * them in the my_spectra structure. We can
     * then either copy these over to myGasVars
     * to use the parameters from the tables,
     * or we can set our own parameters in
     * myGasVars, as required. */
    array_buffer_float = (float *)malloc(sizeof(float));

    dataset = H5Dopen1(file_id, "isotropic_photon_density");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    my_spectra->isotropic_photon_density[i] =
        (ChimesFloat)array_buffer_float[0];

    dataset = H5Dopen1(file_id, "G0_parameter");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    my_spectra->G0_parameter[i] = (ChimesFloat)array_buffer_float[0];

    dataset = H5Dopen1(file_id, "H2_dissocJ");
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            array_buffer_float);
    H5Dclose(dataset);

    my_spectra->H2_dissocJ[i] = (ChimesFloat)array_buffer_float[0];

    free(array_buffer_float);

    H5Fclose(file_id);
  }

  return;
}

/**
 * @brief Set the species index array.
 *
 * This routine determines which species are included
 * in the network, based on the element_included flags.
 * Note that if an element is included but has zero
 * abundance, its corresponding species WILL be given
 * as included by this routine, because they will be
 * present in the abundance array. However, they will
 * NOT be included in the actual rate vector that CVODE
 * will integrate.
 *
 * @param myGlobalVars The #globalVariables struct.
 */
int set_species_index_array(struct globalVariables *myGlobalVars) {
  int i;
  int current_index = 0;

  /* First, include electrons and all
   * H & He ions. */
  for (i = sp_elec; i <= sp_HeIII; i++) {
    myGlobalVars->speciesIndices[i] = current_index;
    current_index += 1;
  }

  /* For each metal, check whether it
   * is included. If it is, add these
   * the speciesIndices, otherwise
   * make that index -1. */
  for (i = sp_CI; i <= sp_Cm; i++) {
    if (myGlobalVars->element_included[0] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_NI; i <= sp_NVIII; i++) {
    if (myGlobalVars->element_included[1] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_OI; i <= sp_Om; i++) {
    if (myGlobalVars->element_included[2] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_NeI; i <= sp_NeXI; i++) {
    if (myGlobalVars->element_included[3] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_MgI; i <= sp_MgXIII; i++) {
    if (myGlobalVars->element_included[4] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_SiI; i <= sp_SiXV; i++) {
    if (myGlobalVars->element_included[5] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_SI; i <= sp_SXVII; i++) {
    if (myGlobalVars->element_included[6] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_CaI; i <= sp_CaXXI; i++) {
    if (myGlobalVars->element_included[7] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  for (i = sp_FeI; i <= sp_FeXXVII; i++) {
    if (myGlobalVars->element_included[8] == 1) {
      myGlobalVars->speciesIndices[i] = current_index;
      current_index += 1;
    } else
      myGlobalVars->speciesIndices[i] = -1;
  }

  /* Determine which molecules are present. */
  for (i = sp_H2; i <= sp_H3p; i++) {
    myGlobalVars->speciesIndices[i] = current_index;
    current_index += 1;
  }

  if (myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_OH] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_H2O] = current_index;
    current_index += 1;
  } else {
    myGlobalVars->speciesIndices[sp_OH] = -1;
    myGlobalVars->speciesIndices[sp_H2O] = -1;
  }

  if (myGlobalVars->element_included[0] == 1) {
    myGlobalVars->speciesIndices[sp_C2] = current_index;
    current_index += 1;
  } else
    myGlobalVars->speciesIndices[sp_C2] = -1;

  if (myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_O2] = current_index;
    current_index += 1;
  } else
    myGlobalVars->speciesIndices[sp_O2] = -1;

  if (myGlobalVars->element_included[0] == 1 &&
      myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_HCOp] = current_index;
    current_index += 1;
  } else
    myGlobalVars->speciesIndices[sp_HCOp] = -1;

  if (myGlobalVars->element_included[0] == 1) {
    myGlobalVars->speciesIndices[sp_CH] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_CH2] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_CH3p] = current_index;
    current_index += 1;
  } else {
    myGlobalVars->speciesIndices[sp_CH] = -1;
    myGlobalVars->speciesIndices[sp_CH2] = -1;
    myGlobalVars->speciesIndices[sp_CH3p] = -1;
  }

  if (myGlobalVars->element_included[0] == 1 &&
      myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_CO] = current_index;
    current_index += 1;
  } else
    myGlobalVars->speciesIndices[sp_CO] = -1;

  if (myGlobalVars->element_included[0] == 1) {
    myGlobalVars->speciesIndices[sp_CHp] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_CH2p] = current_index;
    current_index += 1;
  } else {
    myGlobalVars->speciesIndices[sp_CHp] = -1;
    myGlobalVars->speciesIndices[sp_CH2p] = -1;
  }

  if (myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_OHp] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_H2Op] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_H3Op] = current_index;
    current_index += 1;
  } else {
    myGlobalVars->speciesIndices[sp_OHp] = -1;
    myGlobalVars->speciesIndices[sp_H2Op] = -1;
    myGlobalVars->speciesIndices[sp_H3Op] = -1;
  }

  if (myGlobalVars->element_included[0] == 1 &&
      myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_COp] = current_index;
    current_index += 1;
    myGlobalVars->speciesIndices[sp_HOCp] = current_index;
    current_index += 1;
  } else {
    myGlobalVars->speciesIndices[sp_COp] = -1;
    myGlobalVars->speciesIndices[sp_HOCp] = -1;
  }

  if (myGlobalVars->element_included[2] == 1) {
    myGlobalVars->speciesIndices[sp_O2p] = current_index;
    current_index += 1;
  } else
    myGlobalVars->speciesIndices[sp_O2p] = -1;

  return current_index;
}

/**
 * @brief Initialise the CHIMES network.
 *
 * Calls the various routines that read in the reaction
 * data tables and allocate memory to the various arrays
 * that store these tables.
 *
 * @param myGlobalVars The #globalVariables struct.
 */
void init_chimes(struct globalVariables *myGlobalVars) {
  char fname[520];

  /* Print the float precision used in CHIMES 
   * and in the Sundials library. */ 
  if (sizeof(ChimesFloat) == sizeof(float)) 
    printf("CHIMES built in single precision.\n"); 
  else if (sizeof(ChimesFloat) == sizeof(double)) 
    printf("CHIMES built in double precision.\n"); 
  else
    printf("CHIMES precision not recognised.\n"); 

  if (sizeof(realtype) == sizeof(float)) 
    printf("Sundials built in single precision.\n"); 
  else if (sizeof(realtype) == sizeof(double)) 
    printf("Sundials built in double precision.\n"); 
  else 
    printf("Sundials precision not recognised.\n"); 

  /* Check that both CHIMES and the Sundials
   * library have been built with the same
   * float precision. */
  if (sizeof(ChimesFloat) != sizeof(realtype)) {
    printf("CHIMES ERROR: CHIMES and Sundials have to use the same float precision.\n"); 
    exit(EXIT_FAILURE);
  }

  myGlobalVars->totalNumberOfSpecies = set_species_index_array(myGlobalVars);

  /* Read in the main CHIMES data file and
   * construct the network */
  initialise_main_data(myGlobalVars);

  // Read in tables of equilibrium abundances
  if (myGlobalVars->use_redshift_dependent_eqm_tables == 1) {
    if (myGlobalVars->redshift_dependent_UVB_index < 0) {
      printf(
          "CHIMES ERROR: use_redshift_dependent_eqm_tables = %d, but "
          "redshift_dependent_UVB_index = %d. The redshift dependent UVB needs "
          "to be enabled before you can use the redshift-dependent eqm tables. "
          "\n",
          myGlobalVars->use_redshift_dependent_eqm_tables,
          myGlobalVars->redshift_dependent_UVB_index);
      exit(EXIT_FAILURE);
    }
    sprintf(fname, "%s/z%.3f_eqm.hdf5", myGlobalVars->EqAbundanceTablePath,
            chimes_table_redshift_dependent_UVB.redshift_bins[0]);
    allocate_eqm_table_memory(fname, &chimes_table_eqm_abundances,
                              myGlobalVars);
    allocate_eqm_table_memory(
        fname, &(chimes_table_redshift_dependent_UVB.eqm_abundances[0]),
        myGlobalVars);
    allocate_eqm_table_memory(
        fname, &(chimes_table_redshift_dependent_UVB.eqm_abundances[1]),
        myGlobalVars);
  } else {
    sprintf(fname, "%s", myGlobalVars->EqAbundanceTablePath);
    allocate_eqm_table_memory(fname, &chimes_table_eqm_abundances,
                              myGlobalVars);
    load_eqm_table(fname, &chimes_table_eqm_abundances, myGlobalVars);
  }

  /* Interpolate UVB to current redshift,
   *  if required. */
  if (myGlobalVars->redshift_dependent_UVB_index >= 0)
    interpolate_redshift_dependent_UVB(myGlobalVars);

  /* Initialise the hybrid cooling
   * pointers to NULL. */
  myGlobalVars->hybrid_data = NULL;
  myGlobalVars->hybrid_cooling_fn = NULL;
  myGlobalVars->allocate_gas_hybrid_data_fn = NULL;
  myGlobalVars->free_gas_hybrid_data_fn = NULL;
}

/**
 * @brief Allocate memory for the gasVars arrays.
 *
 * Allocates memory for the abundances, isotropic_photon_density,
 * G0_parameter and H2_dissocJ arrays in the #gasVariables struct.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void allocate_gas_abundances_memory(struct gasVariables *myGasVars,
                                    struct globalVariables *myGlobalVars) {
  myGasVars->abundances = (ChimesFloat *)malloc(
      myGlobalVars->totalNumberOfSpecies * sizeof(ChimesFloat));

  /* We also allocate memory for radiation fields here. */
  if (myGlobalVars->N_spectra > 0) {
    myGasVars->isotropic_photon_density =
        (ChimesFloat *)malloc(myGlobalVars->N_spectra * sizeof(ChimesFloat));
    myGasVars->G0_parameter =
        (ChimesFloat *)malloc(myGlobalVars->N_spectra * sizeof(ChimesFloat));
    myGasVars->H2_dissocJ =
        (ChimesFloat *)malloc(myGlobalVars->N_spectra * sizeof(ChimesFloat));
  }

  if (myGlobalVars->hybrid_cooling_mode == 1)
    (*myGlobalVars->allocate_gas_hybrid_data_fn)(myGasVars);
  else
    myGasVars->hybrid_data = NULL;
}

/**
 * @brief Free memory for the gasVars arrays.
 *
 * Frees the memory for the abundances, isotropic_photon_density,
 * G0_parameter and H2_dissocJ arrays in the #gasVariables struct.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void free_gas_abundances_memory(struct gasVariables *myGasVars,
                                struct globalVariables *myGlobalVars) {
  free(myGasVars->abundances);

  if (myGlobalVars->N_spectra > 0) {
    free(myGasVars->isotropic_photon_density);
    free(myGasVars->G0_parameter);
    free(myGasVars->H2_dissocJ);
  }

  if (myGlobalVars->hybrid_cooling_mode == 1)
    (*myGlobalVars->free_gas_hybrid_data_fn)(myGasVars);
}

/**
 * @brief Initialise the abundance array.
 *
 * Sets the initial abundances in the abundance array, according
 * to the InitIonState parameter in the #gasVariables struct.
 * Each element is set to be in the ionisation state given by this
 * parameter. If the parameter exceeds the maximum ionisation state
 * of that element, then it is set to the highest possible state.
 * For example, if InitIonState == 0, all elements are set to be
 * fully neutral. If InitIonState == 2, all elements are set to be
 * doubly ionised, except Hydrogen which is set to be singly ionised.
 * Molecules are always initialised to zero in this routine.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void initialise_gas_abundances(struct gasVariables *myGasVars,
                               struct globalVariables *myGlobalVars) {
  int i, init_ion_state;
  /* mode determines in what state
   * we initially set the gas. */
  if (myGasVars->InitIonState >= 0 && myGasVars->InitIonState <= 26)
    init_ion_state = myGasVars->InitIonState;
  else {
    printf(
        "WARNING: initialise_gas_abundances() mode not recognised. Assuming "
        "fully neutral\n");
    fflush(stdout);
    init_ion_state = 0;
  }

  /* First, set all abundances to zero */
  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    myGasVars->abundances[i] = 0.0;

  /* Now set the abundances of the initial
   * species.  */
  myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
      sp_HI + init_ion_state, sp_HII)]] = 1.0;
  myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] =
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]];

  myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
      sp_HeI + init_ion_state, sp_HeIII)]] = myGasVars->element_abundances[0];
  for (i = 1; i <= 2; i++)
    myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
        myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI + i]] * i;

  if (myGlobalVars->element_included[0] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_CI + init_ion_state, sp_CVII)]] = myGasVars->element_abundances[1];
    for (i = 1; i <= 6; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_CI + i]] * i;
  }

  if (myGlobalVars->element_included[1] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_NI + init_ion_state, sp_NVIII)]] = myGasVars->element_abundances[2];
    for (i = 1; i <= 7; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_NI + i]] * i;
  }

  if (myGlobalVars->element_included[2] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_OI + init_ion_state, sp_OIX)]] = myGasVars->element_abundances[3];
    for (i = 1; i <= 8; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_OI + i]] * i;
  }

  if (myGlobalVars->element_included[3] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_NeI + init_ion_state, sp_NeXI)]] = myGasVars->element_abundances[4];
    for (i = 1; i <= 10; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_NeI + i]] * i;
  }

  if (myGlobalVars->element_included[4] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_MgI + init_ion_state, sp_MgXIII)]] =
        myGasVars->element_abundances[5];
    for (i = 1; i <= 12; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_MgI + i]] * i;
  }

  if (myGlobalVars->element_included[5] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_SiI + init_ion_state, sp_SiXV)]] = myGasVars->element_abundances[6];
    for (i = 1; i <= 14; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_SiI + i]] * i;
  }

  if (myGlobalVars->element_included[6] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_SI + init_ion_state, sp_SXVII)]] = myGasVars->element_abundances[7];
    for (i = 1; i <= 16; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_SI + i]] * i;
  }

  if (myGlobalVars->element_included[7] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_CaI + init_ion_state, sp_CaXXI)]] = myGasVars->element_abundances[8];
    for (i = 1; i <= 20; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_CaI + i]] * i;
  }

  if (myGlobalVars->element_included[8] == 1) {
    myGasVars->abundances[myGlobalVars->speciesIndices[chimes_min(
        sp_FeI + init_ion_state, sp_FeXXVII)]] =
        myGasVars->element_abundances[9];
    for (i = 1; i <= 26; i++)
      myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] +=
          myGasVars->abundances[myGlobalVars->speciesIndices[sp_FeI + i]] * i;
  }

  /* Make sure that the abundances of any elements
   * not included are set to zero. */
  for (i = 0; i < 9; i++)
    if (myGlobalVars->element_included[i] != 1)
      myGasVars->element_abundances[i + 1] = 0.0;

  // Set constant heating rate to zero.
  myGasVars->constant_heating_rate = 0.0;
}

/**
 * @brief Determines size of the current rates buffer.
 *
 * Counts the total number of reactions in the network,
 * which determines the size of the data buffers that will
 * be needed to store the current values of all the
 * reaction rate coefficients. For some reactions we need
 * to store multiple numbers (e.g. the rate coefficient itself,
 * and the rate after multiplying by the reactant abundances).
 *
 * @param buffer_size number of single/double floats in the data buffer.
 * @param myGlobalVars The #globalVariables struct.
 */
void determine_current_rates_buffer_size(int *buffer_size,
                                         struct globalVariables *myGlobalVars) {
  *buffer_size = 0;

  *buffer_size += chimes_table_T_dependent.N_reactions[1] * 2;
  *buffer_size += chimes_table_constant.N_reactions[1];
  *buffer_size += chimes_table_recombination_AB.N_reactions[1] * 2;
  *buffer_size += chimes_table_grain_recombination.N_reactions[1] * 2;
  *buffer_size += chimes_table_cosmic_ray.N_reactions[1];
  *buffer_size += chimes_table_CO_cosmic_ray.N_reactions[1] * 2;
  *buffer_size += chimes_table_H2_collis_dissoc.N_reactions[1] * 4;
  *buffer_size += chimes_table_photoion_fuv.N_reactions[1] * 4;
  *buffer_size += chimes_table_photoion_euv.N_reactions[1] * 3;
  *buffer_size +=
      chimes_table_photoion_euv.N_reactions[1] * myGlobalVars->N_spectra * 2;
  *buffer_size += chimes_table_photoion_auger_fuv.N_reactions[1] * 2;
  *buffer_size += chimes_table_photoion_auger_euv.N_reactions[1] * 2;
  *buffer_size += chimes_table_photodissoc_group1.N_reactions[1] * 3;
  *buffer_size += chimes_table_photodissoc_group2.N_reactions[1] * 2;
  *buffer_size += chimes_table_H2_photodissoc.N_reactions[1] * 3;
  *buffer_size += chimes_table_CO_photodissoc.N_reactions[1] * 3;
  *buffer_size += chimes_table_cooling.N_coolants;
  *buffer_size += chimes_table_cooling.N_coolants_2d;
  *buffer_size += chimes_table_cooling.N_coolants_4d;

  return;
}

/**
 * @brief Allocates memory for the current rates.
 *
 * Allocates memory for the various arrays in the
 * #chimes_current_rates_struct struct that will
 * record the current values of the various
 * reaction rates.
 *
 * @param chimes_current_rates The #chimes_current_rates_struct struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void allocate_current_rates_memory(
    struct chimes_current_rates_struct *chimes_current_rates,
    struct globalVariables *myGlobalVars) {
  int buffer_size;
  int buffer_position;
  int N_reactions;

  determine_current_rates_buffer_size(&buffer_size, myGlobalVars);
  buffer_position = 0;

  chimes_current_rates->data_buffer =
      (ChimesFloat *)malloc(buffer_size * sizeof(ChimesFloat));

  // T_dependent
  N_reactions = chimes_table_T_dependent.N_reactions[1];

  chimes_current_rates->T_dependent_rate_coefficient =
      chimes_current_rates->data_buffer;
  buffer_position += N_reactions;

  chimes_current_rates->T_dependent_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // constant
  N_reactions = chimes_table_constant.N_reactions[1];

  chimes_current_rates->constant_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // recombination_AB
  N_reactions = chimes_table_recombination_AB.N_reactions[1];

  chimes_current_rates->recombination_AB_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->recombination_AB_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // grain_recombination
  N_reactions = chimes_table_grain_recombination.N_reactions[1];

  chimes_current_rates->grain_recombination_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->grain_recombination_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // cosmic_ray
  N_reactions = chimes_table_cosmic_ray.N_reactions[1];

  chimes_current_rates->cosmic_ray_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // CO_cosmic_ray
  N_reactions = chimes_table_CO_cosmic_ray.N_reactions[1];

  chimes_current_rates->CO_cosmic_ray_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->CO_cosmic_ray_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // H2_collis_dissoc
  N_reactions = chimes_table_H2_collis_dissoc.N_reactions[1];

  chimes_current_rates->H2_collis_dissoc_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->H2_collis_dissoc_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->H2_collis_dissoc_log_k0 =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->H2_collis_dissoc_log_kLTE =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // photoion_fuv
  N_reactions = chimes_table_photoion_fuv.N_reactions[1];

  chimes_current_rates->photoion_fuv_shield_factor =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_fuv_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_fuv_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_fuv_heat_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // photoion_euv
  N_reactions = chimes_table_photoion_euv.N_reactions[1];

  chimes_current_rates->photoion_euv_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_euv_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_euv_heat_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_euv_shield_factor =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions * myGlobalVars->N_spectra;

  chimes_current_rates->photoion_euv_epsilon =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions * myGlobalVars->N_spectra;

  // photoion_auger_fuv
  N_reactions = chimes_table_photoion_auger_fuv.N_reactions[1];

  chimes_current_rates->photoion_auger_fuv_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_auger_fuv_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // photoion_auger_euv
  N_reactions = chimes_table_photoion_auger_euv.N_reactions[1];

  chimes_current_rates->photoion_auger_euv_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photoion_auger_euv_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // photodissoc_group1
  N_reactions = chimes_table_photodissoc_group1.N_reactions[1];

  chimes_current_rates->photodissoc_group1_shield_factor =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photodissoc_group1_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photodissoc_group1_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // photodissoc_group2
  N_reactions = chimes_table_photodissoc_group2.N_reactions[1];

  chimes_current_rates->photodissoc_group2_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->photodissoc_group2_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // H2_photodissoc
  N_reactions = chimes_table_H2_photodissoc.N_reactions[1];

  chimes_current_rates->H2_photodissoc_shield_factor =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->H2_photodissoc_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->H2_photodissoc_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // CO photodissoc
  N_reactions = chimes_table_CO_photodissoc.N_reactions[1];

  chimes_current_rates->CO_photodissoc_shield_factor =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->CO_photodissoc_rate_coefficient =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  chimes_current_rates->CO_photodissoc_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // Cooling
  N_reactions = chimes_table_cooling.N_coolants;

  chimes_current_rates->cooling_rate =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // Cooling 2D
  N_reactions = chimes_table_cooling.N_coolants_2d;

  chimes_current_rates->cooling_rate_2d =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  // Cooling 4D
  N_reactions = chimes_table_cooling.N_coolants_4d;

  chimes_current_rates->cooling_rate_4d =
      chimes_current_rates->data_buffer + buffer_position;
  buffer_position += N_reactions;

  /* Check that the end position in the
   * data buffers matches the expected
   * buffer sizes. */
  if (buffer_position != buffer_size) {
    printf(
        "CHIMES ERROR: in allocate_current_rates_memory(), buffer_position = "
        "%d, buffer_size = %d. \n",
        buffer_position, buffer_size);
    exit(EXIT_FAILURE);
  }

  return;
}

/**
 * @brief Frees memory for the current rates.
 *
 * Frees the memory for the various arrays in the
 * #chimes_current_rates_struct struct.
 *
 * @param chimes_current_rates The #chimes_current_rates_struct struct.
 * @param myGlobalVars The #globalVariables struct.
 */
void free_current_rates_memory(
    struct chimes_current_rates_struct *chimes_current_rates,
    struct globalVariables *myGlobalVars) {
  free(chimes_current_rates->data_buffer);

  return;
}

/**
 * @brief Allocates memory for the redshift-dependent UVB.
 *
 * Allocates memory for the various arrays in the
 * #chimes_redshift_dependent_UVB_struct struct, which will store
 * the cross sections, shielding factors etc. for the two redshift bins
 * bracketing the current redshift.
 *
 * @param myGlobalVars The #globalVariables struct.
 */
void allocate_redshift_dependent_UVB_memory(
    struct globalVariables *myGlobalVars) {
  int N_reactions_all;

  /* Allocate memory to arrays
   * in redshift dependent UVB. */
  N_reactions_all = chimes_table_photoion_fuv.N_reactions[1];
  chimes_table_redshift_dependent_UVB.photoion_fuv_sigmaPhot =
      (ChimesFloat *)malloc(2 * N_reactions_all * sizeof(ChimesFloat));
  chimes_table_redshift_dependent_UVB.photoion_fuv_epsilonPhot =
      (ChimesFloat *)malloc(2 * N_reactions_all * sizeof(ChimesFloat));

  N_reactions_all = chimes_table_photoion_euv.N_reactions[1];
  chimes_table_redshift_dependent_UVB.photoion_euv_sigmaPhot =
      (ChimesFloat *)malloc(2 * N_reactions_all * sizeof(ChimesFloat));
  chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_1D =
      (ChimesFloat *)malloc(N_reactions_all * 2 * 3 *
                            chimes_table_bins.N_Column_densities *
                            sizeof(ChimesFloat));
  chimes_table_redshift_dependent_UVB.photoion_euv_shieldFactor_2D =
      (ChimesFloat *)malloc(
          N_reactions_all * 2 * 6 * chimes_table_bins.N_Column_densities *
          chimes_table_bins.N_Column_densities * sizeof(ChimesFloat));

  N_reactions_all = chimes_table_photoion_auger_fuv.N_reactions[1];
  chimes_table_redshift_dependent_UVB.photoion_auger_fuv_sigmaPhot =
      (ChimesFloat *)malloc(2 * N_reactions_all * sizeof(ChimesFloat));

  N_reactions_all = chimes_table_photoion_auger_euv.N_reactions[1];
  chimes_table_redshift_dependent_UVB.photoion_auger_euv_sigmaPhot =
      (ChimesFloat *)malloc(2 * N_reactions_all * sizeof(ChimesFloat));

  /* Initialise the current redshift
   * indices to -1, as they haven't
   * yet been read in. */
  chimes_table_redshift_dependent_UVB.z_index_low = -1;
  chimes_table_redshift_dependent_UVB.z_index_hi = -1;
}

/**
 * @brief Load the redshift-dependent UVB.
 *
 * Loads the photoionisation cross-sections data tables for
 * the UVB at the specified redshift.
 *
 * @param redshift The redshift of the required UVB table.
 * @param bin_index Indicates the low- (0) or high- (1) redshift bin.
 * @param myGlobalVars The #globalVariables struct.
 */
void load_redshift_dependent_UVB(ChimesFloat redshift, int bin_index,
                                 struct globalVariables *myGlobalVars) {
  char fname[520];
  hid_t file_id, dataset, dataspace_id, memspace_id;
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int i, j, k, l, rank, N_reactions_all, incl_index;

  float *array_buffer_float;
  int spectrum_index = myGlobalVars->redshift_dependent_UVB_index;

  if (!((bin_index == 0) || (bin_index == 1))) {
    printf(
        "CHIMES ERROR: load_redshift_dependent_UVB() called with bin_index == "
        "%d. Allowed values are 0 or 1.\n",
        bin_index);
    exit(EXIT_FAILURE);
  }

  sprintf(fname, "%s/z%.3f_cross_sections.hdf5",
          myGlobalVars->PhotoIonTablePath[spectrum_index], redshift);
  printf("CHIMES: reading redshift-dependent UVB cross sections table: %s\n",
         fname);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0) {
    printf("CHIMES ERROR: unable to open cross sections file: %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // photoion_fuv
  dataset = H5Dopen1(file_id, "photoion_fuv/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &N_reactions_all);
  H5Dclose(dataset);

  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  dataset = H5Dopen1(file_id, "photoion_fuv/sigmaPhot");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  incl_index = 0;
  for (i = 0; i < N_reactions_all; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_table_redshift_dependent_UVB.photoion_fuv_element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_redshift_dependent_UVB
          .photoion_fuv_sigmaPhot[chimes_flatten_index_2d(
              bin_index, incl_index,
              chimes_table_photoion_fuv.N_reactions[1])] =
          (ChimesFloat)array_buffer_float[i];
      incl_index += 1;
    }
  }

  dataset = H5Dopen1(file_id, "photoion_fuv/epsilonPhot");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  incl_index = 0;
  for (i = 0; i < N_reactions_all; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_table_redshift_dependent_UVB.photoion_fuv_element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_redshift_dependent_UVB
          .photoion_fuv_epsilonPhot[chimes_flatten_index_2d(
              bin_index, incl_index,
              chimes_table_photoion_fuv.N_reactions[1])] =
          (ChimesFloat)array_buffer_float[i];
      incl_index += 1;
    }
  }

  free(array_buffer_float);

  // photoion_euv
  dataset = H5Dopen1(file_id, "photoion_euv/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &N_reactions_all);
  H5Dclose(dataset);

  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  dataset = H5Dopen1(file_id, "photoion_euv/sigmaPhot");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  incl_index = 0;
  for (i = 0; i < N_reactions_all; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_table_redshift_dependent_UVB.photoion_euv_element_incl, i,
            myGlobalVars->element_included)) {
      chimes_table_redshift_dependent_UVB
          .photoion_euv_sigmaPhot[chimes_flatten_index_2d(
              bin_index, incl_index,
              chimes_table_photoion_euv.N_reactions[1])] =
          (ChimesFloat)array_buffer_float[i];
      incl_index += 1;
    }
  }

  dataset = H5Dopen1(file_id, "photoion_euv/shieldFactor_1D");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
  dims3D[0] = N_reactions_all;
  dims3D[1] = 1;
  dims3D[2] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims3D, NULL);

  for (j = 0; j < 3; j++) {
    for (k = 0; k < chimes_table_bins.N_Column_densities; k++) {
      offset3D[0] = 0;
      offset3D[1] = j;
      offset3D[2] = k;
      count3D[0] = N_reactions_all;
      count3D[1] = 1;
      count3D[2] = 1;
      H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D,
                          NULL);

      H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
              array_buffer_float);

      incl_index = 0;
      for (i = 0; i < N_reactions_all; i++) {
        if (compare_element_incl_arrays_flatten(
                chimes_table_redshift_dependent_UVB.photoion_euv_element_incl,
                i, myGlobalVars->element_included)) {
          chimes_table_redshift_dependent_UVB
              .photoion_euv_shieldFactor_1D[chimes_flatten_index_4d(
                  incl_index, bin_index, j, k, 2, 3,
                  chimes_table_bins.N_Column_densities)] =
              (ChimesFloat)array_buffer_float[i];
          incl_index += 1;
        }
      }
    }
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  dataset = H5Dopen1(file_id, "photoion_euv/shieldFactor_2D");
  dataspace_id = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
  dims4D[0] = N_reactions_all;
  dims4D[1] = 1;
  dims4D[2] = 1;
  dims4D[3] = 1;
  rank = 1;
  memspace_id = H5Screate_simple(rank, dims4D, NULL);

  for (j = 0; j < 6; j++) {
    for (k = 0; k < chimes_table_bins.N_Column_densities; k++) {
      for (l = 0; l < chimes_table_bins.N_Column_densities; l++) {
        offset4D[0] = 0;
        offset4D[1] = j;
        offset4D[2] = k;
        offset4D[3] = l;
        count4D[0] = N_reactions_all;
        count4D[1] = 1;
        count4D[2] = 1;
        count4D[3] = 1;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL,
                            count4D, NULL);

        H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id,
                H5P_DEFAULT, array_buffer_float);

        incl_index = 0;
        for (i = 0; i < N_reactions_all; i++) {
          if (compare_element_incl_arrays_flatten(
                  chimes_table_redshift_dependent_UVB.photoion_euv_element_incl,
                  i, myGlobalVars->element_included)) {
            chimes_table_redshift_dependent_UVB
                .photoion_euv_shieldFactor_2D[chimes_flatten_index_5d(
                    incl_index, bin_index, j, k, l, 2, 6,
                    chimes_table_bins.N_Column_densities,
                    chimes_table_bins.N_Column_densities)] =
                (ChimesFloat)array_buffer_float[i];
            incl_index += 1;
          }
        }
      }
    }
  }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  free(array_buffer_float);

  // photoion_auger_fuv
  dataset = H5Dopen1(file_id, "photoion_auger_fuv/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &N_reactions_all);
  H5Dclose(dataset);

  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  dataset = H5Dopen1(file_id, "photoion_auger_fuv/sigmaPhot");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  incl_index = 0;
  for (i = 0; i < N_reactions_all; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_table_redshift_dependent_UVB.photoion_auger_fuv_element_incl,
            i, myGlobalVars->element_included)) {
      chimes_table_redshift_dependent_UVB
          .photoion_auger_fuv_sigmaPhot[chimes_flatten_index_2d(
              bin_index, incl_index,
              chimes_table_photoion_auger_fuv.N_reactions[1])] =
          (ChimesFloat)array_buffer_float[i];
      incl_index += 1;
    }
  }

  free(array_buffer_float);

  // photoion_auger_euv
  dataset = H5Dopen1(file_id, "photoion_auger_euv/N_reactions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          &N_reactions_all);
  H5Dclose(dataset);

  array_buffer_float = (float *)malloc(N_reactions_all * sizeof(float));

  dataset = H5Dopen1(file_id, "photoion_auger_euv/sigmaPhot");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  incl_index = 0;
  for (i = 0; i < N_reactions_all; i++) {
    if (compare_element_incl_arrays_flatten(
            chimes_table_redshift_dependent_UVB.photoion_auger_euv_element_incl,
            i, myGlobalVars->element_included)) {
      chimes_table_redshift_dependent_UVB
          .photoion_auger_euv_sigmaPhot[chimes_flatten_index_2d(
              bin_index, incl_index,
              chimes_table_photoion_auger_euv.N_reactions[1])] =
          (ChimesFloat)array_buffer_float[i];
      incl_index += 1;
    }
  }

  free(array_buffer_float);

  // General spectrum info
  array_buffer_float = (float *)malloc(sizeof(float));

  dataset = H5Dopen1(file_id, "isotropic_photon_density");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  chimes_table_redshift_dependent_UVB.isotropic_photon_density[bin_index] =
      (ChimesFloat)array_buffer_float[0];

  dataset = H5Dopen1(file_id, "G0_parameter");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  chimes_table_redshift_dependent_UVB.G0_parameter[bin_index] =
      (ChimesFloat)array_buffer_float[0];

  dataset = H5Dopen1(file_id, "H2_dissocJ");
  H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          array_buffer_float);
  H5Dclose(dataset);

  chimes_table_redshift_dependent_UVB.H2_dissocJ[bin_index] =
      (ChimesFloat)array_buffer_float[0];

  free(array_buffer_float);

  if (myGlobalVars->use_redshift_dependent_eqm_tables == 1) {
    sprintf(fname, "%s/z%.3f_eqm.hdf5", myGlobalVars->EqAbundanceTablePath,
            redshift);
    load_eqm_table(
        fname, &(chimes_table_redshift_dependent_UVB.eqm_abundances[bin_index]),
        myGlobalVars);
  }
}

/**
 * Array of species names
 * for the full network.
 */
const char chimes_species_names[CHIMES_TOTSIZE][CHIMES_NAME_STR_LENGTH] = {
    "elec",   "HI",      "HII",    "Hm",      "HeI",   "HeII",   "HeIII",
    "CI",     "CII",     "CIII",   "CIV",     "CV",    "CVI",    "CVII",
    "Cm",     "NI",      "NII",    "NIII",    "NIV",   "NV",     "NVI",
    "NVII",   "NVIII",   "OI",     "OII",     "OIII",  "OIV",    "OV",
    "OVI",    "OVII",    "OVIII",  "OIX",     "Om",    "NeI",    "NeII",
    "NeIII",  "NeIV",    "NeV",    "NeVI",    "NeVII", "NeVIII", "NeIX",
    "NeX",    "NeXI",    "MgI",    "MgII",    "MgIII", "MgIV",   "MgV",
    "MgVI",   "MgVII",   "MgVIII", "MgIX",    "MgX",   "MgXI",   "MgXII",
    "MgXIII", "SiI",     "SiII",   "SiIII",   "SiIV",  "SiV",    "SiVI",
    "SiVII",  "SiVIII",  "SiIX",   "SiX",     "SiXI",  "SiXII",  "SiXIII",
    "SiXIV",  "SiXV",    "SI",     "SII",     "SIII",  "SIV",    "SV",
    "SVI",    "SVII",    "SVIII",  "SIX",     "SX",    "SXI",    "SXII",
    "SXIII",  "SXIV",    "SXV",    "SXVI",    "SXVII", "CaI",    "CaII",
    "CaIII",  "CaIV",    "CaV",    "CaVI",    "CaVII", "CaVIII", "CaIX",
    "CaX",    "CaXI",    "CaXII",  "CaXIII",  "CaXIV", "CaXV",   "CaXVI",
    "CaXVII", "CaXVIII", "CaXIX",  "CaXX",    "CaXXI", "FeI",    "FeII",
    "FeIII",  "FeIV",    "FeV",    "FeVI",    "FeVII", "FeVIII", "FeIX",
    "FeX",    "FeXI",    "FeXII",  "FeXIII",  "FeXIV", "FeXV",   "FeXVI",
    "FeXVII", "FeXVIII", "FeXIX",  "FeXX",    "FeXXI", "FeXXII", "FeXXIII",
    "FeXXIV", "FeXXV",   "FeXXVI", "FeXXVII", "H2",    "H2p",    "H3p",
    "OH",     "H2O",     "C2",     "O2",      "HCOp",  "CH",     "CH2",
    "CH3p",   "CO",      "CHp",    "CH2p",    "OHp",   "H2Op",   "H3Op",
    "COp",    "HOCp",    "O2p"};
