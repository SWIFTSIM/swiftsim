/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
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
 ******************************************************************************/
#ifndef SWIFT_EAGLE_COOL_TABLES_H
#define SWIFT_EAGLE_COOL_TABLES_H

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function
 */

/* Config parameters. */
#include "../config.h"

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cooling_struct.h"
#include "eagle_cool_tables.h"
#include "error.h"
#include "interpolate.h"

/* Names of the elements in the order they are storeed in the files */
static const char* eagle_tables_element_names[9] = {
    "Carbon",  "Nitrogen", "Oxygen",  "Neon", "Magnesium",
    "Silicon", "Sulphur",  "Calcium", "Iron"};


/*
 * @brief Reads in EAGLE table of redshift values
 *
 * @param cooling Cooling data structure
 */
void GetCoolingRedshifts(struct cooling_function_data *cooling) {
  FILE *infile;

  int i = 0;

  char buffer[500], redfilename[500];

  sprintf(redfilename, "%s/redshifts.dat", cooling->cooling_table_path);
  infile = fopen(redfilename, "r");
  if (infile == NULL) puts("GetCoolingRedshifts can't open a file");

  if (fscanf(infile, "%s", buffer) != EOF) {
    cooling->N_Redshifts = atoi(buffer);
    cooling->Redshifts = (float *)malloc(cooling->N_Redshifts * sizeof(float));

    while (fscanf(infile, "%s", buffer) != EOF) {
      cooling->Redshifts[i] = atof(buffer);
      i += 1;
    }
  }
  fclose(infile);
}

/*
 * @brief Reads in EAGLE cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, helium fraction
 * solar element abundances, and elements used to index the cooling tables.
 *
 * @param fname Filepath for cooling table from which to read header
 * @param cooling Cooling data structure
 */
void ReadCoolingHeader(char *fname, struct cooling_function_data *cooling) {
#ifdef HAVE_HDF5
  int i;

  hid_t tempfile_id, dataset;

  herr_t status;

  // read sizes of array dimensions
  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (tempfile_id < 0) {
    error("[ReadCoolingHeader()]: unable to open file %s\n", fname);
  }

  // read size of each table of values
  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_nH);
  status = H5Dclose(dataset);

  dataset =
      H5Dopen(tempfile_id, "/Header/Number_of_helium_fractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_He);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_SolarAbundances);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_metals", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &cooling->N_Elements);
  status = H5Dclose(dataset);

  // allocate arrays of values for each of the above quantities
  cooling->Temp = malloc(cooling->N_Temp * sizeof(float));
  cooling->nH = malloc(cooling->N_nH * sizeof(float));
  cooling->HeFrac = malloc(cooling->N_He * sizeof(float));
  cooling->SolarAbundances = malloc(cooling->N_SolarAbundances * sizeof(float));
  cooling->Therm = malloc(cooling->N_Temp * sizeof(float));

  // read in values for each of the arrays
  dataset = H5Dopen(tempfile_id, "/Solar/Temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Temp);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Solar/Hydrogen_density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->nH);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Helium_mass_fraction_bins",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->HeFrac);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Solar_mass_fractions",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->SolarAbundances);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Metal_free/Temperature/Energy_density_bins",
                    H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Therm);
  status = H5Dclose(dataset);

  // Convert to temperature, density and internal energy arrays to log10
  for (i = 0; i < cooling->N_Temp; i++) {
    cooling->Temp[i] = log10(cooling->Temp[i]);
    cooling->Therm[i] = log10(cooling->Therm[i]);
  }

  for (i = 0; i < cooling->N_nH; i++) cooling->nH[i] = log10(cooling->nH[i]);
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
 * @brief Get the table of cooling rates for photoionized cooling (before
 * redshift ~9) Reads in table of cooling rates and electron abundances due to
 * metals (depending on temperature, hydrogen number density), cooling rates and
 * electron abundances due to hydrogen and helium (depending on temperature,
 * hydrogen number density and helium fraction), and temperatures (depending on
 * internal energy, hydrogen number density and helium fraction; note: this is
 * distinct from table of temperatures read in ReadCoolingHeader, as that table
 * is used to index the cooling, electron abundance tables, whereas this one is
 * used to obtain temperature of particle)
 *
 * @param cooling_table_path Directory containing cooling table
 * @param cooling Cooling data structure
 */

struct cooling_tables_redshift_invariant get_no_compt_table(
    char *cooling_table_path,
    const struct cooling_function_data *restrict cooling) {
#ifdef HAVE_HDF5

  struct cooling_tables_redshift_invariant cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  // allocate temporary arrays (needed to change order of dimensions
  // of arrays)
  net_cooling_rate =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  electron_abundance =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  temperature = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                cooling->N_nH * sizeof(float));
  he_net_cooling_rate = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                        cooling->N_nH * sizeof(float));
  he_electron_abundance = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                          cooling->N_nH * sizeof(float));

  // allocate arrays that the values will be assigned to
  cooling_table.metal_heating = (float *)malloc(
      cooling->N_Elements * cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.electron_abundance =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.temperature = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                              cooling->N_nH * sizeof(float));
  cooling_table.H_plus_He_heating = (float *)malloc(
      cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.H_plus_He_electron_abundance = (float *)malloc(
      cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float));

  sprintf(fname, "%sz_8.989nocompton.hdf5", cooling_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  // read in cooling rates due to metals
  for (specs = 0; specs < cooling->N_Elements; specs++) {
    sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     net_cooling_rate);
    status = H5Dclose(dataset);

    for (j = 0; j < cooling->N_Temp; j++) {
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_2d(j, k, cooling->N_Temp, cooling->N_nH);
        cooling_index = row_major_index_3d(
            k, j, specs, cooling->N_nH, cooling->N_Temp,
            cooling->N_Elements);  
        cooling_table.metal_heating[cooling_index] =
            -net_cooling_rate[table_index];
      }
    }
  }

  // read in cooling rates due to hydrogen and helium, H + He electron
  // abundances, temperatures
  strcpy(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_net_cooling_rate);
  status = H5Dclose(dataset);

  strcpy(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  status = H5Dclose(dataset);

  strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_He; i++) {
    for (j = 0; j < cooling->N_Temp; j++) {
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_3d(i, j, k, cooling->N_He,
                                         cooling->N_Temp, cooling->N_nH);
        cooling_index =
            row_major_index_3d(k, i, j, cooling->N_nH, cooling->N_He,
                               cooling->N_Temp);  
        cooling_table.H_plus_He_heating[cooling_index] =
            -he_net_cooling_rate[table_index];
        cooling_table.H_plus_He_electron_abundance[cooling_index] =
            he_electron_abundance[table_index];
        cooling_table.temperature[cooling_index] =
            log10(temperature[table_index]);
      }
    }
  }

  // read in electron densities due to metals
  strcpy(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_Temp; i++) {
    for (j = 0; j < cooling->N_nH; j++) {
      table_index = row_major_index_2d(i, j, cooling->N_Temp, cooling->N_nH);
      cooling_index =
          row_major_index_2d(j, i, cooling->N_nH, cooling->N_Temp);
      cooling_table.electron_abundance[cooling_index] =
          electron_abundance[table_index];
    }
  }

  status = H5Fclose(file_id);

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

#ifdef SWIFT_DEBUG_CHECKS
  message("done reading in no compton table");
#endif

  return cooling_table;
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
 * @brief Get the table of cooling rates for collisional cooling
 * Reads in table of cooling rates and electron abundances due to metals
 * (depending on temperature, hydrogen number density), cooling rates and
 * electron abundances due to hydrogen and helium (depending on temperature,
 * hydrogen number density and helium fraction), and temperatures (depending on
 * internal energy, hydrogen number density and helium fraction; note: this is
 * distinct from table of temperatures read in ReadCoolingHeader, as that table
 * is used to index the cooling, electron abundance tables, whereas this one is
 * used to obtain temperature of particle)
 *
 * @param cooling_table_path Directory containing cooling table
 * @param cooling Cooling data structure
 */
struct cooling_tables_redshift_invariant get_collisional_table(
    char *cooling_table_path,
    const struct cooling_function_data *restrict cooling) {
#ifdef HAVE_HDF5

  struct cooling_tables_redshift_invariant cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  // allocate temporary arrays (needed to change order of dimensions
  // of arrays)
  net_cooling_rate = (float *)malloc(cooling->N_Temp * sizeof(float));
  electron_abundance = (float *)malloc(cooling->N_Temp * sizeof(float));
  temperature =
      (float *)malloc(cooling->N_He * cooling->N_Temp * sizeof(float));
  he_net_cooling_rate =
      (float *)malloc(cooling->N_He * cooling->N_Temp * sizeof(float));
  he_electron_abundance =
      (float *)malloc(cooling->N_He * cooling->N_Temp * sizeof(float));

  // allocate arrays that the values will be assigned to
  cooling_table.metal_heating =
      (float *)malloc(cooling->N_Elements * cooling->N_Temp * sizeof(float));
  cooling_table.electron_abundance =
      (float *)malloc(cooling->N_Temp * sizeof(float));
  cooling_table.temperature =
      (float *)malloc(cooling->N_He * cooling->N_Temp * sizeof(float));
  cooling_table.H_plus_He_heating =
      (float *)malloc(cooling->N_He * cooling->N_Temp * sizeof(float));
  cooling_table.H_plus_He_electron_abundance =
      (float *)malloc(cooling->N_He * cooling->N_Temp * sizeof(float));

  sprintf(fname, "%sz_collis.hdf5", cooling_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  // read in cooling rates due to metals
  for (specs = 0; specs < cooling->N_Elements; specs++) {
    sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     net_cooling_rate);
    status = H5Dclose(dataset);

    for (j = 0; j < cooling->N_Temp; j++) {
      cooling_index =
          row_major_index_2d(specs, j, cooling->N_Elements,
                             cooling->N_Temp);  
      cooling_table.metal_heating[cooling_index] = -net_cooling_rate[j];
    }
  }

  // read in cooling rates due to hydrogen and helium, H + He electron
  // abundances, temperatures
  strcpy(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_net_cooling_rate);
  status = H5Dclose(dataset);

  strcpy(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  status = H5Dclose(dataset);

  strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_He; i++) {
    for (j = 0; j < cooling->N_Temp; j++) {
      table_index = row_major_index_2d(i, j, cooling->N_He, cooling->N_Temp);
      cooling_index =
          row_major_index_2d(i, j, cooling->N_He,
                             cooling->N_Temp);  
      cooling_table.H_plus_He_heating[cooling_index] =
          -he_net_cooling_rate[table_index];
      cooling_table.H_plus_He_electron_abundance[cooling_index] =
          he_electron_abundance[table_index];
      cooling_table.temperature[cooling_index] =
          log10(temperature[table_index]);
    }
  }

  // read in electron densities due to metals
  strcpy(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_Temp; i++) {
    cooling_table.electron_abundance[i] = electron_abundance[i];
  }

  status = H5Fclose(file_id);

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

#ifdef SWIFT_DEBUG_CHECKS
  message("done reading in collisional table");
#endif

  return cooling_table;
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
 * @brief Get the table of cooling rates for photodissociation cooling
 * Reads in table of cooling rates and electron abundances due to metals
 * (depending on temperature, hydrogen number density), cooling rates and
 * electron abundances due to hydrogen and helium (depending on temperature,
 * hydrogen number density and helium fraction), and temperatures (depending on
 * internal energy, hydrogen number density and helium fraction; note: this is
 * distinct from table of temperatures read in ReadCoolingHeader, as that table
 * is used to index the cooling, electron abundance tables, whereas this one is
 * used to obtain temperature of particle)
 *
 * @param cooling_table_path Directory containing cooling table
 * @param cooling Cooling data structure
 */
struct cooling_tables_redshift_invariant get_photodis_table(
    char *cooling_table_path,
    const struct cooling_function_data *restrict cooling) {
#ifdef HAVE_HDF5

  struct cooling_tables_redshift_invariant cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  // allocate temporary arrays (needed to change order of dimensions
  // of arrays)
  net_cooling_rate =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  electron_abundance =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  temperature = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                cooling->N_nH * sizeof(float));
  he_net_cooling_rate = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                        cooling->N_nH * sizeof(float));
  he_electron_abundance = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                          cooling->N_nH * sizeof(float));

  // allocate arrays that the values will be assigned to
  cooling_table.metal_heating = (float *)malloc(
      cooling->N_Elements * cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.electron_abundance =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.temperature = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                              cooling->N_nH * sizeof(float));
  cooling_table.H_plus_He_heating = (float *)malloc(
      cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.H_plus_He_electron_abundance = (float *)malloc(
      cooling->N_He * cooling->N_Temp * cooling->N_nH * sizeof(float));

  sprintf(fname, "%sz_photodis.hdf5", cooling_table_path);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  // read in cooling rates due to metals
  for (specs = 0; specs < cooling->N_Elements; specs++) {
    sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     net_cooling_rate);
    status = H5Dclose(dataset);

    for (j = 0; j < cooling->N_Temp; j++) {
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_2d(j, k, cooling->N_Temp, cooling->N_nH);
        cooling_index = row_major_index_3d(
            k, j, specs, cooling->N_nH, cooling->N_Temp,
            cooling->N_Elements);  
        cooling_table.metal_heating[cooling_index] =
            -net_cooling_rate[table_index];
      }
    }
  }

  // read in cooling rates due to hydrogen and helium, H + He electron
  // abundances, temperatures
  strcpy(set_name, "/Metal_free/Net_Cooling");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_net_cooling_rate);
  status = H5Dclose(dataset);

  strcpy(set_name, "/Metal_free/Temperature/Temperature");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temperature);
  status = H5Dclose(dataset);

  strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   he_electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_He; i++) {
    for (j = 0; j < cooling->N_Temp; j++) {
      for (k = 0; k < cooling->N_nH; k++) {
        table_index = row_major_index_3d(i, j, k, cooling->N_He,
                                         cooling->N_Temp, cooling->N_nH);
        cooling_index =
            row_major_index_3d(k, i, j, cooling->N_nH, cooling->N_He,
                               cooling->N_Temp);  
        cooling_table.H_plus_He_heating[cooling_index] =
            -he_net_cooling_rate[table_index];
        cooling_table.H_plus_He_electron_abundance[cooling_index] =
            he_electron_abundance[table_index];
        cooling_table.temperature[cooling_index] =
            log10(temperature[table_index]);
      }
    }
  }

  // read in electron densities due to metals
  strcpy(set_name, "/Solar/Electron_density_over_n_h");
  dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   electron_abundance);
  status = H5Dclose(dataset);

  for (i = 0; i < cooling->N_Temp; i++) {
    for (j = 0; j < cooling->N_nH; j++) {
      table_index = row_major_index_2d(i, j, cooling->N_Temp, cooling->N_nH);
      cooling_index = row_major_index_2d(
          j, i, cooling->N_nH, cooling->N_Temp);  
      cooling_table.electron_abundance[cooling_index] =
          electron_abundance[table_index];
    }
  }

  status = H5Fclose(file_id);

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

#ifdef SWIFT_DEBUG_CHECKS
  message("done reading in photodissociation table");
#endif

  return cooling_table;
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
 * @brief Get the cooling tables dependent on redshift
 *
 * @param cooling_table_path Directory containing cooling table
 * @param cooling Cooling data structure
 */

struct cooling_tables get_cooling_table(
    char *cooling_table_path,
    const struct cooling_function_data *restrict cooling) {
#ifdef HAVE_HDF5

  struct cooling_tables cooling_table;
  hid_t file_id, dataset;

  herr_t status;

  char fname[500], set_name[500];

  int specs, i, j, k, table_index, cooling_index;

  float *net_cooling_rate;
  float *electron_abundance;
  float *temperature;
  float *he_net_cooling_rate;
  float *he_electron_abundance;

  // allocate temporary arrays (needed to change order of dimensions
  // of arrays)
  net_cooling_rate =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  electron_abundance =
      (float *)malloc(cooling->N_Temp * cooling->N_nH * sizeof(float));
  temperature = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                cooling->N_nH * sizeof(float));
  he_net_cooling_rate = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                        cooling->N_nH * sizeof(float));
  he_electron_abundance = (float *)malloc(cooling->N_He * cooling->N_Temp *
                                          cooling->N_nH * sizeof(float));

  // allocate arrays that the values will be assigned to
  cooling_table.metal_heating =
      (float *)malloc(cooling->N_Redshifts * cooling->N_Elements *
                      cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.electron_abundance = (float *)malloc(
      cooling->N_Redshifts * cooling->N_Temp * cooling->N_nH * sizeof(float));
  cooling_table.temperature =
      (float *)malloc(cooling->N_Redshifts * cooling->N_He * cooling->N_Temp *
                      cooling->N_nH * sizeof(float));
  cooling_table.H_plus_He_heating =
      (float *)malloc(cooling->N_Redshifts * cooling->N_He * cooling->N_Temp *
                      cooling->N_nH * sizeof(float));
  cooling_table.H_plus_He_electron_abundance =
      (float *)malloc(cooling->N_Redshifts * cooling->N_He * cooling->N_Temp *
                      cooling->N_nH * sizeof(float));

  // repeat for each redshift
  for (int z_index = 0; z_index < cooling->N_Redshifts; z_index++) {
    sprintf(fname, "%sz_%1.3f.hdf5", cooling_table_path,
            cooling->Redshifts[z_index]);
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0) {
      error("unable to open file %s", fname);
    }

    // read in cooling rates due to metals
    for (specs = 0; specs < cooling->N_Elements; specs++) {
      sprintf(set_name, "/%s/Net_Cooling", eagle_tables_element_names[specs]);
      dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       net_cooling_rate);
      status = H5Dclose(dataset);

      for (i = 0; i < cooling->N_nH; i++) {
        for (j = 0; j < cooling->N_Temp; j++) {
          table_index =
              row_major_index_2d(j, i, cooling->N_Temp, cooling->N_nH);
          cooling_index = row_major_index_4d(
              z_index, i, j, specs, cooling->N_Redshifts, cooling->N_nH,
              cooling->N_Temp, cooling->N_Elements);
          cooling_table.metal_heating[cooling_index] =
              -net_cooling_rate[table_index];
        }
      }
    }

    // read in cooling rates due to hydrogen and helium, H + He electron
    // abundances, temperatures
    strcpy(set_name, "/Metal_free/Net_Cooling");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_net_cooling_rate);
    status = H5Dclose(dataset);

    strcpy(set_name, "/Metal_free/Temperature/Temperature");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     temperature);
    status = H5Dclose(dataset);

    strcpy(set_name, "/Metal_free/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     he_electron_abundance);
    status = H5Dclose(dataset);

    for (i = 0; i < cooling->N_He; i++) {
      for (j = 0; j < cooling->N_Temp; j++) {
        for (k = 0; k < cooling->N_nH; k++) {
          table_index = row_major_index_3d(i, j, k, cooling->N_He,
                                           cooling->N_Temp, cooling->N_nH);
          cooling_index =
              row_major_index_4d(z_index, k, i, j, cooling->N_Redshifts,
                                 cooling->N_nH, cooling->N_He, cooling->N_Temp);
          cooling_table.H_plus_He_heating[cooling_index] =
              -he_net_cooling_rate[table_index];
          cooling_table.H_plus_He_electron_abundance[cooling_index] =
              he_electron_abundance[table_index];
          cooling_table.temperature[cooling_index] =
              log10(temperature[table_index]);
        }
      }
    }

    // read in electron densities due to metals
    strcpy(set_name, "/Solar/Electron_density_over_n_h");
    dataset = H5Dopen(file_id, set_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     electron_abundance);
    status = H5Dclose(dataset);

    for (i = 0; i < cooling->N_Temp; i++) {
      for (j = 0; j < cooling->N_nH; j++) {
        table_index = row_major_index_2d(i, j, cooling->N_Temp, cooling->N_nH);
        cooling_index = row_major_index_3d(z_index, j, i, cooling->N_Redshifts,
                                           cooling->N_nH, cooling->N_Temp);
        cooling_table.electron_abundance[cooling_index] =
            electron_abundance[table_index];
      }
    }

    status = H5Fclose(file_id);
  }

  free(net_cooling_rate);
  free(electron_abundance);
  free(temperature);
  free(he_net_cooling_rate);
  free(he_electron_abundance);

#ifdef SWIFT_DEBUG_CHECKS
  message("done reading in general cooling table");
#endif

  return cooling_table;
#else
  error("Need HDF5 to read cooling tables");
#endif
}

/*
 * @brief Constructs the data structure containting all the cooling tables
 *
 * @param cooling_table_path Directory containing cooling table
 * @param cooling Cooling data structure
 */
struct eagle_cooling_table eagle_readtable(
    char *cooling_table_path,
    const struct cooling_function_data *restrict cooling) {

  struct eagle_cooling_table table;

  table.no_compton_cooling = get_no_compt_table(cooling_table_path, cooling);
  table.photodissociation_cooling =
      get_photodis_table(cooling_table_path, cooling);
  table.collisional_cooling =
      get_collisional_table(cooling_table_path, cooling);
  table.element_cooling = get_cooling_table(cooling_table_path, cooling);

  return table;
}

#endif
