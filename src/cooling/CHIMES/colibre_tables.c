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

/**
 * @file src/cooling/CHIMES/colibre_tables.c
 * @brief Functions to read COLIBRE tables
 */

/* Config parameters. */
#include "../config.h"

/* This file's header */
#include "colibre_tables.h"

/* Standard includes */
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "colibre_tables.h"
#include "error.h"
#include "exp10.h"

/**
 * @brief Reads in COLIBRE cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, metallicity,
 * abundance ratios, and elements used to index the cooling tables.
 *
 * @param fname Filepath for cooling table from which to read header
 * @param cooling Cooling data structure
 */
void read_cooling_header(struct colibre_cooling_tables *table) {

#ifdef HAVE_HDF5

  hid_t dataset;
  herr_t status;

  /* read sizes of array dimensions */
  hid_t tempfile_id =
      H5Fopen(table->cooling_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0)
    error("unable to open file %s\n", table->cooling_table_path);

  /* allocate arrays of bins */
  if (posix_memalign((void **)&table->Temp, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_temperature * sizeof(float)) != 0)
    error("Failed to allocate temperature table\n");

  if (posix_memalign((void **)&table->Redshifts, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_redshifts * sizeof(float)) != 0)
    error("Failed to allocate redshift table\n");

  if (posix_memalign((void **)&table->nH, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate density table\n");

  if (posix_memalign((void **)&table->Metallicity, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity * sizeof(float)) != 0)
    error("Failed to allocate metallicity table\n");

  if (posix_memalign((void **)&table->LogAbundances, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity *
                         colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&table->Abundances, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity *
                         colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&table->Abundances_inv, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity *
                         colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&table->atomicmass, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate atomic masses array\n");

  if (posix_memalign((void **)&table->atomicmass_inv, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate inverse atomic masses array\n");

  if (posix_memalign((void **)&table->Zsol, SWIFT_STRUCT_ALIGNMENT,
                     1 * sizeof(float)) != 0)
    error("Failed to allocate solar metallicity array\n");

  if (posix_memalign((void **)&table->Zsol_inv, SWIFT_STRUCT_ALIGNMENT,
                     1 * sizeof(float)) != 0)
    error("Failed to allocate inverse solar metallicity array\n");

  if (posix_memalign((void **)&table->LogMassFractions,
                     SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity *
                         colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate log mass fraction array\n");

  if (posix_memalign((void **)&table->MassFractions, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_metallicity *
                         colibre_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate mass fraction array\n");

  /* read in bins and misc information */
  dataset = H5Dopen(tempfile_id, "/TableBins/TemperatureBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Temp);
  if (status < 0) error("error reading temperature bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/RedshiftBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Redshifts);
  if (status < 0) error("error reading redshift bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/DensityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->nH);
  if (status < 0) error("error reading density bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/MetallicityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Metallicity);
  if (status < 0) error("error reading metallicity bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TotalAbundances", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->LogAbundances);
  if (status < 0) error("error reading total abundances\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TotalMassFractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->LogMassFractions);
  if (status < 0) error("error reading total mass fractions\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/ElementMasses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->atomicmass);
  if (status < 0) error("error reading element masses\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/SolarMetallicity", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Zsol);
  if (status < 0) error("error reading solar metallicity \n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  table->Zsol_inv[0] = 1.f / table->Zsol[0];

  /* find the metallicity bin that refers to solar metallicity */
  const float tol = 1.e-3;
  for (int i = 0; i < colibre_cooling_N_metallicity; i++) {
    if (fabsf(table->Metallicity[i]) < tol) {
      table->indxZsol = i;
    }
  }

#if defined(__ICC)
#pragma novector
#endif
  for (int i = 0; i < colibre_cooling_N_elementtypes; i++) {
    table->atomicmass_inv[i] = 1.f / table->atomicmass[i];
  }

  /* set some additional useful abundance arrays */
  for (int i = 0; i < colibre_cooling_N_metallicity; i++) {

#if defined(__ICC)
#pragma novector
#endif
    for (int j = 0; j < colibre_cooling_N_elementtypes; j++) {
      const int indx1d = cooling_row_major_index_2d(i, j, colibre_cooling_N_metallicity,
                                            colibre_cooling_N_elementtypes);
      table->Abundances[indx1d] = exp10f(table->LogAbundances[indx1d]);
      table->Abundances_inv[indx1d] = 1.f / table->Abundances[indx1d];
      table->MassFractions[indx1d] =
          exp10f(table->LogMassFractions[indx1d]);
    }
  }

#else
  error("Need HDF5 to read cooling tables");
#endif
}

/**
 * @brief Allocate space for cooling tables and read them
 *
 * @param cooling #cooling_function_data structure
 */
void read_cooling_tables(struct colibre_cooling_tables *restrict table) {

#ifdef HAVE_HDF5
  hid_t dataset;
  herr_t status;

  /* open hdf5 file */
  hid_t tempfile_id =
      H5Fopen(table->cooling_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0)
    error("unable to open file %s\n", table->cooling_table_path);

  /* Allocate and read arrays to store cooling tables. */

  /* Cooling (temperature) */
  if (posix_memalign(
          (void **)&table->Tcooling, SWIFT_STRUCT_ALIGNMENT,
          colibre_cooling_N_redshifts * colibre_cooling_N_temperature *
              colibre_cooling_N_metallicity * colibre_cooling_N_density *
              colibre_cooling_N_cooltypes * sizeof(float)) != 0)
    error("Failed to allocate Tcooling array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/Cooling", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Tcooling);
  if (status < 0) error("error reading Tcooling\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Heating (temperature) */
  if (posix_memalign(
          (void **)&table->Theating, SWIFT_STRUCT_ALIGNMENT,
          colibre_cooling_N_redshifts * colibre_cooling_N_temperature *
              colibre_cooling_N_metallicity * colibre_cooling_N_density *
              colibre_cooling_N_heattypes * sizeof(float)) != 0)
    error("Failed to allocate Theating array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/Heating", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Theating);
  if (status < 0) error("error reading Theating\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Electron fraction (temperature) */
  if (posix_memalign(
          (void **)&table->Telectron_fraction, SWIFT_STRUCT_ALIGNMENT,
          colibre_cooling_N_redshifts * colibre_cooling_N_temperature *
              colibre_cooling_N_metallicity * colibre_cooling_N_density *
              colibre_cooling_N_electrontypes * sizeof(float)) != 0)
    error("Failed to allocate Telectron_fraction array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/ElectronFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->Telectron_fraction);
  if (status < 0) error("error reading electron_fraction (temperature)\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

#ifdef SWIFT_DEBUG_CHECKS
  message("Done reading in general cooling table");
#endif

#else
  error("Need HDF5 to read cooling tables");
#endif
}
