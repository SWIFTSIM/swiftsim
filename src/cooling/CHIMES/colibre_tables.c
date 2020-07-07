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

/* Standard includes */
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "error.h"
#include "exp10.h"

/**
 * @brief Reads in COLIBRE cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, metallicity,
 * abundance ratios, and elements used to index the cooling tables.
 *
 * @param table Colibre cooling table structure.
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

  if (posix_memalign((void **)&table->LogMassFractions, SWIFT_STRUCT_ALIGNMENT,
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

  /* Close the file */
  H5Fclose(tempfile_id);

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
      const int indx1d = cooling_row_major_index_2d(
          i, j, colibre_cooling_N_metallicity, colibre_cooling_N_elementtypes);
      table->Abundances[indx1d] = exp10f(table->LogAbundances[indx1d]);
      table->Abundances_inv[indx1d] = 1.f / table->Abundances[indx1d];
      table->MassFractions[indx1d] = exp10f(table->LogMassFractions[indx1d]);
    }
  }

#else
  error("Need HDF5 to read cooling tables");
#endif
}

/**
 * @brief Allocate space for cooling tables and read them
 *
 * @param table Colibre cooling table structure.
 */
void read_cooling_tables(struct colibre_cooling_tables *table /*, struct dustevo_props *dp*/) {

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

  /* Thermal equilibrium temperature */
  if (posix_memalign((void **)&table->logTeq, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_redshifts *
                         colibre_cooling_N_metallicity *
                         colibre_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate logTeq array\n");

  dataset = H5Dopen(tempfile_id, "/ThermEq/Temperature", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->logTeq);
  if (status < 0) error("error reading Teq array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing logTeq dataset");

  /* Mean particle mass at thermal equilibrium temperature */
  if (posix_memalign((void **)&table->meanpartmass_Teq, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_redshifts *
                         colibre_cooling_N_metallicity *
                         colibre_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate mu array\n");

  dataset = H5Dopen(tempfile_id, "/ThermEq/MeanParticleMass", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   table->meanpartmass_Teq);
  if (status < 0) error("error reading mu array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing mu dataset");

  /** If I could #include "dust.h" in colibre_tables.h, would run 
   * read_colibre_depletion() from src/dust/T20/dust_yield_tables.h 
   * here **/

  /* Close the file */
  H5Fclose(tempfile_id);

  /* Pressure at thermal equilibrium temperature */
  if (posix_memalign((void **)&table->logPeq, SWIFT_STRUCT_ALIGNMENT,
                     colibre_cooling_N_redshifts *
                         colibre_cooling_N_metallicity *
                         colibre_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate logPeq array\n");

  const double log10_kB_cgs = table->log10_kB_cgs;

  /* Compute the pressures at thermal eq. */
  for (int ired = 0; ired < colibre_cooling_N_redshifts; ired++) {
    for (int imet = 0; imet < colibre_cooling_N_metallicity; imet++) {

      const int index_XH =
          cooling_row_major_index_2d(imet, 0, colibre_cooling_N_metallicity,
                                     colibre_cooling_N_elementtypes);

      const float log10_XH = table->LogMassFractions[index_XH];

      for (int iden = 0; iden < colibre_cooling_N_density; iden++) {

        const int index_Peq = cooling_row_major_index_3d(
            ired, imet, iden, colibre_cooling_N_redshifts,
            colibre_cooling_N_metallicity, colibre_cooling_N_density);

        table->logPeq[index_Peq] =
            table->nH[iden] + table->logTeq[index_Peq] - log10_XH -
            log10(table->meanpartmass_Teq[index_Peq]) + log10_kB_cgs;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  message("Done reading in general cooling table");
#endif

#else
  error("Need HDF5 to read cooling tables");
#endif
}

/**
 * @brief Computes net metal heating rate from Colibre tables.
 *
 * Computes the net heating rate (heating - cooling) for a given element
 * abundance ratio, temperature, redshift, and density. The unit of the net
 * cooling rate is Lambda / nH**2 [erg cm^3 s-1] and all input values are in
 * cgs.
 *
 * @param myGasVars The #gasVariables struct.
 * @param myGlobalVars The #globalVariables struct.
 */
double colibre_metal_cooling_rate_temperature(
    struct gasVariables *myGasVars, struct globalVariables *myGlobalVars) {
  struct global_hybrid_data_struct *myGlobalData;
  myGlobalData = (struct global_hybrid_data_struct *)myGlobalVars->hybrid_data;
  struct colibre_cooling_tables *table = myGlobalData->table;

  double log_T_cgs = log10(myGasVars->temperature);
  double redshift = myGlobalVars->redshift;
  double n_H_cgs = myGasVars->nH_tot;
  double Z_absolute = myGasVars->metallicity * myGlobalData->Zsol;
  double noneq_electron_fraction = myGasVars->abundances[sp_elec];

  const float *abundance_ratio;
  struct gas_hybrid_data_struct *myGasData;
  myGasData = (struct gas_hybrid_data_struct *)myGasVars->hybrid_data;
  abundance_ratio = myGasData->abundance_ratio;

  /* Set weights for cooling rates */
  float weights_cooling[colibre_cooling_N_cooltypes - 2];
  for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {

    if (i <= element_He) {
      /* H and He are in CHIMES. */
      weights_cooling[i] = 0.f;
    } else if (i < element_OA) {
      /* Only include metals that aren't
       * already included in CHIMES. */
      if (myGlobalVars->element_included[i - 2] == 0)
        weights_cooling[i] = abundance_ratio[i];
      else
        weights_cooling[i] = 0.f;
    } else if (i == element_OA) {
      weights_cooling[i] = abundance_ratio[i];
    } else if ((i == cooltype_H2) || (i == cooltype_NetFFH) ||
               (i == cooltype_Compton) || (i == cooltype_Dust)) {
      /* These channels are in CHIMES. */
      weights_cooling[i] = 0.f;
    } else {
      /* use same abundances as in the tables */
      weights_cooling[i] = 1.f;
    }
  }

  /* Set weights for heating rates */
  float weights_heating[colibre_cooling_N_heattypes - 2];
  for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
    if (i <= element_He) {
      /* H and He are in CHIMES. */
      weights_heating[i] = 0.f;
    } else if (i < element_OA) {
      /* Only include metals that aren't
       * already included in CHIMES. */
      if (myGlobalVars->element_included[i - 2] == 0)
        weights_heating[i] = abundance_ratio[i];
      else
        weights_heating[i] = 0.f;
    } else if (i == element_OA) {
      weights_heating[i] = abundance_ratio[i];
    } else if ((i == heattype_H2) || (i == heattype_CosmicRay) ||
               (i == heattype_HFF) || (i == heattype_Compton) ||
               (i == heattype_Dust)) {
      /* These channels are in CHIMES. */
      weights_heating[i] = 0.f;
    } else {
      weights_heating[i] = 1.f; /* use same abundances as in the tables */
    }
  }

  // Set weights for electron densities.
  float weights_electron[colibre_cooling_N_electrontypes];
  for (int i = 0; i < colibre_cooling_N_electrontypes; i++) {
    if (i < colibre_cooling_N_elementtypes - 1)
      weights_electron[i] = abundance_ratio[i];
    else
      weights_electron[i] = 1.f;
  }

  /* Get indices of T, nH, metallicity and redshift */
  int T_index, n_H_index, met_index, red_index;
  float d_T, d_n_H, d_met, d_red;
  float logZZsol = log10((Z_absolute / table->Zsol[0]) + FLT_MIN);

  cooling_get_index_1d(table->Redshifts, colibre_cooling_N_redshifts, redshift,
                       &red_index, &d_red);
  cooling_get_index_1d(table->Temp, colibre_cooling_N_temperature, log_T_cgs,
                       &T_index, &d_T);
  cooling_get_index_1d(table->Metallicity, colibre_cooling_N_metallicity,
                       logZZsol, &met_index, &d_met);
  cooling_get_index_1d(table->nH, colibre_cooling_N_density, log10(n_H_cgs),
                       &n_H_index, &d_n_H);

  // n_e / n_H from Colibre

  // From H + He
  double colibre_electron_fraction_prim = interpolation4d_plus_summation(
      table->Telectron_fraction, weights_electron,
      colibre_cooling_N_electrontypes - 3, colibre_cooling_N_electrontypes - 3,
      red_index, T_index, met_index, n_H_index, d_red, d_T, d_met, d_n_H,
      colibre_cooling_N_redshifts, colibre_cooling_N_temperature,
      colibre_cooling_N_metallicity, colibre_cooling_N_density,
      colibre_cooling_N_electrontypes);

  // From metals, with table metal ratios
  double colibre_electron_fraction_metal_table = interpolation4d_plus_summation(
      table->Telectron_fraction, weights_electron,
      colibre_cooling_N_electrontypes - 2, colibre_cooling_N_electrontypes - 2,
      red_index, T_index, met_index, n_H_index, d_red, d_T, d_met, d_n_H,
      colibre_cooling_N_redshifts, colibre_cooling_N_temperature,
      colibre_cooling_N_metallicity, colibre_cooling_N_density,
      colibre_cooling_N_electrontypes);

  /* From metals, with actual metal ratios.
   * We also need to exclude any metals
   * that are already included in CHIMES. */
  for (int i = element_C; i < element_OA; i++) {
    if (myGlobalVars->element_included[i - 2] == 1) weights_electron[i] = 0.f;
  }

  double colibre_electron_fraction_metal_actual =
      interpolation4d_plus_summation(
          table->Telectron_fraction, weights_electron, element_C,
          colibre_cooling_N_electrontypes - 4, red_index, T_index, met_index,
          n_H_index, d_red, d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
          colibre_cooling_N_temperature, colibre_cooling_N_metallicity,
          colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  /* Adjust the weights based on the non-eq
   * electron fraction from CHIMES, compared
   * to the electron fraction that was used
   * in the Colibre tables. */
  double electron_fraction_ratio =
      (noneq_electron_fraction + colibre_electron_fraction_metal_actual) /
      (colibre_electron_fraction_prim + colibre_electron_fraction_metal_table +
       FLT_MIN);

  for (int i = element_C; i < colibre_cooling_N_elementtypes; i++)
    weights_cooling[i] *= electron_fraction_ratio;

  weights_cooling[cooltype_molecules] *= electron_fraction_ratio;
  weights_cooling[cooltype_HD] *= electron_fraction_ratio;
  weights_cooling[cooltype_NetFFM] *= electron_fraction_ratio;
  weights_cooling[cooltype_eeBrems] *=
      electron_fraction_ratio * electron_fraction_ratio;

  /* Lambda / n_H**2 */
  const double cooling_rate = interpolation4d_plus_summation(
      table->Tcooling, weights_cooling,           /* */
      element_C, colibre_cooling_N_cooltypes - 3, /* */
      red_index, T_index, met_index, n_H_index,   /* */
      d_red, d_T, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_temperature,              /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_cooltypes);               /* */

  /* Gamma / n_H**2 */
  const double heating_rate = interpolation4d_plus_summation(
      table->Theating, weights_heating,           /* */
      element_C, colibre_cooling_N_heattypes - 3, /* */
      red_index, T_index, met_index, n_H_index,   /* */
      d_red, d_T, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_temperature,              /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_heattypes);               /* */

  /* Return the net heating rate (Lambda_heat - Lambda_cool) */
  return heating_rate - cooling_rate;
}

/**
 * @brief Allocate gas_hybrid_data struct in gasVars.
 *
 * @param myGasVars The #gasVariables struct.
 */
void chimes_allocate_gas_hybrid_data(struct gasVariables *myGasVars) {
  myGasVars->hybrid_data =
      (void *)malloc(sizeof(struct gas_hybrid_data_struct));
}

/**
 * @brief Allocate gas_hybrid_data struct in gasVars.
 *
 * @param myGasVars The #gasVariables struct.
 */
void chimes_free_gas_hybrid_data(struct gasVariables *myGasVars) {
  free(myGasVars->hybrid_data);
}
