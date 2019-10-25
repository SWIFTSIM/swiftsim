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
#ifndef SWIFT_COLIBRE_TABLES_H
#define SWIFT_COLIBRE_TABLES_H

/**
 * @file src/cooling/CHIMES/colibre_tables.h
 * @brief COLIBRE cooling tables
 */

/* Config parameters. */
#include "config.h"

/* Some standard headers. */ 
#include <float.h> 
#include <math.h> 

/* Local includes. */
#include "align.h"
#include "chemistry_struct.h"
#include "error.h"
#include "inline.h"

#define colibre_table_path_name_length 500

/**
 * @brief struct containing cooling tables
 */
struct colibre_cooling_tables {

  /*! Filepath to the directory containing the HDF5 cooling tables */
  char cooling_table_path[colibre_table_path_name_length];

  /* array of all cooling processes (temperature) */
  float *Tcooling;

  /* array of all heating processes (temperature) */
  float *Theating;

  /* array of all electron abundances (temperature) */
  float *Telectron_fraction;

  /*! Redshift bins */
  float *Redshifts;

  /*! Hydrogen number density bins */
  float *nH;

  /*! Temperature bins */
  float *Temp;

  /*! Metallicity bins */
  float *Metallicity;

  /*! Abundance ratios for each metallicity bin and for each included element */
  float *LogAbundances;
  float *Abundances;
  float *Abundances_inv;

  /*! Atomic masses for all included elements */
  float *atomicmass;
  float *atomicmass_inv;

  /*! Mass fractions of all included elements */
  float *LogMassFractions;
  float *MassFractions;

  /*! Index for solar metallicity in the metallicity dimension */
  int indxZsol;

  /*! Solar metallicity (metal mass fraction) */
  float *Zsol;

  /*! Inverse of solar metallicity (metal mass fraction) */
  float *Zsol_inv;

  /*! Ca over Si abundance divided by the solar ratio for these elements */
  float Ca_over_Si_ratio_in_solar;

  /*! S over Si abundance divided by the solar ratio for these elements */
  float S_over_Si_ratio_in_solar;

  /*! Redshift of H reionization */
  float H_reion_z;
};

/*! Number of different bins along the temperature axis of the tables */
#define colibre_cooling_N_temperature 86

/*! Number of different bins along the redshift axis of the tables */
#define colibre_cooling_N_redshifts 10

/*! Number of different bins along the density axis of the tables */
#define colibre_cooling_N_density 71

/*! Number of different bins along the metallicity axis of the tables */
#define colibre_cooling_N_metallicity 11

/*! Number of different bins along the internal energy axis of the tables */
#define colibre_cooling_N_internalenergy 191

/*! Number of different cooling channels in the tables */
#define colibre_cooling_N_cooltypes 22

/*! Number of different heating channels in the tables */
#define colibre_cooling_N_heattypes 24

/*! Number of different electron fractions (each element - other atoms
 *  + tot prim + tot metal + tot)  in the tables */
#define colibre_cooling_N_electrontypes 14

/*! Number of different elements in the tables */
#define colibre_cooling_N_elementtypes 12

/**
 * @brief Elements present in the tables
 */
enum colibre_cooling_element {
  element_H,
  element_He,
  element_C,
  element_N,
  element_O,
  element_Ne,
  element_Mg,
  element_Si,
  element_S,
  element_Ca,
  element_Fe,
  element_OA
};

/**
 * @brief Hydrogen species
 */
enum colibre_hydrogen_species { neutral, ionized, molecular };

/**
 * @brief Cooling channels beyond the metal lines
 */
enum colibre_cooling_channels {
  cooltype_H2 = element_OA + 1,
  cooltype_molecules,
  cooltype_HD,
  cooltype_NetFFH,
  cooltype_NetFFM,
  cooltype_eeBrems,
  cooltype_Compton,
  cooltype_Dust
};

/**
 * @brief Heating channels beyond the metal lines
 */
enum colibre_heating_channels {
  heattype_H2 = element_OA + 1,
  heattype_COdiss,
  heattype_CosmicRay,
  heattype_UTA,
  heattype_line,
  heattype_Hlin,
  heattype_ChaT,
  heattype_HFF,
  heattype_Compton,
  heattype_Dust
};

void read_cooling_header(struct colibre_cooling_tables *table);
void read_cooling_tables(struct colibre_cooling_tables *table);

/**
 * @brief Returns the 1d index of element with 2d indices x,y
 * from a flattened 2d array in row major order
 *
 * @param x, y Indices of element of interest
 * @param Nx, Ny Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int cooling_row_major_index_2d(const int x,
                                                             const int y,
                                                             const int Nx,
                                                             const int Ny) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
#endif
  return x * Ny + y;
}

/**
 * @brief Returns the 1d index of element with 5d indices x,y,z,w
 * from a flattened 5d array in row major order
 *
 * @param x, y, z, v, w Indices of element of interest
 * @param Nx, Ny, Nz, Nv, Nw Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int cooling_row_major_index_5d(
    const int x, const int y, const int z, const int w, const int v,
    const int Nx, const int Ny, const int Nz, const int Nw, const int Nv) {

#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
  assert(w < Nw);
  assert(v < Nv);
#endif

  return x * Ny * Nz * Nw * Nv + y * Nz * Nw * Nv + z * Nw * Nv + w * Nv + v;
}

/**
 * @brief Finds the index of a value in a table and compute delta to nearest
 * element.
 *
 * This function assumes the table is monotonically increasing with a constant
 * difference between adjacent values.
 *
 * The returned difference is expressed in units of the table separation. This
 * means dx = (x - table[i]) / (table[i+1] - table[i]). It is always between
 * 0 and 1.
 *
 * We use a small epsilon of 1e-4 to avoid out-of-range accesses due to
 * rounding errors.
 *
 * @param table The table to search in.
 * @param size The number of elements in the table.
 * @param x The value to search for.
 * @param i (return) The index in the table of the element.
 * @param *dx (return) The difference between x and table[i]
 */
__attribute__((always_inline)) INLINE void cooling_get_index_1d(
    const float *restrict table, const int size, const float x, int *i,
    float *restrict dx) {

  /* Small epsilon to avoid rounding issues leading to out-of-bound
   * access when using the indices later to read data from the tables. */
  const float epsilon = 1e-4f;

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  /* Distance between elements in the array */
  /* Do not use first or last entry, might be an extra bin with uneven spacing
   */
  const float delta = (size - 3) / (table[size - 2] - table[1]);

  /* Check for an extra entry at the beginning (e.g. metallicity) */
  int istart = 0;
  int iend = size - 1;

  if (fabsf(table[1] - table[0]) > delta + epsilon) {
    istart = 1;
  }
  if (fabsf(table[size - 1] - table[size - 2]) > delta + epsilon) {
    iend = size - 2;
  }

  /*extra array at the beginning */
  if (x < table[istart] + epsilon) {
    /* We are before the first element */
    *i = 0;
    *dx = 0.f;
  } else if (x < table[iend] - epsilon) {
    *i = (x - table[1]) * delta + 1;
    *dx = (x - table[*i]) * delta;
  } else {
    /* We are after the last element */
    *i = iend - 1;
    *dx = 1.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (*dx < -0.001f || *dx > 1.001f) error("Invalid distance found dx=%e", *dx);
#endif
}

/**
 * @brief Interpolates a 5 dimensional array in the first 4 dimensions and
 * adds the individual contributions from the 5th dimension according to their
 * weights
 *
 * @param table The table to interpolate
 * @param weights The weights for summing up the individual contributions
 * @param istart, iend Start and stop index for 5th dimension
 * @param xi, yi, zi, wi Indices of table element
 * @param dx, dy, dz, dw Distance between the point and the index in units of
 * the grid spacing.
 * @param Nx, Ny, Nz, Nw, Nv Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE double interpolation4d_plus_summation(
    const float *table, const float *weights, const int istart, const int iend,
    const int xi, const int yi, const int zi, const int wi, const float dx,
    const float dy, const float dz, const float dw, const int Nx, const int Ny,
    const int Nz, const int Nw, const int Nv) {

  const float tx = 1.f - dx;
  const float ty = 1.f - dy;
  const float tz = 1.f - dz;
  const float tw = 1.f - dw;

  float result;
  double result_global = 0.;

  for (int i = istart; i <= iend; i++) {

    /* Linear interpolation along each axis. We read the table 2^4=16 times */
    result = tx * ty * tz * tw *
             table[cooling_row_major_index_5d(xi + 0, yi + 0, zi + 0, wi + 0, i, Nx, Ny,
                                      Nz, Nw, Nv)];

    result += tx * ty * tz * dw *
              table[cooling_row_major_index_5d(xi + 0, yi + 0, zi + 0, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];

    result += tx * ty * dz * tw *
              table[cooling_row_major_index_5d(xi + 0, yi + 0, zi + 1, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += tx * dy * tz * tw *
              table[cooling_row_major_index_5d(xi + 0, yi + 1, zi + 0, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += dx * ty * tz * tw *
              table[cooling_row_major_index_5d(xi + 1, yi + 0, zi + 0, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];

    result += tx * ty * dz * dw *
              table[cooling_row_major_index_5d(xi + 0, yi + 0, zi + 1, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += tx * dy * tz * dw *
              table[cooling_row_major_index_5d(xi + 0, yi + 1, zi + 0, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += dx * ty * tz * dw *
              table[cooling_row_major_index_5d(xi + 1, yi + 0, zi + 0, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += tx * dy * dz * tw *
              table[cooling_row_major_index_5d(xi + 0, yi + 1, zi + 1, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += dx * ty * dz * tw *
              table[cooling_row_major_index_5d(xi + 1, yi + 0, zi + 1, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += dx * dy * tz * tw *
              table[cooling_row_major_index_5d(xi + 1, yi + 1, zi + 0, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];

    result += dx * dy * dz * tw *
              table[cooling_row_major_index_5d(xi + 1, yi + 1, zi + 1, wi + 0, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += dx * dy * tz * dw *
              table[cooling_row_major_index_5d(xi + 1, yi + 1, zi + 0, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += dx * ty * dz * dw *
              table[cooling_row_major_index_5d(xi + 1, yi + 0, zi + 1, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];
    result += tx * dy * dz * dw *
              table[cooling_row_major_index_5d(xi + 0, yi + 1, zi + 1, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];

    result += dx * dy * dz * dw *
              table[cooling_row_major_index_5d(xi + 1, yi + 1, zi + 1, wi + 1, i, Nx,
                                       Ny, Nz, Nw, Nv)];

    result_global += weights[i] * pow(10.0, result);
  }

  return result_global;
}

/**
 * @brief Returns the element index of the particle-carried chemistry field
 * corresponding to a given cooling element.
 */
__attribute__((always_inline)) INLINE int element_from_table_to_code(
    const enum colibre_cooling_element i) {

#ifdef SWIFT_DEBUG_CHECKS
  if ((i >= colibre_cooling_N_elementtypes) || (i < 0))
    error("Outside range of elements in cooling tables");
#endif

  switch (i) {
    case element_H:
      return chemistry_element_H;
    case element_He:
      return chemistry_element_He;
    case element_C:
      return chemistry_element_C;
    case element_N:
      return chemistry_element_N;
    case element_O:
      return chemistry_element_O;
    case element_Ne:
      return chemistry_element_Ne;
    case element_Mg:
      return chemistry_element_Mg;
    case element_Si:
      return chemistry_element_Si;
      /* S and Ca are not tracked individually; their abundance is
       * assumed to be the same as Si (with respect to solar) */
    case element_S:
      return chemistry_element_Si;
    case element_Ca:
      return chemistry_element_Si;
    case element_Fe:
      return chemistry_element_Fe;
      /* other elements, if used, scale with metallicity */
    case element_OA:
      return -1;
  }

  return -1;
}

/**
 * @brief Computes the net heating rate (heating - cooling) for a given element
 * abundance ratio, temperature, redshift, and density. The unit of the net
 * cooling rate is Lambda / nH**2 [erg cm^3 s-1] and all input values are in
 * cgs. 
 *
 * @param log_T_cgs Log base 10 of temperature in K
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param Z_absolute absolute metal mass fraction 
 * @param noneq_electron_fraction ne/nH from CHIMES 
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param cooling #colibre_cooling_function_data structure
 */
INLINE static double colibre_metal_cooling_rate_temperature(
    double log_T_cgs, double redshift, double n_H_cgs, 
    double Z_absolute, double noneq_electron_fraction, 
    const float abundance_ratio[colibre_cooling_N_elementtypes], 
    const struct colibre_cooling_tables *restrict table) {

  /* Set weights for cooling rates */
  float weights_cooling[colibre_cooling_N_cooltypes - 2];
  for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {

    if (i <= element_He) {
      /* H and He are in CHIMES. */ 
      weights_cooling[i] = 0.f; 
    } else if (i < colibre_cooling_N_elementtypes) {
      /* Use abundance ratios */
      weights_cooling[i] = abundance_ratio[i];
    } else if ((i == cooltype_H2) || 
	       (i == cooltype_NetFFH) || 
	       (i == cooltype_Compton) || 
	       (i == cooltype_Dust)) {
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
    } else if (i < colibre_cooling_N_elementtypes) {
      weights_heating[i] = abundance_ratio[i];
    } else if ((i == heattype_H2) || 
	       (i == heattype_CosmicRay) || 
	       (i == heattype_HFF) || 
	       (i == heattype_Compton) || 
	       (i == heattype_Dust)){ 
      /* These channels are in CHIMES. */ 
      weights_heating[i] = 0.f; 
    } else {
      weights_heating[i] = 1.f; /* use same abundances as in the tables */
    }
  }

  /* Get indices of T, nH, metallicity and redshift */
  int T_index, n_H_index, met_index, red_index;
  float d_T, d_n_H, d_met, d_red;
  float logZZsol = log10((Z_absolute / table->Zsol[0]) + FLT_MIN); 

  if (redshift < table->H_reion_z) { 
    cooling_get_index_1d(table->Redshifts, colibre_cooling_N_redshifts, redshift,
		 &red_index, &d_red);
  } else {
    red_index = colibre_cooling_N_redshifts - 2;
    d_red = 1.0;
  }

  cooling_get_index_1d(table->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);
  cooling_get_index_1d(table->Metallicity, colibre_cooling_N_metallicity, logZZsol,
	       &met_index, &d_met);
  cooling_get_index_1d(table->nH, colibre_cooling_N_density, log10(n_H_cgs),
	       &n_H_index, &d_n_H);

  /* n_e / n_H from metals */
  double electron_fraction = interpolation4d_plus_summation(
      table->Telectron_fraction, abundance_ratio, /* */
      element_H, colibre_cooling_N_electrontypes - 4,     /* */
      red_index, T_index, met_index, n_H_index,           /* */
      d_red, d_T, d_met, d_n_H,                           /* */
      colibre_cooling_N_redshifts,                        /* */
      colibre_cooling_N_internalenergy,                   /* */
      colibre_cooling_N_metallicity,                      /* */
      colibre_cooling_N_density,                          /* */
      colibre_cooling_N_electrontypes);                   /* */

  /* Add on non-eq electron 
   * fraction from CHIMES. */ 
  electron_fraction += noneq_electron_fraction; 

  /* Lambda / n_H**2 */
  const double cooling_rate = interpolation4d_plus_summation(
      table->Tcooling, weights_cooling,   /* */
      element_H, colibre_cooling_N_cooltypes - 3, /* */
      red_index, T_index, met_index, n_H_index,   /* */
      d_red, d_T, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_internalenergy,           /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_cooltypes);               /* */

  /* Gamma / n_H**2 */
  const double heating_rate = interpolation4d_plus_summation(
      table->Theating, weights_heating,   /* */
      element_H, colibre_cooling_N_heattypes - 3, /* */
      red_index, T_index, met_index, n_H_index,   /* */
      d_red, d_T, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_internalenergy,           /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_heattypes);               /* */

  /* Return the net heating rate (Lambda_heat - Lambda_cool) */
  return heating_rate - cooling_rate;
}

#endif
