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
#ifndef SWIFT_COOLING_EAGLE_H
#define SWIFT_COOLING_EAGLE_H

/**
 * @file src/cooling/none/cooling.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <hdf5.h>
#include <time.h>

/* Local includes. */
#include "cooling_struct.h"
#include "error.h"
#include "hydro.h"
#include "chemistry.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"
#include "eagle_cool_tables.h"

/* number of calls to eagle cooling rate */
extern int n_eagle_cooling_rate_calls_1;
extern int n_eagle_cooling_rate_calls_2;
extern int n_eagle_cooling_rate_calls_3;
extern int n_eagle_cooling_rate_calls_4;

static int get_redshift_index_first_call = 0;
static int get_redshift_index_previous = -1;

enum hdf5_allowed_types {
  hdf5_short,
  hdf5_int,
  hdf5_long,
  hdf5_float,
  hdf5_double,
  hdf5_char
};

__attribute__((always_inline)) INLINE int row_major_index_2d(int i, int j,
                                                               int nx, int ny) {
  int index = i * ny + j;
#ifdef SWIFT_DEBUG_CHECKS
  assert(i < nx);
  assert(j < ny);
#endif
  return index;
}

__attribute__((always_inline)) INLINE int row_major_index_3d(int i, int j,
                                                              int k, int nx, 
							      int ny, int nz) {
  int index = i * ny * nz + j * nz + k;
#ifdef SWIFT_DEBUG_CHECKS
  assert(i < nx);
  assert(j < ny);
  assert(k < nz);
#endif
  return index;
}

__attribute__((always_inline)) INLINE int row_major_index_4d(int i, int j,
                                                              int k, int l, 
							      int nx, int ny, 
							      int nz, int nw) {
  int index = i * ny * nz * nw + j * nz * nw + k * nw + l;
#ifdef SWIFT_DEBUG_CHECKS
  //printf("Eagle cooling.h j, ny %d %d\n",j,ny);
  assert(i < nx);
  assert(j < ny);
  assert(k < nz);
  assert(l < nw);
#endif
  return index;
}


/*
 * ----------------------------------------------------------------------
 * This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation.  The table is assumed to
 * be evenly spaced.
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE void get_index_1d(float *table, int ntable, double x, int *i, float *dx) {
  float dxm1;
  const float EPS = 1.e-4;

  dxm1 = (float)(ntable - 1) / (table[ntable - 1] - table[0]);

  if ((float)x <= table[0] + EPS) {
    *i = 0;
    *dx = 0;
  } else if ((float)x >= table[ntable - 1] - EPS) {
    *i = ntable - 2;
    *dx = 1;
  } else {
    *i = (int)floor(((float)x - table[0]) * dxm1);
    if (*i >= ntable || *i < 0){
      printf("Eagle cooling.h i, ntable, x, table[0], dxm1 %d %d %.5e %.5e %.5e \n", *i, ntable, x, table[0], dxm1);
      fflush(stdout);
    }
    *dx = ((float)x - table[*i]) * dxm1;
  }
}

/*
 * ----------------------------------------------------------------------
 * Get cooling table redshift index
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE void get_redshift_index(float z, int *z_index, float *dz, const struct cooling_function_data* restrict cooling) {
  int i, iz;

  if (get_redshift_index_first_call == 0) {
    get_redshift_index_first_call = 1;
    get_redshift_index_previous = cooling->N_Redshifts - 2;

    /* this routine assumes cooling_redshifts table is in increasing order. Test
     * this. */
    for (i = 0; i < cooling->N_Redshifts - 2; i++)
      if (cooling->Redshifts[i + 1] < cooling->Redshifts[i]) {
        error("[get_redshift_index]: table should be in increasing order\n");
      }
  }

  /* before the earliest redshift or before hydrogen reionization, flag for
   * collisional cooling */
  if (z > cooling->reionisation_redshift) {
    *z_index = cooling->N_Redshifts;
    *dz = 0.0;
  }
  /* from reionization use the cooling tables */
  else if (z > cooling->Redshifts[cooling->N_Redshifts - 1] &&
           z <= cooling->reionisation_redshift) {
    *z_index = cooling->N_Redshifts + 1;
    *dz = 0.0;
  }
  /* at the end, just use the last value */
  else if (z <= cooling->Redshifts[0]) {
    *z_index = 0;
    *dz = 0.0;
  } else {
    /* start at the previous index and search */
    for (iz = get_redshift_index_previous; iz >= 0; iz--) {
      if (z > cooling->Redshifts[iz]) {
        *dz = (z - cooling->Redshifts[iz]) /
              (cooling->Redshifts[iz + 1] - cooling->Redshifts[iz]);

        get_redshift_index_previous = *z_index = iz;

        break;
      }
    }
  }
}


/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_1d(float *table, int i, float dx) {
  float result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE double interpol_1d_dbl(double *table, int i, float dx) {
  double result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_2d(float *table, int i, int j, float dx, float dy, int nx, int ny) {
  float result;
  int index[4];

  index[0] = row_major_index_2d(i,j,nx,ny);
  index[1] = row_major_index_2d(i,j+1,nx,ny);
  index[2] = row_major_index_2d(i+1,j,nx,ny);
  index[3] = row_major_index_2d(i+1,j+1,nx,ny);
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j %d, %d, table size %d\n", index[0],i,j,nx*ny);
  if(index[1] >= nx*ny || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j %d, %d, table size %d\n", index[1],i,j+1,nx*ny);
  if(index[2] >= nx*ny || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j %d, %d, table size %d\n", index[2],i+1,j,nx*ny);
  if(index[3] >= nx*ny || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j %d, %d, table size %d\n", index[3],i+1,j+1,nx*ny);
#endif

  result = (1 - dx) * (1 - dy) * table[index[0]] + (1 - dx) * dy * table[index[1]] +
           dx * (1 - dy) * table[index[2]] + dx * dy * table[index[3]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE double interpol_2d_dbl(double *table, int i, int j, double dx, double dy, int nx, int ny) {
  double result;
  int index[4];

  index[0] = row_major_index_2d(i,j,nx,ny);
  index[1] = row_major_index_2d(i,j+1,nx,ny);
  index[2] = row_major_index_2d(i+1,j,nx,ny);
  index[3] = row_major_index_2d(i+1,j+1,nx,ny);
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j %d, %d, table size %d\n", index[0],i,j,nx*ny);
  if(index[1] >= nx*ny || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j %d, %d, table size %d\n", index[1],i,j+1,nx*ny);
  if(index[2] >= nx*ny || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j %d, %d, table size %d\n", index[2],i+1,j,nx*ny);
  if(index[3] >= nx*ny || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j %d, %d, table size %d\n", index[3],i+1,j+1,nx*ny);
#endif

  result = (1 - dx) * (1 - dy) * table[index[0]] + (1 - dx) * dy * table[index[1]] +
           dx * (1 - dy) * table[index[2]] + dx * dy * table[index[3]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a tri-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_3d(float *table, int i, int j, int k, float dx, float dy,
                  float dz, int nx, int ny, int nz) {
  float result;
  int index[8];

  index[0] = row_major_index_3d(i,j,k,nx,ny,nz);
  index[1] = row_major_index_3d(i,j,k+1,nx,ny,nz);
  index[2] = row_major_index_3d(i,j+1,k,nx,ny,nz);
  index[3] = row_major_index_3d(i,j+1,k+1,nx,ny,nz);
  index[4] = row_major_index_3d(i+1,j,k,nx,ny,nz);
  index[5] = row_major_index_3d(i+1,j,k+1,nx,ny,nz);
  index[6] = row_major_index_3d(i+1,j+1,k,nx,ny,nz);
  index[7] = row_major_index_3d(i+1,j+1,k+1,nx,ny,nz);
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny*nz || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[0],i,j,k,nx*ny*nz);
  if(index[1] >= nx*ny*nz || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[1],i,j,k+1,nx*ny*nz);
  if(index[2] >= nx*ny*nz || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[2],i,j+1,k,nx*ny*nz);
  if(index[3] >= nx*ny*nz || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[3],i,j+1,k+1,nx*ny*nz);
  if(index[4] >= nx*ny*nz || index[4] < 0) fprintf(stderr,"index 4 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[4],i+1,j,k,nx*ny*nz);
  if(index[5] >= nx*ny*nz || index[5] < 0) fprintf(stderr,"index 5 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[5],i+1,j,k+1,nx*ny*nz);
  if(index[6] >= nx*ny*nz || index[6] < 0) fprintf(stderr,"index 6 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[6],i+1,j+1,k,nx*ny*nz);
  if(index[7] >= nx*ny*nz || index[7] < 0) fprintf(stderr,"index 7 out of bounds %d, i,j,k %d, %d, %d, table size %d\n", index[7],i+1,j+1,k+1,nx*ny*nz);
#endif

  result = (1 - dx) * (1 - dy) * (1 - dz) * table[index[0]] +
           (1 - dx) * (1 - dy) * dz * table[index[1]] +
           (1 - dx) * dy * (1 - dz) * table[index[2]] +
           (1 - dx) * dy * dz * table[index[3]] +
           dx * (1 - dy) * (1 - dz) * table[index[4]] +
           dx * (1 - dy) * dz * table[index[5]] +
           dx * dy * (1 - dz) * table[index[6]] +
           dx * dy * dz * table[index[7]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_4d(float *table, int i, int j, int k, int l, float dx,
                        float dy, float dz, float dw, int nx, int ny, int nz, int nw) {
  float result;
  int index[16];

  index[0]  = row_major_index_4d(i,j,k,l,nx,ny,nz,nw);
  index[1]  = row_major_index_4d(i,j,k,l+1,nx,ny,nz,nw);
  index[2]  = row_major_index_4d(i,j,k+1,l,nx,ny,nz,nw);
  index[3]  = row_major_index_4d(i,j,k+1,l+1,nx,ny,nz,nw);
  index[4]  = row_major_index_4d(i,j+1,k,l,nx,ny,nz,nw);
  index[5]  = row_major_index_4d(i,j+1,k,l+1,nx,ny,nz,nw);
  index[6]  = row_major_index_4d(i,j+1,k+1,l,nx,ny,nz,nw);
  index[7]  = row_major_index_4d(i,j+1,k+1,l+1,nx,ny,nz,nw);
  index[8]  = row_major_index_4d(i+1,j,k,l,nx,ny,nz,nw);
  index[9]  = row_major_index_4d(i+1,j,k,l+1,nx,ny,nz,nw);
  index[10] = row_major_index_4d(i+1,j,k+1,l,nx,ny,nz,nw);
  index[11] = row_major_index_4d(i+1,j,k+1,l+1,nx,ny,nz,nw);
  index[12] = row_major_index_4d(i+1,j+1,k,l,nx,ny,nz,nw);
  index[13] = row_major_index_4d(i+1,j+1,k,l+1,nx,ny,nz,nw);
  index[14] = row_major_index_4d(i+1,j+1,k+1,l,nx,ny,nz,nw);
  index[15] = row_major_index_4d(i+1,j+1,k+1,l+1,nx,ny,nz,nw);
#ifdef SWIFT_DEBUG_CHECKS
  if(index[0] >= nx*ny*nz*nw || index[0] < 0) fprintf(stderr,"index 0 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[0],i,j,k,l,nx*ny*nz*nw);
  if(index[1] >= nx*ny*nz*nw || index[1] < 0) fprintf(stderr,"index 1 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[1],i,j,k,l+1,nx*ny*nz*nw);
  if(index[2] >= nx*ny*nz*nw || index[2] < 0) fprintf(stderr,"index 2 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[2],i,j,k+1,l,nx*ny*nz*nw);
  if(index[3] >= nx*ny*nz*nw || index[3] < 0) fprintf(stderr,"index 3 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[3],i,j,k+1,l+1,nx*ny*nz*nw);
  if(index[4] >= nx*ny*nz*nw || index[4] < 0) fprintf(stderr,"index 4 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[4],i,j+1,k,l,nx*ny*nz*nw);
  if(index[5] >= nx*ny*nz*nw || index[5] < 0) fprintf(stderr,"index 5 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[5],i,j+1,k,l+1,nx*ny*nz*nw);
  if(index[6] >= nx*ny*nz*nw || index[6] < 0) fprintf(stderr,"index 6 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[6],i,j+1,k+1,l,nx*ny*nz*nw);
  if(index[7] >= nx*ny*nz*nw || index[7] < 0) fprintf(stderr,"index 7 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[7],i,j+1,k+1,l+1,nx*ny*nz*nw);
  if(index[8] >= nx*ny*nz*nw || index[8] < 0) fprintf(stderr,"index 8 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[8],i+1,j,k,l,nx*ny*nz*nw);
  if(index[9] >= nx*ny*nz*nw || index[9] < 0) fprintf(stderr,"index 9 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[9],i+1,j,k,l+1,nx*ny*nz*nw);
  if(index[10] >= nx*ny*nz*nw || index[10] < 0) fprintf(stderr,"index 10 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[10],i+1,j,k+1,l,nx*ny*nz*nw);
  if(index[11] >= nx*ny*nz*nw || index[11] < 0) fprintf(stderr,"index 11 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[11],i+1,j,k+1,l+1,nx*ny*nz*nw);
  if(index[12] >= nx*ny*nz*nw || index[12] < 0) fprintf(stderr,"index 12 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[12],i+1,j+1,k,l,nx*ny*nz*nw);
  if(index[13] >= nx*ny*nz*nw || index[13] < 0) fprintf(stderr,"index 13 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[13],i+1,j+1,k,l+1,nx*ny*nz*nw);
  if(index[14] >= nx*ny*nz*nw || index[14] < 0) fprintf(stderr,"index 14 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[14],i+1,j+1,k+1,l,nx*ny*nz*nw);
  if(index[15] >= nx*ny*nz*nw || index[15] < 0) fprintf(stderr,"index 15 out of bounds %d, i,j,k,l, %d, %d, %d, %d, table size %d\n", index[15],i+1,j+1,k+1,l+1,nx*ny*nz*nw);
#endif

  result = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table[index[0]] +
           (1 - dx) * (1 - dy) * (1 - dz) * dw * table[index[1]] +
           (1 - dx) * (1 - dy) * dz * (1 - dw) * table[index[2]] +
           (1 - dx) * (1 - dy) * dz * dw * table[index[3]] +
           (1 - dx) * dy * (1 - dz) * (1 - dw) * table[index[4]] +
           (1 - dx) * dy * (1 - dz) * dw * table[index[5]] +
           (1 - dx) * dy * dz * (1 - dw) * table[index[6]] +
           (1 - dx) * dy * dz * dw * table[index[7]] +
           dx * (1 - dy) * (1 - dz) * (1 - dw) * table[index[8]] +
           dx * (1 - dy) * (1 - dz) * dw * table[index[9]] +
           dx * (1 - dy) * dz * (1 - dw) * table[index[10]] +
           dx * (1 - dy) * dz * dw * table[index[11]] +
           dx * dy * (1 - dz) * (1 - dw) * table[index[12]] +
           dx * dy * (1 - dz) * dw * table[index[13]] +
           dx * dy * dz * (1 - dw) * table[index[14]] +
           dx * dy * dz * dw * table[index[15]];

  return result;
}

/*
 * Constructs 1d table dependent on temperature to interpolate from 3d table
 */
__attribute__((always_inline)) INLINE static void construct_1d_table_from_3d(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,int z_i, float d_z,int He_i,float d_He,int n_h_i,float d_n_h,double *result_table){
  int index[4];
  //int He_i, n_h_i, z_i;
  //float d_He, d_n_h, d_z;
  //float inHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  //float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;
  
  //get_redshift_index(cosmo->z,&z_i,&d_z,cooling);	
  //get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  //get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  for(int i = 0; i < cooling->N_Temp; i++){
    index[0] = row_major_index_3d(z_i,   n_h_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    index[1] = row_major_index_3d(z_i,   n_h_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    index[2] = row_major_index_3d(z_i+1, n_h_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    index[3] = row_major_index_3d(z_i+1, n_h_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_Temp);
    
    result_table[i] = (1 - d_z) * (1 - d_n_h) * table[index[0]] +
                      (1 - d_z) * d_n_h       * table[index[1]] +
                      d_z       * (1 - d_n_h) * table[index[2]] +
                      d_z       * d_n_h       * table[index[3]];
  }
}

/*
 * Constructs 1d table dependent on temperature to interpolate from 4d table
 * Used for getting the contribution from H and He to cooling and finding the
 * location of maximum gradient of the f(u) = u - u_0 - lambda*nh^2*dt
 */
__attribute__((always_inline)) INLINE static void construct_1d_table_from_4d_H_He(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,int z_i, float d_z,int He_i,float d_He,int n_h_i,float d_n_h,double *result_table, double *x_max_grad, double u_ini, float dt){
  int index[8];
  //int He_i, n_h_i, z_i;
  //float d_He, d_n_h, d_z;
  float x0, x1, y0, y1, old_grad, new_grad;
  //float inHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float ratefact = inn_h * (XH / eagle_proton_mass_cgs);

  //get_redshift_index(cosmo->z,&z_i,&d_z,cooling);
  //get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  //get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  old_grad = 0.0;

  for(int i = 0; i < cooling->N_Temp; i++){
    index[0] = row_major_index_4d(z_i,   n_h_i,   He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[1] = row_major_index_4d(z_i,   n_h_i,   He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[2] = row_major_index_4d(z_i,   n_h_i+1, He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[3] = row_major_index_4d(z_i,   n_h_i+1, He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[4] = row_major_index_4d(z_i+1, n_h_i,   He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[5] = row_major_index_4d(z_i+1, n_h_i,   He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[6] = row_major_index_4d(z_i+1, n_h_i+1, He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[7] = row_major_index_4d(z_i+1, n_h_i+1, He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);

    result_table[i] = (1 - d_z) * (1 - d_n_h) * (1 - d_He) * table[index[0]] +
                      (1 - d_z) * (1 - d_n_h) * d_He       * table[index[1]] +
                      (1 - d_z) * d_n_h       * (1 - d_He) * table[index[2]] +
                      (1 - d_z) * d_n_h       * d_He       * table[index[3]] +
                      d_z       * (1 - d_n_h) * (1 - d_He) * table[index[4]] +
                      d_z       * (1 - d_n_h) * d_He       * table[index[5]] +
                      d_z       * d_n_h       * (1 - d_He) * table[index[6]] +
                      d_z       * d_n_h       * d_He       * table[index[7]];
#if SWIFT_DEBUG_CHECKS
    if (isnan(result_table[i])) printf("Eagle cooling.h 1 i dz dnh dHe table values %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n",i,d_z,d_n_h,d_He,table[index[0]],table[index[1]],table[index[2]],table[index[3]],table[index[4]],table[index[5]],table[index[6]],table[index[7]]);
#endif
    if (i > 0 && dt > 0 && x_max_grad != NULL){
      x0 = pow(10.0,cooling->Therm[i-1]);
      x1 = pow(10.0,cooling->Therm[i]);
      y0 = x0 - u_ini - result_table[i-1]*ratefact*dt;
      y1 = x1 - u_ini - result_table[i]*ratefact*dt;
      //new_grad = (y1 - y0)/(x1 - x0);
      new_grad = log10(y1/y0)/log10(x1/x0);
      //printf("Eagle cooling.h gradient max gradient y1 y0 x1 x0 u_ini lambda %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n", new_grad, old_grad, y1, y0, x1, x0, u_ini, result_table[i-1]*ratefact*dt);
      if (new_grad > old_grad) {
        *x_max_grad = 0.5*(x0 + x1);
        old_grad = new_grad;
      }
    }
  }
    //printf("Eagle cooling.h max gradient, guess, u_ini %.5e %.5e %.5e\n", old_grad, *x_max_grad, u_ini);
}

/*
 * Constructs 1d table dependent on temperature to interpolate from 4d table
 */
__attribute__((always_inline)) INLINE static void construct_1d_table_from_4d(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,int z_i, float d_z,int He_i,float d_He,int n_h_i,float d_n_h,double *result_table){
  int index[8];
  //int He_i, n_h_i, z_i;
  //float d_He, d_n_h, d_z;
  //float inHe = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  //float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;
  
  //get_redshift_index(cosmo->z,&z_i,&d_z,cooling);	
  //get_index_1d(cooling->HeFrac, cooling->N_He, inHe, &He_i, &d_He);
  //get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  for(int i = 0; i < cooling->N_Temp; i++){
    index[0] = row_major_index_4d(z_i,   n_h_i,   He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[1] = row_major_index_4d(z_i,   n_h_i,   He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[2] = row_major_index_4d(z_i,   n_h_i+1, He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[3] = row_major_index_4d(z_i,   n_h_i+1, He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[4] = row_major_index_4d(z_i+1, n_h_i,   He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[5] = row_major_index_4d(z_i+1, n_h_i,   He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[6] = row_major_index_4d(z_i+1, n_h_i+1, He_i,   i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    index[7] = row_major_index_4d(z_i+1, n_h_i+1, He_i+1, i, cooling->N_Redshifts, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    
    result_table[i] = (1 - d_z) * (1 - d_n_h) * (1 - d_He) * table[index[0]] +
                      (1 - d_z) * (1 - d_n_h) * d_He       * table[index[1]] +
                      (1 - d_z) * d_n_h       * (1 - d_He) * table[index[2]] +
                      (1 - d_z) * d_n_h       * d_He       * table[index[3]] +
                      d_z       * (1 - d_n_h) * (1 - d_He) * table[index[4]] +
                      d_z       * (1 - d_n_h) * d_He       * table[index[5]] +
                      d_z       * d_n_h       * (1 - d_He) * table[index[6]] +
                      d_z       * d_n_h       * d_He       * table[index[7]];
#if SWIFT_DEBUG_CHECKS
    if (isnan(result_table[i])) printf("Eagle cooling.h 2 i dz dnh dHe table values %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n",i,d_z,d_n_h,d_He,table[index[0]],table[index[1]],table[index[2]],table[index[3]],table[index[4]],table[index[5]],table[index[6]],table[index[7]]);
#endif
  }
}

/*
 * Constructs 1d table dependent on temperature to interpolate from 4d table
 * Used for tabulating contribution to cooling from metals
 */
__attribute__((always_inline)) INLINE static void construct_1d_table_from_4d_elements(const struct part* restrict p,const struct cooling_function_data* restrict cooling,const struct cosmology* restrict cosmo, const struct phys_const *internal_const, float *table,int z_i, float d_z,int n_h_i,float d_n_h,double *result_table){
  int index[4];
  //int n_h_i, z_i;
  //float d_n_h, d_z;
  //float inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale;

  //get_redshift_index(cosmo->z,&z_i,&d_z,cooling);
  //get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

    for(int j = 0; j < cooling->N_Elements; j++){
  for(int i = 0; i < cooling->N_Temp; i++){
  if (i == 0) result_table[i] = 0.0;
      index[0] = row_major_index_4d(z_i,   j, n_h_i,   i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      index[1] = row_major_index_4d(z_i,   j, n_h_i+1, i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      index[2] = row_major_index_4d(z_i+1, j, n_h_i,   i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);
      index[3] = row_major_index_4d(z_i+1, j, n_h_i+1, i, cooling->N_Redshifts, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);

      result_table[i] += (1 - d_z) * (1 - d_n_h) * table[index[0]] +
                         (1 - d_z) * d_n_h       * table[index[1]] +
                         d_z       * (1 - d_n_h) * table[index[2]] +
                         d_z       * d_n_h       * table[index[3]];
// #if SWIFT_DEBUG_CHECKS
    if (isnan(result_table[i])) printf("Eagle cooling.h 3 i partial sums %d %.5e %.5e %.5e %.5e %.5e \n",i, result_table[i],(1 - d_z) * (1 - d_n_h) * table[index[0]],
                       (1 - d_z) * (1 - d_n_h) * table[index[0]] + (1 - d_z) * d_n_h * table[index[1]],
                       (1 - d_z) * (1 - d_n_h) * table[index[0]] + (1 - d_z) * d_n_h * table[index[1]] + d_z * (1 - d_n_h) * table[index[2]],
                       (1 - d_z) * (1 - d_n_h) * table[index[0]] + (1 - d_z) * d_n_h * table[index[1]] + d_z * (1 - d_n_h) * table[index[2]] + d_z * d_n_h * table[index[3]]);
// #endif
    }
  }
}

/*
 * Function which calculates ratio of particle abundances to solar
 * Imported from EAGLE, needs rewriting to run fast
 */
inline int set_cooling_SolarAbundances(const float *element_abundance,
                                double *cooling_element_abundance,
                                const struct cooling_function_data* restrict cooling,
				const struct part* restrict p) {
  int i, index;
  int Silicon_SPH_Index = -1;
  int Calcium_SPH_Index = -1;
  int Sulphur_SPH_Index = -1;
  
  int Silicon_CoolHeat_Index = -1;
  int Calcium_CoolHeat_Index = -1;
  int Sulphur_CoolHeat_Index = -1;

  //float *cooling->ElementAbundance_SOLARM1 = malloc(cooling->N_SolarAbundances*sizeof(float));

    /* determine (inverse of) solar abundance of these elements */
    for (i = 0; i < cooling->N_Elements; i++) {
      index =
          get_element_index(cooling->SolarAbundanceNames,
                            cooling->N_SolarAbundances, cooling->ElementNames[i]);

      if (index < 0) error("Eagle cooling.h index out of bounds");

      index = cooling->SolarAbundanceNamePointers[i];

      if(cooling->SolarAbundances[index] != 0) cooling->ElementAbundance_SOLARM1[i] = 1. / cooling->SolarAbundances[index];
      else cooling->ElementAbundance_SOLARM1[i] = 0.0;

    }

    /* Sulphur tracks Silicon: may choose not to follow Sulphur as SPH element
     */
    /* Same is true for Calcium */
    /* We will assume the code tracks Silicon, and may need to scale Calcium and
     * Sulphur accordingly */

    Silicon_SPH_Index = element_index("Silicon",cooling);
    Calcium_SPH_Index = element_index("Calcium",cooling);
    Sulphur_SPH_Index = element_index("Sulphur",cooling);

    Silicon_CoolHeat_Index =
        get_element_index(cooling->ElementNames, cooling->N_Elements, "Silicon");
    Calcium_CoolHeat_Index =
        get_element_index(cooling->ElementNames, cooling->N_Elements, "Calcium");
    Sulphur_CoolHeat_Index =
        get_element_index(cooling->ElementNames, cooling->N_Elements, "Sulphur");

    if (Silicon_CoolHeat_Index == -1 || Calcium_CoolHeat_Index == -1 ||
        Sulphur_CoolHeat_Index == -1) {
        error("Si, Ca, or S index out of bound\n");
    }

  int sili_index;
  for (i = 0; i < cooling->N_Elements; i++) {
    if (strcmp(chemistry_get_element_name((enum chemistry_element) i), "Silicon") == 0) sili_index = i;
  }

  // Eagle way of identifying and assigning element abundance with strange workaround for calcium and sulphur
  //for (i = 0; i < cooling->N_Elements; i++) {
  //  if (i == Calcium_CoolHeat_Index && Calcium_SPH_Index == -1)
  //    /* SPH does not track Calcium: use Si abundance */
  //    if (Silicon_SPH_Index == -1)
  //      cooling_element_abundance[i] = 0.0;
  //    else{
  //      cooling_element_abundance[i] =
  //          element_abundance[Silicon_SPH_Index] *
  //          cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
  //    }
  //  else if (i == Sulphur_CoolHeat_Index && Sulphur_SPH_Index == -1)
  //    /* SPH does not track Sulphur: use Si abundance */
  //    if (Silicon_SPH_Index == -1)
  //      cooling_element_abundance[i] = 0.0;
  //    else{
  //      cooling_element_abundance[i] =
  //          element_abundance[Silicon_SPH_Index] *
  //          cooling_ElementAbundance_SOLARM1[Silicon_CoolHeat_Index];
  //    }
  //  else{
  //    cooling_element_abundance[i] = element_abundance[cooling->ElementNamePointers[i]] *
  //                                   cooling_ElementAbundance_SOLARM1[i];
  //    //printf ("Eagle cooling.h element, name, abundance, solar abundance, solarm1, cooling abundance %d %s %.5e %.5e %.5e %.5e\n",cooling->ElementNamePointers[i],cooling->ElementNames[i], element_abundance[cooling->ElementNamePointers[i]],cooling->SolarAbundances[i], cooling_ElementAbundance_SOLARM1[i], cooling_element_abundance[i]);
  //  }
  //}
  
  for (i = 0; i < cooling->N_Elements; i++) {
    if (i == Calcium_CoolHeat_Index && Calcium_SPH_Index != -1)
      if (Silicon_SPH_Index == -1)
        cooling_element_abundance[i] = 0.0;
      else{
        cooling_element_abundance[i] =
            element_abundance[sili_index] * cooling->calcium_over_silicon_ratio *
            cooling->ElementAbundance_SOLARM1[Calcium_CoolHeat_Index];
      }
    else if (i == Sulphur_CoolHeat_Index && Sulphur_SPH_Index != -1)
      /* SPH does not track Sulphur: use Si abundance */
      if (Silicon_SPH_Index == -1)
        cooling_element_abundance[i] = 0.0;
      else{
        cooling_element_abundance[i] =
            element_abundance[sili_index] * cooling->sulphur_over_silicon_ratio *
            cooling->ElementAbundance_SOLARM1[Sulphur_CoolHeat_Index];
      }
    else{
      cooling_element_abundance[i] = element_abundance[cooling->ElementNamePointers[i]] *
                                     cooling->ElementAbundance_SOLARM1[i];
    }
    //printf ("Eagle cooling.h i, solar abundance name pointer, element, name, abundance, solarm1, cooling abundance %d %d %d %s %.5e %.5e %.5e\n", i, cooling->SolarAbundanceNamePointers[i],cooling->ElementNamePointers[i],cooling->ElementNames[i], element_abundance[cooling->ElementNamePointers[i]], cooling->ElementAbundance_SOLARM1[i], cooling_element_abundance[i]);
  }

  return 0;
}

/*
 * Replacement for set_cooling_SolarAbundances, work in progress
 */
__attribute__((always_inline)) INLINE float eagle_solar_abundance_factor(enum chemistry_element elem,
									 const struct part* restrict p,
									 const struct cooling_function_data* restrict cooling){
  float element_mass_fraction, solar_abundance;

  //if (elem == chemistry_element_S){
  //  element_mass_fraction = p->chemistry_data.metal_mass_fraction[chemistry_element_Si]*cooling->sulphur_over_silicon_ratio;
  //} else if (elem == chemistry_element_Ca){
  //  element_mass_fraction = p->chemistry_data.metal_mass_fraction[chemistry_element_Si]*cooling->calcium_over_silicon_ratio;
  //} else {
    element_mass_fraction = p->chemistry_data.metal_mass_fraction[elem];
    solar_abundance = cooling->SolarAbundances[elem];
  //}
  
  float element_abundance_factor = element_mass_fraction/solar_abundance;

  return element_abundance_factor;
}

/*
 * @brief interpolates temperature from internal energy based on table
 *
 */
__attribute__((always_inline)) INLINE static double eagle_convert_temp_to_u_1d_table(double temp,
										float *temperature_table,
										const struct part* restrict p,
										const struct cooling_function_data* restrict cooling,
										const struct cosmology* restrict cosmo,
										const struct phys_const *internal_const) {
  
  int temp_i;
  float d_temp, u;

  get_index_1d(temperature_table, cooling->N_Temp, log10(temp), &temp_i, &d_temp);

  u = pow(10.0,interpol_1d(cooling->Therm,temp_i,d_temp));

  return u;
}

/*
 * @brief interpolates temperature from internal energy based on table
 *
 */
__attribute__((always_inline)) INLINE static double eagle_convert_u_to_temp_1d_table(double u,
										float *delta_u,
										double *temperature_table,
										const struct part* restrict p,
										const struct cooling_function_data* restrict cooling,
										const struct cosmology* restrict cosmo,
										const struct phys_const *internal_const) {
  
  int u_i;
  float d_u, logT, T;

  get_index_1d(cooling->Therm, cooling->N_Temp, log10(u), &u_i, &d_u);

  logT = interpol_1d_dbl(temperature_table,u_i,d_u);
  T = pow(10.0, logT);

  if (u_i == 0 && d_u == 0) T *= u / pow(10.0, cooling->Therm[0]);

  *delta_u = pow(10.0,cooling->Therm[u_i + 1]) - pow(10.0,cooling->Therm[u_i]);

  return T;
}

/*
 * @brief calculates cooling rate from multi-d tables
 *
 */
__attribute__((always_inline)) INLINE static double eagle_metal_cooling_rate(double u,
									     double *temp_table,
									     int z_index,
									     float dz,
									     int n_h_i,
									     float d_n_h,
									     int He_i,
									     float d_He,
									     const struct part* restrict p, 
									     const struct cooling_function_data* restrict cooling, 
									     const struct cosmology* restrict cosmo, 
									     const struct phys_const *internal_const, 
									     double* element_lambda) {
  double T_gam, solar_electron_abundance;
  double n_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale; // chemistry_data
  double z = cosmo->z;
  double cooling_rate = 0.0, temp_lambda;
  float du;
  float h_plus_he_electron_abundance;

  int i;
  double temp;
  int temp_i;
  float d_temp;
  //float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  //get_redshift_index(z,&z_index,&dz,cooling);	
  
#if SWIFT_DEBUG_CHECKS
  if (isnan(u)) printf("u is nan id %llu\n", p->id);
#endif
  temp = eagle_convert_u_to_temp_1d_table(u,&du,temp_table,p,cooling,cosmo,internal_const);

  get_index_1d(cooling->Temp, cooling->N_Temp, log10(temp), &temp_i, &d_temp);
  //get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  //get_index_1d(cooling->nH, cooling->N_nH, log10(n_h), &n_h_i, &d_n_h);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

    /* Collisional cooling */
    //element_lambda[0] =
    //    interpol_2d(cooling->table.collisional_cooling.H_plus_He_heating, He_i,
    //                 temp_i, d_He, d_temp,cooling->N_He,cooling->N_Temp);
    //h_plus_he_electron_abundance =
    //    interpol_2d(cooling->table.collisional_cooling.H_plus_He_electron_abundance, He_i,
    //                 temp_i, d_He, d_temp,cooling->N_He,cooling->N_Temp);
    
    /* Photodissociation */
    //element_lambda[0] =
    //    interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_heating, He_i,
    //                 temp_i, n_h_i, d_He, d_temp, d_n_h,cooling->N_He,cooling->N_Temp,cooling->N_nH);
    //h_plus_he_electron_abundance =
    //    interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_electron_abundance, He_i,
    //                 temp_i, n_h_i, d_He, d_temp, d_n_h,cooling->N_He,cooling->N_Temp,cooling->N_nH);

    /* redshift tables */
    temp_lambda = interpol_4d(cooling->table.element_cooling.H_plus_He_heating, z_index, n_h_i, He_i, temp_i, dz, d_n_h, d_He, d_temp,cooling->N_Redshifts,cooling->N_nH,cooling->N_He,cooling->N_Temp);
    h_plus_he_electron_abundance = interpol_4d(cooling->table.element_cooling.H_plus_He_electron_abundance, z_index, n_h_i, He_i, temp_i, dz, d_n_h, d_He, d_temp,cooling->N_Redshifts,cooling->N_nH,cooling->N_He,cooling->N_Temp);
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[0] = temp_lambda;

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {
    /* inverse Compton cooling is not in collisional table
       before reionisation so add now */

    T_gam = eagle_cmb_temperature * (1 + z);

    temp_lambda = -eagle_compton_rate * (temp - eagle_cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                 h_plus_he_electron_abundance / n_h;
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

    /* for each element, find the abundance and multiply it
       by the interpolated heating-cooling */

    //set_cooling_SolarAbundances(p->chemistry_data.metal_mass_fraction, cooling->solar_abundances, cooling, p);

    /* Collisional cooling */
    //solar_electron_abundance =
    //    interpol_1d(cooling->table.collisional_cooling.electron_abundance, temp_i, d_temp); /* ne/n_h */

    //for (i = 0; i < cooling->N_Elements; i++){
    //    element_lambda[i+2] = interpol_2d(cooling->table.collisional_cooling.metal_heating, i,
    //                    temp_i, 0.0, d_temp,cooling->N_Elements,cooling->N_Temp) *
    //        (h_plus_he_electron_abundance / solar_electron_abundance) *
    //        cooling->solar_abundances[i];
    //}
    
    /* Photodissociation */
    //solar_electron_abundance =
    //    interpol_2d(cooling->table.photodissociation_cooling.electron_abundance, temp_i, n_h_i, d_temp, d_n_h, cooling->N_Temp, cooling->N_nH); /* ne/n_h */
      
    //for (i = 0; i < cooling->N_Elements; i++){
    //    element_lambda[i+2] = interpol_3d(cooling->table.photodissociation_cooling.metal_heating, i,
    //                    temp_i, n_h_i, 0.0, d_temp, d_n_h,cooling->N_Elements,cooling->N_Temp,cooling->N_nH) *
    //        (h_plus_he_electron_abundance / solar_electron_abundance) *
    //        cooling->solar_abundances[i];
    //}
    
    /* redshift tables */
    solar_electron_abundance = interpol_3d(cooling->table.element_cooling.electron_abundance, z_index,n_h_i,dz,temp_i,d_n_h,d_temp,cooling->N_Redshifts,cooling->N_nH,cooling->N_Temp);
    
    for (i = 0; i < cooling->N_Elements; i++){
    	// WARNING: before doing proper runs need to 
	// multiply by ratio of particle abundances to solar
	temp_lambda = interpol_4d(cooling->table.element_cooling.metal_heating,z_index,i,n_h_i,temp_i,dz,0.0,d_n_h,d_temp,cooling->N_Redshifts,cooling->N_Elements,cooling->N_nH,cooling->N_Temp) *
            (h_plus_he_electron_abundance / solar_electron_abundance);
        cooling_rate += temp_lambda;
        if (element_lambda != NULL) element_lambda[i+2] = temp_lambda;
    }

    return cooling_rate;
}


/*
 * @brief calculates cooling rate from 1d tables
 *
 */
__attribute__((always_inline)) INLINE static double eagle_metal_cooling_rate_1d_table(double u,
										      double *dlambda_du,
									              double *H_plus_He_heat_table,
										      double *H_plus_He_electron_abundance_table,
										      double *element_cooling_table,
										      double *element_electron_abundance_table,
										      double *temp_table,
										      int z_index,
										      float dz,
										      int n_h_i,
										      float d_n_h,
										      int He_i,
										      float d_He,
										      const struct part* restrict p, 
										      const struct cooling_function_data* restrict cooling, 
										      const struct cosmology* restrict cosmo, 
										      const struct phys_const *internal_const, 
										      double* element_lambda) {
  double T_gam, solar_electron_abundance;
  double n_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,internal_const)*cooling->number_density_scale; // chemistry_data
  double z = cosmo->z;
  double cooling_rate = 0.0, temp_lambda;
  float du;
  float h_plus_he_electron_abundance;

  int i;
  double temp;
  int temp_i;
  float d_temp;
  //float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He]/(p->chemistry_data.metal_mass_fraction[chemistry_element_H]+p->chemistry_data.metal_mass_fraction[chemistry_element_He]);

  *dlambda_du = 0.0;
  
  //get_redshift_index(z,&z_index,&dz,cooling);	
  
#if SWIFT_DEBUG_CHECKS
  if (isnan(u)) printf("u is nan id %llu\n", p->id);
#endif
  temp = eagle_convert_u_to_temp_1d_table(u,&du,temp_table,p,cooling,cosmo,internal_const);

  get_index_1d(cooling->Temp, cooling->N_Temp, log10(temp), &temp_i, &d_temp);
  //get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  //get_index_1d(cooling->nH, cooling->N_nH, log10(n_h), &n_h_i, &d_n_h);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

    /* Collisional cooling */
    //element_lambda[0] =
    //    interpol_2d(cooling->table.collisional_cooling.H_plus_He_heating, He_i,
    //                 temp_i, d_He, d_temp,cooling->N_He,cooling->N_Temp);
    //h_plus_he_electron_abundance =
    //    interpol_2d(cooling->table.collisional_cooling.H_plus_He_electron_abundance, He_i,
    //                 temp_i, d_He, d_temp,cooling->N_He,cooling->N_Temp);
    
    /* Photodissociation */
    //element_lambda[0] =
    //    interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_heating, He_i,
    //                 temp_i, n_h_i, d_He, d_temp, d_n_h,cooling->N_He,cooling->N_Temp,cooling->N_nH);
    //h_plus_he_electron_abundance =
    //    interpol_3d(cooling->table.photodissociation_cooling.H_plus_He_electron_abundance, He_i,
    //                 temp_i, n_h_i, d_He, d_temp, d_n_h,cooling->N_He,cooling->N_Temp,cooling->N_nH);

    /* redshift tables */
    temp_lambda = interpol_1d_dbl(H_plus_He_heat_table, temp_i, d_temp);
    h_plus_he_electron_abundance = interpol_1d_dbl(H_plus_He_electron_abundance_table, temp_i, d_temp);
    cooling_rate += temp_lambda;
    *dlambda_du += (((double) H_plus_He_heat_table[temp_i+1]) - ((double) H_plus_He_heat_table[temp_i]))/((double) du);
    if (element_lambda != NULL) element_lambda[0] = temp_lambda;

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {
    /* inverse Compton cooling is not in collisional table
       before reionisation so add now */

    T_gam = eagle_cmb_temperature * (1 + z);

    temp_lambda = -eagle_compton_rate * (temp - eagle_cmb_temperature * (1 + z)) * pow((1 + z), 4) *
                 h_plus_he_electron_abundance / n_h;
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

    /* for each element, find the abundance and multiply it
       by the interpolated heating-cooling */

    //set_cooling_SolarAbundances(p->chemistry_data.metal_mass_fraction, cooling->solar_abundances, cooling, p);

    /* Collisional cooling */
    //solar_electron_abundance =
    //    interpol_1d(cooling->table.collisional_cooling.electron_abundance, temp_i, d_temp); /* ne/n_h */

    //for (i = 0; i < cooling->N_Elements; i++){
    //    element_lambda[i+2] = interpol_2d(cooling->table.collisional_cooling.metal_heating, i,
    //                    temp_i, 0.0, d_temp,cooling->N_Elements,cooling->N_Temp) *
    //        (h_plus_he_electron_abundance / solar_electron_abundance) *
    //        cooling->solar_abundances[i];
    //}
    
    /* Photodissociation */
    //solar_electron_abundance =
    //    interpol_2d(cooling->table.photodissociation_cooling.electron_abundance, temp_i, n_h_i, d_temp, d_n_h, cooling->N_Temp, cooling->N_nH); /* ne/n_h */
      
    //for (i = 0; i < cooling->N_Elements; i++){
    //    element_lambda[i+2] = interpol_3d(cooling->table.photodissociation_cooling.metal_heating, i,
    //                    temp_i, n_h_i, 0.0, d_temp, d_n_h,cooling->N_Elements,cooling->N_Temp,cooling->N_nH) *
    //        (h_plus_he_electron_abundance / solar_electron_abundance) *
    //        cooling->solar_abundances[i];
    //}
    
    /* redshift tables */
    solar_electron_abundance = interpol_1d_dbl(element_electron_abundance_table, temp_i, d_temp);
    
    for (i = 0; i < cooling->N_Elements; i++){
    	// WARNING: before doing proper runs need to 
	// multiply by ratio of particle abundances to solar
	temp_lambda = interpol_1d_dbl(element_cooling_table, temp_i, d_temp) *
            (h_plus_he_electron_abundance / solar_electron_abundance);
        cooling_rate += temp_lambda;
	*dlambda_du += (((double) element_cooling_table[temp_i+1])*((double) H_plus_He_electron_abundance_table[temp_i+1])/((double) element_electron_abundance_table[temp_i+1]) - 
			((double) element_cooling_table[temp_i])*((double) H_plus_He_electron_abundance_table[temp_i])/((double) element_electron_abundance_table[temp_i]))/((double) du);
        if (element_lambda != NULL) element_lambda[i+2] = temp_lambda;
    }

    return cooling_rate;
}

/*
 * Wrapper function previously used to calculate cooling rate and dLambda_du
 * Can be used to check whether calculation of dLambdaa_du is done correctly.
 * Should be removed when finished
 */
__attribute__((always_inline)) INLINE static double eagle_cooling_rate_1d_table(double u,
										double *dLambdaNet_du,
									        double *H_plus_He_heat_table,
										double *H_plus_He_electron_abundance_table,
										double *element_cooling_table,
										double *element_electron_abundance_table,
										double *temp_table,
										int z_index,
										float dz,
										int n_h_i,
										float d_n_h,
										int He_i,
										float d_He,
										const struct part* restrict p, 
										const struct cooling_function_data* restrict cooling, 
										const struct cosmology* restrict cosmo, 
										const struct phys_const *internal_const) {
  double *element_lambda = NULL, lambda_net1 = 0.0;//, lambda_net2 = 0.0, delta, dLambdaNet_du_calc;
  //float d_u;
  //int u_i;
  //get_index_1d(cooling->Therm, cooling->N_Temp, log10(u), &u_i, &d_u);
  //if (u_i >= cooling->N_Temp-2){
  //  delta = pow(10.0,cooling->Therm[cooling->N_Temp-2]) - u;
  //} else {
  //  delta = pow(10.0,cooling->Therm[u_i+1]) - pow(10.0,cooling->Therm[u_i]);
  //}

  //lambda_net1 = eagle_metal_cooling_rate_1d_table(pow(10.0,cooling->Therm[u_i]), dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);
  //lambda_net2 = eagle_metal_cooling_rate_1d_table(pow(10.0,cooling->Therm[u_i+1]), dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, internal_const, element_lambda);
  //
  //dLambdaNet_du_calc = (lambda_net2 - lambda_net1)/delta;

  //printf( "Eagle cooling.h 2 dlambda_du, dlambda, du, index %0.5e %0.5e %.5e %d\n", dLambdaNet_du_calc, lambda_net2 - lambda_net1, delta, u_i);
  
  lambda_net1 = eagle_metal_cooling_rate_1d_table(u, dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, internal_const, element_lambda);

  //printf("Eagle cooling.h double check dlambda_du endpoints, lambda_nets %.5e %.5e\n", *dLambdaNet_du, dLambdaNet_du_calc);


  return lambda_net1;
}

/* 
 * Same as eagle_metal_cooling_rate_1d_table, but prints contribution to
 * cooling rate from each element for confirmation against Wiersma et al 2009
 * idl routines
 */
__attribute__((always_inline)) INLINE static double eagle_print_metal_cooling_rate_1d_table(double *H_plus_He_heat_table,
											    double *H_plus_He_electron_abundance_table,
											    double *element_cooling_table,
											    double *element_electron_abundance_table,
											    double *temp_table,
							       				    int z_index,
							       				    float dz,
							       				    int n_h_i,
							       				    float d_n_h,
							       				    int He_i,
							       				    float d_He,
											    const struct part* restrict p, 
											    const struct cooling_function_data* restrict cooling, 
											    const struct cosmology* restrict cosmo, 
											    const struct phys_const *internal_const) {
  double *element_lambda, lambda_net = 0.0;
  element_lambda = malloc((cooling->N_Elements+2)*sizeof(double));
  double u = hydro_get_physical_internal_energy(p,cosmo)*cooling->internal_energy_scale;
  double dLambdaNet_du;
  
  char output_filename[21];
  FILE** output_file = malloc((cooling->N_Elements+2)*sizeof(FILE*));
  for (int element = 0; element < cooling->N_Elements+2; element++){
    sprintf(output_filename, "%s%d%s", "cooling_output_",element,".dat");
    output_file[element] = fopen(output_filename, "a");
    if (output_file == NULL)
    {   
        printf("Error opening file!\n");
        exit(1);
    }
  }

  for (int j = 0; j < cooling->N_Elements+2; j++) element_lambda[j] = 0.0;
  lambda_net = eagle_metal_cooling_rate_1d_table(u,&dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, internal_const, element_lambda);
  for (int j = 0; j < cooling->N_Elements+2; j++) {
    fprintf(output_file[j],"%.5e\n",element_lambda[j]);
  }
  
  for (int i = 0; i < cooling->N_Elements+2; i++) fclose(output_file[i]);

  return lambda_net;
}

/*
 * Work in progress more bisection scheme
 */
__attribute__((always_inline)) INLINE static float bisection_iter(float x_init,
                                                               double u_ini,
                                                               double *H_plus_He_heat_table,
                                                               double *H_plus_He_electron_abundance_table,
                                                               double *element_cooling_table,
                                                               double *element_electron_abundance_table,
                                                               double *temp_table,
							       int z_index,
							       float dz,
							       int n_h_i,
							       float d_n_h,
							       int He_i,
							       float d_He,
                                                               struct part* restrict p,
                                                               const struct cosmology* restrict cosmo,
                                                               const struct cooling_function_data* restrict cooling,
                                                               const struct phys_const* restrict phys_const,
                                                               float dt){
  double x_a, x_b, x_next, f_a, f_b, f_next, dLambdaNet_du;
  int i = 0;
  float shift = 0.3;

  x_a = x_init;
  x_b = cooling->Therm[0]/log10(exp(1.0));
  f_a = eagle_cooling_rate_1d_table(exp(x_a), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
  f_b = eagle_cooling_rate_1d_table(exp(x_b), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);

  while (f_a*f_b >= 0 & i < 10*eagle_max_iterations) {
    if ((x_a+shift)/log(10) < cooling->Therm[cooling->N_Temp-1]) {
      x_a += shift;
      f_a = eagle_cooling_rate_1d_table(exp(x_a), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
    }
    //if ((x_b-shift)/log(10) > cooling->Therm[0]) {
    //  x_b -= shift;
    //  if (isnan(exp(x_b))) printf("4 u is nan id %llu\n",p->id);
    //  f_b = eagle_cooling_rate_1d_table(exp(x_b), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, p, cooling, cosmo, phys_const);
    //}
    i++;
  }
#if SWIFT_DEBUG_CHECKS
  assert(f_a*f_b < 0);
#endif

  i = 0;
  do{
    x_next = 0.5*(x_a + x_b);
    f_next = eagle_cooling_rate_1d_table(exp(x_next), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
    if (f_a*f_next < 0) {
      x_b = x_next;
      f_b = f_next;
    } else {
      x_a = x_next;
      f_a = f_next;
    }
    i++;
    //printf("Eagle cooling.h id x_a x_b f_a f_b i %llu %.5e %.5e %.5e %.5e %d\n", p->id, x_a, x_b, f_a, f_b, i );
  } while (fabs(x_a - x_b) > 1.0e-2 && i < eagle_max_iterations);
  
  return x_b;
}

/*
 * Work in progress more stable integration scheme
 * that may be used for particles which do not
 * converge using newton_iter
 */
__attribute__((always_inline)) INLINE static float brent_iter(float x_init,
                                                               double u_ini,
                                                               double *H_plus_He_heat_table,
                                                               double *H_plus_He_electron_abundance_table,
                                                               double *element_cooling_table,
                                                               double *element_electron_abundance_table,
                                                               double *temp_table,
							       int z_index,
							       float dz,
							       int n_h_i,
							       float d_n_h,
							       int He_i,
							       float d_He,
                                                               struct part* restrict p,
                                                               const struct cosmology* restrict cosmo,
                                                               const struct cooling_function_data* restrict cooling,
                                                               const struct phys_const* restrict phys_const,
                                                               float dt){
  /* this routine does the iteration scheme, call one and it iterates to convergence */

  double x_a,x_b,x_c,x_d,x_s, x_temp, delta = 1.0e-2;
  double dLambdaNet_du, Lambda_a, Lambda_b, Lambda_c, Lambda_s;
  double f_a,f_b,f_c,f_s, f_temp;
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction in Hydrogen number density */
  double inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,phys_const)*cooling->number_density_scale;

  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  double ratefact = inn_h * (XH / eagle_proton_mass_cgs);

  /* set helium and hydrogen reheating term */
  //LambdaTune = eagle_helium_reionization_extraheat(z_index, dz); // INCLUDE HELIUM REIONIZATION????
  //if (zold > z_index) {
  //  LambdaCumul += LambdaTune;
  //  printf(" EagleDoCooling %d %g %g %g\n", z_index, dz, LambdaTune, LambdaCumul);
  //  zold = z_index;
  //}

  int i = 0;
  x_a = x_init;
  //x_b = cooling->Therm[0]*log(10);
  x_b = 0.5*x_init;
  x_c = x_a;
  x_d = 0.0;
  Lambda_a = eagle_cooling_rate_1d_table(exp(x_a), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
  Lambda_b = eagle_cooling_rate_1d_table(exp(x_b), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
  f_a = exp(x_a) - exp(x_init) - Lambda_a*ratefact*dt;
  f_b = exp(x_b) - exp(x_init) - Lambda_b*ratefact*dt;
  assert(f_a*f_b < 0);
  if (f_a < f_b) {
    x_temp = x_a; x_a = x_b; x_b = x_temp;
    f_temp = f_a; f_a = f_b; f_b = f_temp;
  }
  int mflag = 1;

  do /* iterate to convergence */
    {
      Lambda_c = eagle_cooling_rate_1d_table(exp(x_c), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
      f_c = exp(x_c) - exp(x_init) - Lambda_c*ratefact*dt;

      if ((f_a != f_c) && (f_b != f_c)) {
        x_s = x_a*f_b*f_c/((f_a - f_b)*(f_a - f_c)) +
            x_b*f_a*f_c/((f_b - f_a)*(f_b - f_c)) +
            x_c*f_b*f_a/((f_c - f_b)*(f_c - f_a));
        printf("Eagle cooling.h 1 id x_s, u_s, error, %llu %.5e %.5e %.5e %d\n", p->id, x_s, exp(x_s), x_a - x_b, i);
        printf("Eagle cooling.h 1 id u_a, u_b, u_c, f_a, f_b, f_c, %llu %.5e %.5e %.5e %.5e %.5e %.5e %.5e %d\n", p->id, exp(x_a), exp(x_b), exp(x_c), f_a, f_b, f_c, x_a - x_b, i);
      } else {
        x_s = x_b - f_b*(x_b - x_a)/(f_b - f_a);
        printf("Eagle cooling.h 2 id x_s, u_s, error, %llu %.5e %.5e %.5e %d\n", p->id, x_s, exp(x_s), x_a - x_b, i);
      }

      if ( (x_s < (3.0*x_a + x_b)/4.0 || x_s > x_b) ||
           (mflag == 1 && fabs(x_s - x_b) >= 0.5*fabs(x_b - x_c)) ||
           (mflag != 1 && fabs(x_s - x_b) >= 0.5*fabs(x_c - x_d)) ||
           (mflag == 1 && fabs(x_b - x_c) < delta) ||
           (mflag != 1 && fabs(x_c - x_d) < delta) ) {
        x_s = 0.5*(x_a + x_b);
        printf("Eagle cooling.h 3 id x_s, u_s, error, %llu %.5e %.5e %.5e %d\n", p->id, x_s, exp(x_s), x_a - x_b, i);
        mflag = 1;
      } else {
        mflag = 0;
      }
      x_d = x_c;
      x_c = x_b;
      printf("Eagle cooling.h 4 id x_s, u_s, error, %llu %.5e %.5e %.5e %d\n", p->id, x_s, exp(x_s), x_a - x_b, i);
      Lambda_s = eagle_cooling_rate_1d_table(exp(x_s), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
      f_s = exp(x_s) - exp(x_init) - Lambda_s*ratefact*dt;
      if (f_s*f_a < 0) {
        x_b = x_s;
      } else {
        x_a = x_s;
      }
      Lambda_a = eagle_cooling_rate_1d_table(exp(x_a), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
      Lambda_b = eagle_cooling_rate_1d_table(exp(x_b), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
      f_a = exp(x_a) - exp(x_init) - Lambda_a*ratefact*dt;
      f_b = exp(x_b) - exp(x_init) - Lambda_b*ratefact*dt;
      if (f_a < f_b) {
        x_temp = x_a; x_a = x_b; x_b = x_temp;
        f_temp = f_a; f_a = f_b; f_b = f_temp;
      }
      i++;

    } while (fabs(x_b - x_a) > delta && i < eagle_max_iterations);

    if (i >= eagle_max_iterations) {
      n_eagle_cooling_rate_calls_3++;
      printf("Eagle cooling.h particle %llu with density %.5e, initial energy, %.5e and redshift %.5e not converged \n", p->id, inn_h, exp(x_init), cosmo->z);
    }

  return x_s;

}


/*
 * Performs newton iteration until convergence or exceeding
 * allowed number of iterations.
 */
__attribute__((always_inline)) INLINE static float newton_iter(float x_init,
							       double u_ini,
							       double *H_plus_He_heat_table,
							       double *H_plus_He_electron_abundance_table,
							       double *element_cooling_table,
							       double *element_electron_abundance_table,
							       double *temp_table,
							       int z_index,
							       float dz,
							       int n_h_i,
							       float d_n_h,
							       int He_i,
							       float d_He,
    							       struct part* restrict p,
    							       const struct cosmology* restrict cosmo,
    							       const struct cooling_function_data* restrict cooling,
    							       const struct phys_const* restrict phys_const,
							       float dt,
							       int *printflag,
							       float relax_factor){
  /* this routine does the iteration scheme, call one and it iterates to convergence */

  double x, x_old;
  double dLambdaNet_du, LambdaNet;
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction in Hydrogen number density */
  double inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,phys_const)*cooling->number_density_scale;
  
  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  double ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  float residual ;  // added by RGB
  
  /* set helium and hydrogen reheating term */
  //LambdaTune = eagle_helium_reionization_extraheat(z_index, dz); // INCLUDE HELIUM REIONIZATION????
  //if (zold > z_index) {
  //  LambdaCumul += LambdaTune;
  //  printf(" EagleDoCooling %d %g %g %g\n", z_index, dz, LambdaTune, LambdaCumul);
  //  zold = z_index;
  //}
  
  x_old = x_init ;
  x = x_old;
  int i = 0;

  do /* iterate to convergence */
    {
      x_old = x ;
      // this version is needed when reionization is included...
      // LambdaNet = (LambdaTune / (dt * ratefact)) + eagle_cooling_rate_1d_table(exp(x_old), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
      LambdaNet = eagle_cooling_rate_1d_table(exp(x_old), &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);
      n_eagle_cooling_rate_calls_1++;
      
      // Newton iterate the variable.
      x = x_old - relax_factor*(1.0 - u_ini*exp(-x_old) - LambdaNet*ratefact*dt*exp(-x_old))/(1.0 - dLambdaNet_du*ratefact*dt);

      int table_out_of_bounds = 0;
      //if (x > cooling->Therm[cooling->N_Temp - 1]*log(10) + 1) {
      if (exp(x) > 1.0e19) {
        //x = log(u_ini);
        table_out_of_bounds = 1;
      }
      //if (x < cooling->Therm[0]*log(10) - 1) {
      if (exp(x) < 1.0e9) {
        //x = log(u_ini);
        table_out_of_bounds = 1;
      }
      
      residual = exp(x_old) - u_ini - LambdaNet*ratefact*dt ;
      if(dt > 0 && *printflag == 1) {
        printf("z, n_h, x, x_old, u_ini, num, denom, LNet, dL_du  %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d %d\n", 
          	  cosmo->z, inn_h, exp(x), exp(x_old), u_ini,(1.0 - u_ini*exp(-x_old) - LambdaNet*ratefact*dt*exp(-x_old)),(1.0 - dLambdaNet_du*ratefact*dt), LambdaNet, dLambdaNet_du, x - x_old, i, table_out_of_bounds);
        table_out_of_bounds = 0;
      }

      i++;
      if (table_out_of_bounds == 1) i = eagle_max_iterations;
    } while (fabs(x - x_old) > 1.0e-2 && i < eagle_max_iterations);
  if (i >= eagle_max_iterations) {
    n_eagle_cooling_rate_calls_3++;
    if (*printflag == 0) *printflag = 1;
  }

  return x;

}

/**
 * @brief Apply the cooling function to a particle.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo cosmology struct
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void cooling_cool_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    struct part* restrict p, struct xpart* restrict xp, float dt) {
  
  double u_old = hydro_get_physical_internal_energy(p,cosmo)*cooling->internal_energy_scale;
  float dz; 
  int z_index;
  get_redshift_index(cosmo->z,&z_index,&dz,cooling);

  float XH, HeFrac;
  double inn_h;

  double ratefact, u, LambdaNet, LambdaTune = 0, dLambdaNet_du;

  static double zold = 100, LambdaCumul = 0;
  dt *= units_cgs_conversion_factor(us,UNIT_CONV_TIME);

  u = u_old;
  double u_ini = u_old, u_temp;

  XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] / (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  //printf("Eagle cooling.h density %.5e\n", hydro_get_physical_density(p,cosmo)*units_cgs_conversion_factor(us,UNIT_CONV_DENSITY));

  /* convert Hydrogen mass fraction in Hydrogen number density */
  inn_h = chemistry_get_number_density(p,cosmo,chemistry_element_H,phys_const)*cooling->number_density_scale;
  //float inn_h_eagle = hydro_get_physical_density(p,cosmo)*units_cgs_conversion_factor(us,UNIT_CONV_DENSITY) * XH / eagle_proton_mass_cgs;
  /* ratefact = inn_h * inn_h / rho; Might lead to round-off error: replaced by
   * equivalent expression  below */
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  /* set helium and hydrogen reheating term */
  //LambdaTune = eagle_helium_reionization_extraheat(z_index, dz); // INCLUDE HELIUM REIONIZATION????
  if (zold > z_index) {
    LambdaCumul += LambdaTune;
    printf(" EagleDoCooling %d %g %g %g\n", z_index, dz, LambdaTune, LambdaCumul);
    zold = z_index;
  }
  
  int He_i, n_h_i;
  float d_He, d_n_h;
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  double temp_table[176]; 				// WARNING sort out how it is declared/allocated
  construct_1d_table_from_4d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.temperature,z_index,dz,He_i,d_He,n_h_i,d_n_h,temp_table);

  // To be used when amount of cooling is small, put back in when ready
  LambdaNet = eagle_metal_cooling_rate(u_ini, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const,NULL);
  if (fabs(ratefact * LambdaNet * dt) < 0.05 * u_old) {
    /* cooling rate is small, take the explicit solution */
    u = u_old + ratefact * LambdaNet * dt;
  } else {
    // construct 1d table of cooling rates wrt temperature
    double H_plus_He_heat_table[176]; 			// WARNING sort out how it is declared/allocated
    double H_plus_He_electron_abundance_table[176]; 	// WARNING sort out how it is declared/allocated
    double element_cooling_table[176]; 			// WARNING sort out how it is declared/allocated
    double element_electron_abundance_table[176]; 		// WARNING sort out how it is declared/allocated
    construct_1d_table_from_4d_H_He(p,cooling,cosmo,phys_const,cooling->table.element_cooling.H_plus_He_heating,z_index,dz,He_i,d_He,n_h_i,d_n_h,H_plus_He_heat_table, &u_temp, u_ini, dt);
    construct_1d_table_from_4d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.H_plus_He_electron_abundance,z_index,dz,He_i,d_He,n_h_i,d_n_h,H_plus_He_electron_abundance_table);
    construct_1d_table_from_4d_elements(p,cooling,cosmo,phys_const,cooling->table.element_cooling.metal_heating,z_index,dz,n_h_i,d_n_h,element_cooling_table);
    construct_1d_table_from_3d(p,cooling,cosmo,phys_const,cooling->table.element_cooling.electron_abundance,z_index,dz,He_i,d_He,n_h_i,d_n_h,element_electron_abundance_table);

    n_eagle_cooling_rate_calls_2++;

    LambdaNet = eagle_cooling_rate_1d_table(u_ini, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, z_index, dz, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo, phys_const);

    double u_eq = 5.0e12;
    u_temp = u_ini;


    if (dt > 0) {
      float u_temp_guess = u_ini + LambdaNet*ratefact*dt;
      //printf("Eagle cooling.h id u_temp_guess, u_ini, LambdaNet, ratefact, dt %llu %.5e %.5e %.5e %.5e %.5e \n", p->id,u_temp_guess, u_ini, LambdaNet, ratefact, dt);
      //if ((LambdaNet < 0 && u_temp < u_temp_guess) ||
      //    (LambdaNet >= 0 && u_temp > u_temp_guess))
      //  u_temp = u_temp_guess;
      u_temp = u_temp_guess;
      if ((LambdaNet < 0 && u_temp < u_eq) ||
          (LambdaNet >= 0 && u_temp > u_eq))
        u_temp = u_eq;
      if (isnan(u_temp)) u_temp = u_eq;
    
      //printf("Eagle cooling.h id, u_temp, n_h, z %llu %.5e %.5e %.5e\n", p->id, u_temp, inn_h, cosmo->z);

      int printflag = 0;
      float relax_factor = 1.0;
      float x = newton_iter(log(u_temp),u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,z_index,dz,n_h_i,d_n_h,He_i,d_He,p,cosmo,cooling,phys_const,dt,&printflag,relax_factor);
      if (printflag == 1) {
        //relax_factor = 0.75;
        //x = newton_iter(log(u_temp),u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,z_index,dz,n_h_i,d_n_h,He_i,d_He,p,cosmo,cooling,phys_const,dt,&printflag,relax_factor);
        x = bisection_iter(log(u_temp),u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,z_index,dz,n_h_i,d_n_h,He_i,d_He,p,cosmo,cooling,phys_const,dt);
        n_eagle_cooling_rate_calls_4++;
        printflag = 2;
      }
      u = exp(x);
    }
  }

  float cooling_du_dt = 0.0;
  if (dt > 0){ 
    cooling_du_dt = (u - u_ini)/dt/cooling->power_scale;
  }

  const float hydro_du_dt = hydro_get_internal_energy_dt(p);

  /* Update the internal energy time derivative */
  hydro_set_internal_energy_dt(p, hydro_du_dt + cooling_du_dt);

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy += -hydro_get_mass(p) * cooling_du_dt * (dt/units_cgs_conversion_factor(us,UNIT_CONV_TIME));

}


/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Cooling Model", "EAGLE");
}

/**
 * @brief Computes the cooling time-step.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float cooling_timestep(
    const struct cooling_function_data* restrict cooling,
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us, const struct part* restrict p) {

  /* Remember to update when using an implicit integrator!!!*/
  //const float cooling_rate = cooling->cooling_rate;
  //const float internal_energy = hydro_get_comoving_internal_energy(p);
  //return cooling->cooling_tstep_mult * internal_energy / fabsf(cooling_rate);

  return FLT_MAX;
}

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param cooling The properties of the cooling function.
 */
__attribute__((always_inline)) INLINE static void cooling_first_init_part(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct cooling_function_data* cooling) {

  xp->cooling_data.radiated_energy = 0.f; 
}

/**
 * @brief Returns the total radiated energy by this particle.
 *
 * @param xp The extended particle data
 */
__attribute__((always_inline)) INLINE static float cooling_get_radiated_energy(
    const struct xpart* restrict xp) {

  return xp->cooling_data.radiated_energy;
}

/**
 * @brief Initialises the cooling properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cooling The cooling properties to initialize
 */
static INLINE void cooling_init_backend(
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const,
    struct cooling_function_data* cooling) {
  
  char fname[200];

  parser_get_param_string(parameter_file, "EagleCooling:filename",cooling->cooling_table_path);
  cooling->reionisation_redshift = parser_get_param_float(parameter_file, "EagleCooling:reionisation_redshift");
  cooling->calcium_over_silicon_ratio = parser_get_param_float(parameter_file, "EAGLEChemistry:CalciumOverSilicon");
  cooling->sulphur_over_silicon_ratio = parser_get_param_float(parameter_file, "EAGLEChemistry:SulphurOverSilicon");
  GetCoolingRedshifts(cooling);
  sprintf(fname, "%sz_0.000.hdf5", cooling->cooling_table_path);
  ReadCoolingHeader(fname,cooling);
  MakeNamePointers(cooling);
  cooling->table = eagle_readtable(cooling->cooling_table_path,cooling);
  printf("Eagle cooling.h read table \n");

  cooling->ElementAbundance_SOLARM1 = malloc(cooling->N_SolarAbundances*sizeof(float));
  cooling->solar_abundances = malloc(cooling->N_Elements*sizeof(double));

  cooling->delta_u = cooling->Therm[1] - cooling->Therm[0];

  cooling->internal_energy_scale = units_cgs_conversion_factor(us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
  cooling->number_density_scale = units_cgs_conversion_factor(us,UNIT_CONV_DENSITY)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
  cooling->power_scale = units_cgs_conversion_factor(us,UNIT_CONV_POWER)/units_cgs_conversion_factor(us,UNIT_CONV_MASS);
  cooling->temperature_scale = units_cgs_conversion_factor(us,UNIT_CONV_TEMPERATURE);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * @param cooling The properties of the cooling function.
 */
static INLINE void cooling_print_backend(
    const struct cooling_function_data* cooling) {

  message("Cooling function is 'EAGLE'.");
}

#endif /* SWIFT_COOLING_EAGLE_H */
