/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_ROW_MAJOR_ID_H
#define SWIFT_ROW_MAJOR_ID_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "inline.h"

/**
 * @brief Returns 1D index of a 3D NxNxN array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices is >= N
 * or < 0.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline, const)) INLINE static int row_major_id_periodic(
    const int i, const int j, const int k, const int N) {

  return (((i + N) % N) * N * N + ((j + N) % N) * N + ((k + N) % N));
}


/**
 * @brief Returns 1D index of a FFTW-padded 3D NxNxN array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices is >= N
 * or < 0. Note that indexes are of type size_t in this version and the array
 * is assumed to be padded to size 2*(N/2+1) in the last dimension as required
 * by FFTW MPI routines.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline, const)) INLINE static size_t
row_major_id_periodic_size_t_padded(const int i, const int j, const int k,
                                    const int N) {
  /* Find last two dimensions of the padded array */
  const size_t Nj = N;
  const size_t Nk = 2*(N/2+1);
  
  /* Get box-wrapped coordinates (note that we don't use the padded size here) */
  const size_t i_wrap = (size_t) ((i + N) % N);
  const size_t j_wrap = (size_t) ((j + N) % N);
  const size_t k_wrap = (size_t) ((k + N) % N);

  /* Compute index in the padded array */
  return (i_wrap*Nj*Nk) + (j_wrap*Nk) + k_wrap;
}


/**
 * @brief Return i coordinate from an id returned by row_major_id_periodic_size_t_padded
 *
 * This extracts the index in the first dimension from a row major id
 * returned by row_major_id_periodic_size_t_padded. I.e. it finds the
 * 'i' input parameter that was used to generate the id.
 *
 * @param id The padded row major ID.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline, const)) INLINE static int
get_xcoord_from_padded_row_major_id(const size_t id, const int N) {
  const size_t Nj = N;
  const size_t Nk = 2*(N/2+1);
  return (int) (id / (Nj*Nk));
}


/**
 * @brief Convert a global mesh array index to local slice index
 * 
 * Given an index into the padded N*N*2*(N/2+1) array, compute
 * the corresponding index in the slice of the array stored
 * on the local MPI rank.
 *
 * @param id The padded row major ID.
 * @param N Size of the array along one axis.
 * @param slice_offset Index of the first slice on this rank
 */
__attribute__((always_inline, const)) INLINE static size_t
get_index_in_local_slice(const size_t id, const int N, const int slice_offset) {
  const size_t Nj = N;
  const size_t Nk = 2*(N/2+1);
  return id - ((size_t) slice_offset)*Nj*Nk;
}

#endif /* SWIFT_ROW_MAJOR_ID_H */
