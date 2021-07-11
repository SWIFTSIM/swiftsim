/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#ifndef SWIFT_LIGHTCONE_MAP_H
#define SWIFT_LIGHTCONE_MAP_H

/* Standard headers */
#include <math.h>
#include <limits.h>

/* Config parameters. */
#include "../config.h"

/* HDF5 */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/* Define this to compute expected sum of pixels for consistency check:
   this lets us check that the sum of the updates buffered when particles
   cross the lightcone equals the sum over pixels in the final map.
   This is a test of the communication and smoothing routines. */
#define LIGHTCONE_MAP_CHECK_TOTAL

/* Local headers */
#include "error.h"
#include "healpix_smoothing.h"
#include "particle_buffer.h"
#include "parser.h"
#include "threadpool.h"
#include "units.h"

/**
 * @brief Struct to store a contribution to a healpix map pixel
 */
struct lightcone_map_contribution {

  /*! Amount to contribute to the pixel */
  double value;

  /*! Position on the sphere: spherical coordinates encoded as ints */
  int itheta, iphi;
  
  /*! Smoothing radius */
  float radius;

  /*! MPI ranks to send this contribution to */
  unsigned short first_dest, last_dest;

};


/**
 * @brief Struct to store a single lightcone healpix map
 */
struct lightcone_map {

  /*! Healpix nside parameter */
  int nside;

  /*! Pointer to struct with C++ Healpix_base class instance  */
  struct healpix_smoothing_info *smoothing_info;

  /*! Buffer to store contributions to the healpix map */
  struct particle_buffer buffer;

  /*! Block size in the particle buffer */
  size_t elements_per_block;

  /*! Total pixels in the map */
  size_t total_nr_pix;

  /*! Number of pixels stored on this node */
  size_t local_nr_pix;
  
  /*! Offset of the firts pixel stored on this rank */
  size_t local_pix_offset;

  /*! Number of pixels per rank (last node has any extra) */
  size_t pix_per_rank;

  /*! Local healpix map data */
  double *data;

  /*! Inner radius */
  double r_min;

  /*! Outer radius */
  double r_max;

  /*! Units of this map */
  enum unit_conversion_factor units;

  /*! Size of each pixel in steradians */
  double pixel_area_steradians;

  /*! Whether to smooth this map */
  int smooth;

#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  /*! Sum of the quantity accumulated to this map, used for consistency check */
  double sum;
#endif

  /*! MPI communicator info */
  int comm_rank, comm_size;

};


/**
 * @brief Store an angle in the range 0 to 2pi in an int
 *
 * @param theta The angle to store
 *
 */
__attribute__((always_inline)) INLINE static int angle_to_int(double theta) {
  
  if((theta < 0.0) || (theta > 2*M_PI))error("angle must be in range 0 to 2pi");

  const int nmax = (INT_MAX-1); /* -1 so we don't overflow for theta=2pi */
  const double dtheta = (2.0*M_PI)/(((double) nmax)+1.0);
  const double inv_dtheta = 1.0 / dtheta;
  return (int) (theta*inv_dtheta);
}


/**
 * @brief Convert an int back to an angle
 *
 * @param i The int to be interpreted as an angle
 *
 */
__attribute__((always_inline)) INLINE static double int_to_angle(int i) {
  
  const int nmax = (INT_MAX-1);
  const double dtheta = (2.0*M_PI)/(((double) nmax)+1.0);
  return dtheta*i+0.5*dtheta;
}


/**
 * @brief Add a value to the buffer for a healpix map
 *
 * @param map the #lightcone_map to update
 * @param pixel the pixel index to update
 * @param value the value to add
 *
 */
__attribute__((always_inline)) INLINE static void lightcone_map_buffer_update(struct lightcone_map *map,
                                                                              const double *pos,
                                                                              const double radius,
                                                                              const double value) {

  double theta, phi;
  healpix_smoothing_vec2ang(map->smoothing_info, pos, &theta, &phi);

  struct lightcone_map_contribution contr;
  contr.itheta = angle_to_int(theta);
  contr.iphi   = angle_to_int(phi);
  contr.radius = radius;
  contr.value = value;
  particle_buffer_append(&map->buffer, &contr);

  /* Sum values added to map for consistency check */
#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  atomic_add_d(&map->sum, value);
#endif

}


void lightcone_map_init(struct lightcone_map *map, const int nside,
                        const double r_min, const double r_max,
                        const size_t elements_per_block,
                        enum unit_conversion_factor units,
                        const int smooth);

void lightcone_map_clean(struct lightcone_map *map);

void lightcone_map_struct_dump(const struct lightcone_map *map, FILE *stream);

void lightcone_map_struct_restore(struct lightcone_map *map, FILE *stream);

void lightcone_map_allocate_pixels(struct lightcone_map *map, const int zero_pixels);

void lightcone_map_free_pixels(struct lightcone_map *map);

void lightcone_map_update_from_buffer(struct lightcone_map *map, struct threadpool *tp,
                                      const int verbose);

#ifdef HAVE_HDF5
void lightcone_map_write(struct lightcone_map *map, const hid_t loc_id, const char *name,
                         const struct unit_system *internal_units,
                         const struct unit_system *snapshot_units);
#endif

#endif /* #ifndef SWIFT_LIGHTCONE_MAP_H */
