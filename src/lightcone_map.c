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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* This object's header. */
#include "lightcone_map.h"

/* Local headers */
#include "align.h"
#include "common_io.h"
#include "error.h"
#include "exchange_structs.h"
#include "hydro.h"
#include "memuse.h"
#include "particle_buffer.h"
#include "restart.h"
#include "healpix_smoothing.h"

/* HDF5 */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif


#ifdef WITH_MPI
static int pixel_to_rank(struct lightcone_map *map, size_t pixel) {
  int rank = pixel / map->pix_per_rank;
  if(rank >= map->comm_size)rank = map->comm_size-1;
  return rank;
}
#endif

void lightcone_map_init(struct lightcone_map *map, const int nside,
                        const double r_min, const double r_max,
                        const size_t elements_per_block,
                        enum unit_conversion_factor units,
                        const int smooth) {
  
  int comm_rank = 0, comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif
  map->comm_size = comm_size;
  map->comm_rank = comm_rank;

  /* Initialise C++ code for smoothing on the sphere */
  map->smoothing_info = healpix_smoothing_init(nside, kernel_gamma);

  /* Initialise the data buffer for this map */
  particle_buffer_init(&map->buffer, sizeof(struct lightcone_map_contribution),
                       elements_per_block, "lightcone_map");

  /* Determine number of pixels in the map */  
  map->total_nr_pix = healpix_smoothing_get_npix(map->smoothing_info);

  /* Determine the area of a pixel */
  map->pixel_area_steradians = 4*M_PI/(map->total_nr_pix);

  /* Determine which pixels are stored on which rank:
     put pix_per_rank on each node with any extra assigned to
     the last node. This makes it easy to convert a pixel index
     to a node index. */
  map->pix_per_rank = map->total_nr_pix / comm_size;
  if(map->pix_per_rank==0)error("Must have healpix npix > number of MPI ranks!");
  if(comm_rank < comm_size-1)
    map->local_nr_pix = map->pix_per_rank;
  else
    map->local_nr_pix = map->total_nr_pix - (comm_size-1)*map->pix_per_rank;

  /* Store offset from local array index to global healpix pixel index  */
  map->local_pix_offset = map->pix_per_rank * comm_rank;
  
  /* Pixel data is initially not allocated */
  map->data = NULL;

  /* Record block size so we can re-initialise particle_buffer on restarting */
  map->elements_per_block = elements_per_block;
  
  /* Store resolution parameter, shell size, units */
  map->nside = nside;
  map->r_min = r_min;
  map->r_max = r_max;
  map->units = units;
  map->smooth = smooth;
#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  map->sum = 0.0;
#endif
}


/**
 * @brief Deallocate the lightcone_map pixel data
 *
 * @param map the #lightcone_map structure
 */
void lightcone_map_clean(struct lightcone_map *map) {
  
  healpix_smoothing_clean(map->smoothing_info);
  particle_buffer_free(&map->buffer);
  if(map->data)lightcone_map_free_pixels(map);
}


/**
 * @brief Allocate (and maybe initialize) the lightcone_map pixel data
 *
 * @param map the #lightcone_map structure
 * @param zero_pixels if true, set allocated pixels to zero
 */
void lightcone_map_allocate_pixels(struct lightcone_map *map, const int zero_pixels) {
  
  if(swift_memalign("lightcone_map_pixels", (void **) &map->data,
                    SWIFT_STRUCT_ALIGNMENT, sizeof(double)*map->local_nr_pix) != 0)
    error("Failed to allocate lightcone map pixel data");

  if(zero_pixels) {
    for(size_t i=0; i<map->local_nr_pix; i+=1)
      map->data[i] = 0.0;
  }

}


void lightcone_map_free_pixels(struct lightcone_map *map) {
  
  swift_free("lightcone_map_pixels", (void *) map->data);
  map->data = NULL;

}


/**
 * @brief Dump lightcone_map struct to the output stream.
 *
 * @param map the #lightcone_map structure
 * @param stream The stream to write to.
 */
void lightcone_map_struct_dump(const struct lightcone_map *map, FILE *stream) {

  /* Don't write the particle_buffer (must flush before dumping) */
  struct lightcone_map tmp = *map;
  memset(&tmp.buffer, 0, sizeof(struct particle_buffer));

  /* Don't write pointer to C++ smoothing info */
  tmp.smoothing_info = NULL;

  /* Write the struct */
  restart_write_blocks((void *) &tmp, sizeof(struct lightcone_map), 1, stream,
                       "lightcone_map", "lightcone_map");

  /* Write the pixel data if it is allocated */
  if(tmp.data)
    restart_write_blocks((void *) tmp.data, sizeof(double), tmp.local_nr_pix, 
                         stream, "lightcone_map_data", "lightcone_map_data");
}


/**
 * @brief Restore lightcone_map struct from the input stream.
 *
 * @param map the #lightcone_map structure
 * @param stream The stream to read from.
 */
void lightcone_map_struct_restore(struct lightcone_map *map, FILE *stream) {

  /* Read the struct */
  restart_read_blocks((void *)map, sizeof(struct lightcone_map), 1, stream,
                      NULL, "lightcone_map");

  /* Initialise the buffer for this map */
  particle_buffer_init(&map->buffer, sizeof(struct lightcone_map_contribution),
                       map->elements_per_block, "lightcone_map");
  
  /* Read the pixel data if it was allocated.
     map->data from the restart file is not a valid pointer now but we can
     check if it is not null to see if the pixel data block was written out. */
  if(map->data) {
    lightcone_map_allocate_pixels(map, /* zero_pixels = */ 0);
    restart_read_blocks((void *)map->data, sizeof(double), map->local_nr_pix,
                        stream, NULL, "lightcone_map");
  }

  /* Re-initialise C++ smoothing info */
  map->smoothing_info = healpix_smoothing_init(map->nside, kernel_gamma);

}

/**
 * @brief Mapper function for updating the healpix map
 *
 * @param map_data Pointer to an array of #lightcone_map_contribution
 * @param num_elements Number of elements in map_data
 * @param extra_data Pointer to the #lightcone_map struct
 *
 */
void healpix_smoothing_mapper(void *map_data, int num_elements,
                              void *extra_data) {
  
  /* Find the array of updates to apply */
  const struct lightcone_map_contribution *contr = 
    (struct lightcone_map_contribution *) map_data;

  /* Get a pointer to the lightcone_map struct to update */
  struct lightcone_map *map = (struct lightcone_map *) extra_data;

  /* Loop over updates to apply */
  for(size_t i=0; i<num_elements; i+=1) {

    /* Unpack angular coordinates */
    const double theta = int_to_angle(contr[i].itheta);
    const double phi   = int_to_angle(contr[i].iphi);

    /* Determine smoothing radius */
    double radius;
    if(map->smooth)
      radius = contr[i].radius;
    else
      radius = 0.0;

    /* Add this contribution to the map */
    healpix_smoothing_add_to_map(map->smoothing_info, theta, phi,
                                 radius, contr[i].value,
                                 map->local_pix_offset, map->local_nr_pix,
                                 map->data);
  }
}

#ifdef WITH_MPI

struct buffer_block_info {

  /*! Pointer to the buffer block */
  struct particle_buffer_block *block;

  /*! Number of elements from this block to go to each MPI rank */
  size_t *count;

  /*! Offsets at which to write elements in the send buffer */
  size_t *offset;
};




static void count_elements_to_send_mapper(void *map_data, int num_elements,
                                          void *extra_data) {
  
  struct buffer_block_info *block_info = (struct buffer_block_info *) map_data;
  struct lightcone_map *map = (struct lightcone_map *) extra_data;
  const int comm_size = map->comm_size;
  
  for(int block_nr=0; block_nr<num_elements; block_nr+=1) {
    
    /* Find the count and offset for this block */
    size_t *count = block_info[block_nr].count;
    size_t *offset = block_info[block_nr].offset;

    /* Get a pointer to the block itself */
    struct particle_buffer_block *block = block_info[block_nr].block;

    /* Initialise count and offset into the send buffer for this block */
    for(int i=0; i<comm_size; i+=1) {
      count[i] = 0;
      offset[i] = 0;
    }
    
    /* Loop over lightcone map contributions in this block */
    struct lightcone_map_contribution *contr =
      (struct lightcone_map_contribution *) block->data;
    for(size_t i=0; i<block->num_elements; i+=1) {

      /* Unpack angular coordinates */
      double theta = int_to_angle(contr[i].itheta);
      double phi   = int_to_angle(contr[i].iphi);

      /* Determine smoothing radius */
      double radius;
      if(map->smooth)
        radius = contr[i].radius;
      else
        radius = 0.0;

      /* Determine which ranks this contribution goes to */
      size_t first_pixel, last_pixel;
      healpix_smoothing_get_pixel_range(map->smoothing_info, theta, phi,
                                        radius, &first_pixel, &last_pixel);
      contr[i].first_dest = pixel_to_rank(map, first_pixel);
      contr[i].last_dest = pixel_to_rank(map, last_pixel);

      /* Update the counts for this block */
      for(int dest=contr[i].first_dest; dest<=contr[i].last_dest; dest+=1)
        count[dest] += 1;
    }
    
    /* Next block */
  }
}


static void store_elements_to_send_mapper(void *map_data, int num_elements,
                                          void *extra_data) {
  
  /* Unpack data we need */
  struct lightcone_map_contribution *sendbuf = (struct lightcone_map_contribution *) extra_data;
  struct buffer_block_info *block_info = (struct buffer_block_info *) map_data;
  
  /* Loop over blocks to process on this call */
  for(int block_nr=0; block_nr<num_elements; block_nr+=1) {
    
    /* Find the count and offset for this block */
    size_t *offset = block_info[block_nr].offset;

    /* Get a pointer to the block itself */
    struct particle_buffer_block *block = block_info[block_nr].block;
    
    /* Loop over lightcone map contributions in this block */
    struct lightcone_map_contribution *contr =
      (struct lightcone_map_contribution *) block->data;
    for(size_t i=0; i<block->num_elements; i+=1) {

      /* Determine which ranks this contribution goes to */
      const int first_dest = contr[i].first_dest;
      const int last_dest = contr[i].last_dest;

      /* Store this contribution to the send buffer (possibly multiple times) */
      for(int dest=first_dest; dest<=last_dest; dest+=1) {
        memcpy(sendbuf+offset[dest], contr+i, sizeof(struct lightcone_map_contribution));
        offset[dest] += 1;
      }

      /* Next element in this block */
    }
    /* Next block */
  }
}

#endif

/**
 * @brief Apply buffered updates to the healpix map
 *
 * @param map the #lightcone_map structure
 */
void lightcone_map_update_from_buffer(struct lightcone_map *map,
                                      struct threadpool *tp,
                                      const int verbose) {
  
#ifdef WITH_MPI

  /* Get MPI rank, number of ranks */
  const int comm_size = map->comm_size;

  /* Count data blocks and ensure number of elements is in range */
  size_t nr_blocks = 0;
  struct particle_buffer_block *block = map->buffer.first_block;
  while(block) {
    if(block->num_elements > map->buffer.elements_per_block)
      block->num_elements = map->buffer.elements_per_block;
    nr_blocks += 1;
    block = block->next;
  }

  /* Allocate array with counts and offsets for each block */
  struct buffer_block_info *block_info = malloc(sizeof(struct buffer_block_info)*nr_blocks);

  /* Initialize array of blocks */
  nr_blocks = 0;
  block = map->buffer.first_block;
  while(block) {
    block_info[nr_blocks].block = block;
    block_info[nr_blocks].count = malloc(sizeof(size_t)*comm_size);
    block_info[nr_blocks].offset = malloc(sizeof(size_t)*comm_size);
    nr_blocks += 1;
    block = block->next;
  }

  /* For each block, count how many elements are to be sent to each MPI rank */
  threadpool_map(tp, count_elements_to_send_mapper, block_info, nr_blocks,
                 sizeof(struct buffer_block_info), 1, map);
  
  /* Find total number of elements to go to each rank */
  size_t *send_count = malloc(sizeof(size_t)*comm_size);
  for(int i=0; i<comm_size; i+=1)
    send_count[i] = 0;
  for(size_t block_nr=0; block_nr<nr_blocks; block_nr+=1) {
    for(int i=0; i<comm_size; i+=1)
      send_count[i] += block_info[block_nr].count[i];
  }

  /* Find offset to the first element to go to each rank if we sort them by destination */
  size_t *send_offset = malloc(sizeof(size_t)*comm_size);
  send_offset[0] = 0;
  for(int i=1; i<comm_size; i+=1) {
    send_offset[i] = send_offset[i-1] + send_count[i-1];
  }

  /* For each block, find the location in the send buffer where we need to
     place the first element to go to each MPI rank */
  for(size_t block_nr=0; block_nr<nr_blocks; block_nr+=1) {
    for(int i=0; i<map->comm_size; i+=1) {
      if(block_nr==0) {
        /* This is the first block */
        block_info[block_nr].offset[i] = send_offset[i];
      } else {
        /* Not first, so elements are written after those of the previous block */
        block_info[block_nr].offset[i] = block_info[block_nr-1].offset[i] +
          block_info[block_nr-1].count[i];
      }
    }
  }
  
  /* Find the total number of elements to be sent */
  size_t total_nr_send = 0;
  for(int i=0; i<comm_size; i+=1)
    total_nr_send += send_count[i];

  /* Allocate the send buffer */
  struct lightcone_map_contribution *sendbuf = 
    malloc(sizeof(struct lightcone_map_contribution)*total_nr_send);
  
  /* Populate the send buffer */
  threadpool_map(tp, store_elements_to_send_mapper, block_info, nr_blocks,
                 sizeof(struct buffer_block_info), 1, sendbuf);

  /* We no longer need the array of blocks */
  for(size_t block_nr=0; block_nr<nr_blocks; block_nr+=1) {
    free(block_info[block_nr].count);
    free(block_info[block_nr].offset);
  }
  free(block_info);

  /* Empty the particle buffer now that we copied the data from it */
  particle_buffer_empty(&map->buffer);

  /* Determine number of elements to receive */
  size_t *recv_count = malloc(comm_size*sizeof(size_t));
  MPI_Alltoall(send_count, sizeof(size_t), MPI_BYTE, recv_count, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t total_nr_recv = 0;
  for(int i=0; i<comm_size; i+=1)
    total_nr_recv += recv_count[i];
  
  /* Allocate receive buffer */
  struct lightcone_map_contribution *recvbuf = 
    malloc(sizeof(struct lightcone_map_contribution)*total_nr_recv);
  
  /* Exchange data */
  exchange_structs(send_count, sendbuf, recv_count, recvbuf,
                   sizeof(struct lightcone_map_contribution));

  /* Apply received updates to the healpix map */
  threadpool_map(tp, healpix_smoothing_mapper, recvbuf, total_nr_recv,
                 sizeof(struct lightcone_map_contribution), 
                 threadpool_auto_chunk_size, map);

  /* Tidy up */
  free(send_count);
  free(send_offset);
  free(sendbuf);
  free(recv_count);
  free(recvbuf);

#else
  
  /* If not using MPI, we can update the map directly from the buffer */
  struct particle_buffer_block *block = NULL;
  struct lightcone_map_contribution *contr;
  size_t num_elements;
  do {
    particle_buffer_iterate(&map->buffer, &block, &num_elements, (void **) &contr);
    threadpool_map(tp, healpix_smoothing_mapper, contr, num_elements,
                   sizeof(struct lightcone_map_contribution), 
                   threadpool_auto_chunk_size, map);
  } while(block);
  particle_buffer_empty(&map->buffer);
    
#endif

}

#ifdef HAVE_HDF5
/**
 * @brief Write a lightcone map to a HDF5 file
 *
 * @param map the #lightcone_map structure
 * @param loc a HDF5 file or group identifier to write to
 * @param name the name of the dataset to create
 */
void lightcone_map_write(struct lightcone_map *map, const hid_t loc_id, const char *name,
                         const struct unit_system *internal_units,
                         const struct unit_system *snapshot_units) {

  /* Find unit conversion factor for this quantity */
  const double conversion_factor =
    units_conversion_factor(internal_units, snapshot_units, map->units);
  
  /* Convert units if necessary */
  if(conversion_factor != 1.0) {
    for(size_t i=0; i<map->local_nr_pix; i+=1)
      map->data[i] *= conversion_factor;
  }

  /* Create dataspace in memory corresponding to local pixels */
  const hsize_t mem_dims[1] = {(hsize_t) map->local_nr_pix};
  hid_t mem_space_id = H5Screate_simple(1, mem_dims, NULL);
  if(mem_space_id < 0)error("Unable to create memory dataspace");
  
  /* Create dataspace in the file corresponding to the full map */
  const hsize_t file_dims[1] = {(hsize_t) map->total_nr_pix};
  hid_t file_space_id = H5Screate_simple(1, file_dims, NULL);
  if(file_space_id < 0)error("Unable to create file dataspace");

  /* Select the part of the dataset in the file to write to */
#ifdef WITH_MPI
#ifdef HAVE_PARALLEL_HDF5
  const int comm_rank = map->comm_rank;
  const size_t pixel_offset = map->pix_per_rank * comm_rank;
  const hsize_t start[1] = {(hsize_t) pixel_offset};
  const hsize_t count[1] = {(hsize_t) map->local_nr_pix};
  if(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, NULL, count, NULL) < 0)
    error("Unable to select part of file dataspace to write to");
#else
  error("Writing lightcone maps with MPI requires parallel HDF5");
#endif
#endif

  /* Create the dataset */
  hid_t dset_id = H5Dcreate(loc_id, name, H5T_NATIVE_DOUBLE, file_space_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(dset_id < 0)error("Unable to create dataset %s", name);
    
  /* Write attributes */
  io_write_attribute_i(dset_id, "nside", map->nside);
  io_write_attribute_l(dset_id, "number_of_pixels", map->total_nr_pix);
  io_write_attribute_s(dset_id, "pixel_ordering_scheme", "ring");
  io_write_attribute_d(dset_id, "comoving_inner_radius", map->r_min);
  io_write_attribute_d(dset_id, "comoving_outer_radius", map->r_max);

#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  /* For consistency check: write sum of quantities accumulated to this map */
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &map->sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  io_write_attribute_d(dset_id, "expected_sum", map->sum);
#endif

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, map->units, 0.f);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, map->units);
  io_write_attribute_f(dset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(dset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(dset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(dset_id, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(dset_id, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(dset_id, "h-scale exponent", 0.f);
  io_write_attribute_f(dset_id, "a-scale exponent", 0.f);
  io_write_attribute_s(dset_id, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double cgs_factor = units_cgs_conversion_factor(snapshot_units, map->units);
  io_write_attribute_d(dset_id,
                       "Conversion factor to CGS (not including cosmological corrections)",
                       cgs_factor);

  /* Set up property list for the write */
  hid_t h_plist_id = H5Pcreate(H5P_DATASET_XFER);
#if defined(WITH_MPI)
  if(H5Pset_dxpl_mpio(h_plist_id, H5FD_MPIO_COLLECTIVE) < 0)
    error("Unable to set collective transfer mode");
#endif

  /* Write the data */
  if(H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id,
              h_plist_id, map->data) < 0)
    error("Unable to write dataset %s", name);

  /* Tidy up */
  H5Dclose(dset_id);
  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);
  H5Pclose(h_plist_id);
  
}
#endif /* HAVE_HDF5*/
