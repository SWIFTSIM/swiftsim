/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "mesh_gravity_mpi.h"

/* Local includes. */
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "part.h"
#include "space.h"
#include "lock.h"
#include "threadpool.h"
#include "active.h"
#include "mesh_gravity_patch.h"
#include "row_major_id.h"
#include "periodic.h"

/**
 * @brief Accumulate contributions from cell to density field
 *
 * Allocates a temporary mesh which covers the top level cell,
 * accumulates mass contributions to this mesh, and then
 * adds these contributions to the supplied hashmap.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param c The #cell containing the particles.
 * @param map The hashmap in which to store the results
 * @param lock A lock used to prevent concurrent access to map
 *
 */
void accumulate_cell_to_hashmap(const int N, const double fac,
				const double *dim, const struct cell *cell,
				hashmap_t *map, swift_lock_type *lock) {

  /* If the cell is empty, then there's nothing to do
     (and the code to find the extent of the cell would fail) */
  if(cell->grav.count == 0)return;

  /* Allocate the local mesh patch */
  struct pm_mesh_patch patch;
  pm_mesh_patch_init(&patch, cell, N, fac, dim, /*boundary_size=*/1);
  pm_mesh_patch_zero(&patch);

  /* Loop over particles in this cell */
  for (int ipart = 0; ipart < cell->grav.count; ipart += 1) {

    const struct gpart *gp = &(cell->grav.parts[ipart]);

    /* Box wrap the particle's position to the copy nearest the cell centre */
    const double pos_x = box_wrap(gp->x[0], patch.wrap_min[0], patch.wrap_max[0]);
    const double pos_y = box_wrap(gp->x[1], patch.wrap_min[1], patch.wrap_max[1]);
    const double pos_z = box_wrap(gp->x[2], patch.wrap_min[2], patch.wrap_max[2]);

    /* Workout the CIC coefficients */
    int i = (int) floor(fac * pos_x);
    const double dx = fac * pos_x - i;
    const double tx = 1. - dx;

    int j = (int) floor(fac * pos_y);
    const double dy = fac * pos_y - j;
    const double ty = 1. - dy;
    
    int k = (int) floor(fac * pos_z);
    const double dz = fac * pos_z - k;
    const double tz = 1. - dz;

    /* Get coordinates within the mesh patch */
    const int ii = i - patch.mesh_min[0];
    const int jj = j - patch.mesh_min[1];
    const int kk = k - patch.mesh_min[2];

    /* Accumulate contributions to the local mesh patch */
    const double mass = gp->mass;
    pm_mesh_patch_CIC_set(&patch, ii, jj, kk, tx, ty, tz, dx, dy, dz, mass);   
  }

  /* Add contributions from the local mesh patch to the hashmap.
   * Need to use a lock here because the hashmap can't be modified
   * by more than one thread at a time. Keeping the lock while we
   * do all our updates seems to be a bit faster than unlocking
   * after each one (in one test case at least) */
  lock_lock(lock);
  pm_mesh_patch_add_values_to_hashmap(&patch, map);
  lock_unlock(lock);
  
  /* Done*/
  pm_mesh_patch_clean(&patch);
}


/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct accumulate_mapper_data {
  const struct cell* cells;
  int N;
  double fac;
  double dim[3];
  hashmap_t *map;
  swift_lock_type *lock;
};


/**
 * @brief Threadpool mapper function for the mesh CIC assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
void accumulate_cell_to_hashmap_mapper(void* map_data, int num, void* extra) {

  /* Unpack the shared information */
  const struct accumulate_mapper_data* data = (struct accumulate_mapper_data*)extra;
  const struct cell* cells = data->cells;
  const int N = data->N;
  const double fac = data->fac;
  const double dim[3] = {data->dim[0], data->dim[1], data->dim[2]};
  hashmap_t *map = data->map;
  swift_lock_type *lock = data->lock;
  
  /* Pointer to the chunk to be processed */
  int* local_cells = (int*)map_data;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell* c = &cells[local_cells[i]];

    /* Assign this cell's content to the mesh */
    accumulate_cell_to_hashmap(N, fac, dim, c, map, lock);
  }
}


/**
 * @brief Accumulate local contributions to the density field
 *
 * Creates a hashmap with the contributions to the density
 * mesh from local particles. Here we require that hash_key_t
 * can store values up to at least N*N*N.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param s The #space containing the particles.
 * @param map The hashmap in which to store the results
 *
 */
void mpi_mesh_accumulate_gparts_to_hashmap(struct threadpool* tp,
                                           const int N, const double fac,
                                           const struct space *s, hashmap_t *map) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Use the threadpool to parallelize over cells */
  struct accumulate_mapper_data data;
  data.cells = s->cells_top;
  data.N = N;
  data.fac = fac;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];
  data.map = map;
  swift_lock_type lock;
  lock_init(&lock);
  data.lock = &lock;
  threadpool_map(tp, accumulate_cell_to_hashmap_mapper, (void*)local_cells,
                 nr_local_cells, sizeof(int), threadpool_auto_chunk_size,
                 (void*)&data);
  lock_destroy(&lock);
  return;
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

/**
 * @brief Store contributions to the mesh as (index, mass) pairs
 */
struct mesh_key_value {
  hashmap_key_t key;
  double value;
};

/**
 * @brief Comparison function to sort mesh_key_value by key
 *
 * @param a The first #mesh_key_value object.
 * @param b The second #mesh_key_value object.
 * @return 1 if a's key field is greater than b's key field,
 * -1 if a's key field is less than b's key field, and zero otherwise
 *
 */
int cmp_func_mesh_key_value(const void *a, const void *b) {
  struct mesh_key_value *a_key_value = (struct mesh_key_value *)a;
  struct mesh_key_value *b_key_value = (struct mesh_key_value *)b;
  if (a_key_value->key > b_key_value->key)
    return 1;
  else if (a_key_value->key < b_key_value->key)
    return -1;
  else
    return 0;
}

/**
 * @brief Data needed to iterate over a hashmap copying key-value pairs
 */
struct hashmap_mapper_data {
  size_t offset;
  struct mesh_key_value *buf;
};

/**
 * @brief Mapper function to copy elements from the hashmap
 *
 * @param key The key associated with this hashmap entry
 * @param value The value associated with this hashmap entry
 * @param data Contains pointer to output buffer and offset to
 * next element to write
 *
 */
void hashmap_copy_elements_mapper(hashmap_key_t key, hashmap_value_t *value,
                                  void *data) {

  struct hashmap_mapper_data *mapper_data = (struct hashmap_mapper_data *)data;
  struct mesh_key_value *element = &(mapper_data->buf[mapper_data->offset]);
  element->key = key;
  element->value = value->value_dbl;
  mapper_data->offset += 1;
}

/**
 * @brief Copy keys and values from a hashmap into an array of struct
 * mesh_key_value, with the elements sorted by key
 *
 * @param map The hashmap to copy into the array
 * @param array The output array
 * @param n The size of the array
 * @return The number of elements copied into the array
 *
 */
size_t hashmap_to_sorted_array(hashmap_t *map, struct mesh_key_value *array, size_t n) {

  /* Find how many elements are in the hashmap */
  size_t num_in_map = hashmap_size(map);
  if(num_in_map > n) {
    error("Array is to small to contain hash map elements!");
  }

  /* Copy the key-value pairs to the new array */
  struct hashmap_mapper_data mapper_data = {(size_t)0, array};
  hashmap_iterate(map, hashmap_copy_elements_mapper, &mapper_data);

  /* And sort the pairs by key */
  qsort(array, num_in_map, sizeof(struct mesh_key_value),
        cmp_func_mesh_key_value);

  return num_in_map;
}


/**
 * @brief Given an array of structs of size element_size, send 
 * nr_send[i] elements to each node i. Allocates the receive
 * buffer recvbuf to the appropriate size and returns its size
 * in nr_recv_tot.
 *
 * TODO: can/should we replace this with a call to engine_do_redistribute()?
 *
 * @param nr_send Number of elements to send to each node
 * @param nr_recv Number of elements to receive from each node
 * @param sendbuf The elements to send
 * @param recvbuf The output buffer
 *
 */
void exchange_structs(size_t *nr_send, char *sendbuf,
                      size_t *nr_recv, char *recvbuf,
                      size_t element_size) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Compute send offsets */
  size_t *send_offset = malloc(nr_nodes * sizeof(size_t));
  send_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    send_offset[i] = send_offset[i - 1] + nr_send[i - 1];
  }

  /* Compute receive offsets */
  size_t *recv_offset = malloc(nr_nodes * sizeof(size_t));
  recv_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    recv_offset[i] = recv_offset[i - 1] + nr_recv[i - 1];
  }

  /* Allocate request objects (one send and receive per node) */
  MPI_Request *request = malloc(2 * sizeof(MPI_Request) * nr_nodes);

  /* Make type to communicate mesh_key_value struct */
  MPI_Datatype mesh_key_value_mpi_type;
  if (MPI_Type_contiguous(element_size, MPI_BYTE, &mesh_key_value_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&mesh_key_value_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for mesh_key_value struct.");
  }

  /*
   * Post the send operations. This is an alltoallv really but
   * we want to avoid the limits imposed by int counts and offsets
   * in MPI_Alltoallv.
   */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_send[i] > 0) {

      /* TODO: handle very large messages */
      if(nr_send[i] > INT_MAX) error("exchange_structs() fails if nr_send > INT_MAX!");

      MPI_Isend(&(sendbuf[send_offset[i]*element_size]), (int)nr_send[i],
                mesh_key_value_mpi_type, i, 0, MPI_COMM_WORLD, &(request[i]));
    } else {
      request[i] = MPI_REQUEST_NULL;
    }
  }

  /* Post the receives */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_recv[i] > 0) {

      /* TODO: handle very large messages */
      if(nr_recv[i] > INT_MAX) error("exchange_structs() fails if nr_recv > INT_MAX!");

      MPI_Irecv(&(recvbuf[recv_offset[i]*element_size]), (int)nr_recv[i],
                mesh_key_value_mpi_type, i, 0, MPI_COMM_WORLD,
                &(request[i + nr_nodes]));
    } else {
      request[i + nr_nodes] = MPI_REQUEST_NULL;
    }
  }

  /* Wait for everything to complete */
  MPI_Waitall(2 * nr_nodes, request, MPI_STATUSES_IGNORE);

  /* Done with the MPI type */
  MPI_Type_free(&mesh_key_value_mpi_type);

  /* Tidy up */
  free(recv_offset);
  free(send_offset);
  free(request);
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}

/**
 * @brief Convert hashmaps to a slab-distributed 3D mesh
 *
 * For FFTW each rank needs to hold a slice of the full mesh.
 * This routine does the necessary communication to convert
 * the per-rank hashmaps into a slab-distributed mesh.
 *
 * @param N The size of the mesh
 * @param local_n0 The thickness of the slice to store on this rank
 * @param map The hashmap with the local part of the mesh
 * @param mesh Pointer to the output data buffer
 *
 */
void mpi_mesh_hashmaps_to_slices(const int N, const int local_n0, hashmap_t *map,
                                 double *mesh) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Make an array with the (key, value) pairs from the hashmap.
   * The elements are sorted by key, which means they're sorted 
   * by x coordinate, then y coordinate, then z coordinate.
   * We're going to distribute them between ranks according to their
   * x coordinate, so this puts them in order of destination rank.
   */
  size_t map_size = hashmap_size(map);
  struct mesh_key_value *mesh_sendbuf;
  if (swift_memalign("mesh_sendbuf", (void **) &mesh_sendbuf, 32,
                     map_size * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate array for mesh send buffer!");
  size_t nr_send_tot = hashmap_to_sorted_array(map, mesh_sendbuf, map_size);

  /* Get width of the slice on each rank */
  int *slice_width = malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&local_n0, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine offset to the slice on each rank */
  int *slice_offset = malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Compute how many elements are to be sent to each rank */
  size_t *nr_send = malloc(nr_nodes * sizeof(size_t));
  for (int i = 0; i < nr_nodes; i += 1) {
    nr_send[i] = 0;
  }
  int dest_node = 0;
  for (size_t i = 0; i < nr_send_tot; i += 1) {
    /* Get the x coordinate of this mesh cell in the global mesh */
    int mesh_x = get_xcoord_from_padded_row_major_id((size_t) mesh_sendbuf[i].key, N);
    /* Advance to the destination node that is to contain this x coordinate */
    while ((mesh_x >= slice_offset[dest_node] + slice_width[dest_node]) ||
           (slice_width[dest_node] == 0)) {
      dest_node += 1;
    }
    nr_send[dest_node] += 1;
  }

  /* Determine how many requests we'll receive from each MPI rank */
  size_t *nr_recv = malloc(sizeof(size_t)*nr_nodes);
  MPI_Alltoall(nr_send, sizeof(size_t), MPI_BYTE, nr_recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t nr_recv_tot = 0;
  for(int i=0; i<nr_nodes; i+=1) {
    nr_recv_tot += nr_recv[i];
  }

  /* Allocate the receive buffer */
  struct mesh_key_value *mesh_recvbuf;
  if (swift_memalign("mesh_recvbuf", (void **) &mesh_recvbuf, 32,
                     nr_recv_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate receive buffer for constructing MPI FFT mesh");

  /* Carry out the communication */
  exchange_structs(nr_send, (char *) mesh_sendbuf,
                   nr_recv, (char *) mesh_recvbuf, 
                   sizeof(struct mesh_key_value));

  /* Copy received data to the output buffer */
  for (size_t i = 0; i < nr_recv_tot; i += 1) {
#ifdef SWIFT_DEBUG_CHECKS
    const int xcoord = get_xcoord_from_padded_row_major_id(mesh_recvbuf[i].key, N);
    if(xcoord < slice_offset[nodeID])
      error("Received mesh cell is not in the local slice (xcoord too small)");
    if(xcoord >= slice_offset[nodeID] + slice_width[nodeID])
      error("Received mesh cell is not in the local slice (xcoord too large)");
#endif
    mesh[get_index_in_local_slice((size_t) mesh_recvbuf[i].key, N, slice_offset[nodeID])] += mesh_recvbuf[i].value;
  }

  /* Tidy up */
  free(slice_width);
  free(slice_offset);
  free(nr_send);
  free(nr_recv);
  swift_free("mesh_recvbuf", mesh_recvbuf);
  swift_free("mesh_sendbuf", mesh_sendbuf);
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}


/**
 * @brief Retrieve the potential in the mesh cells we need to
 * compute the force on particles on this MPI rank. Result is
 * returned in the supplied hashmap, which should be initially
 * empty.
 *
 * We need all cells containing points -2 and +3 mesh cell widths
 * away from each particle along each axis to compute the
 * potential gradient.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the FFT mesh cell size
 * @param s The #space containing the particles.
 * @param local_0_start Offset to the first mesh x coordinate on this rank
 * @param local_n0 Width of the mesh slab on this rank
 * @param potential_slice Array with the potential on the local slice of the mesh
 * @param potential_map A hashmap in which to store the potential data
 *
 */
void mpi_mesh_fetch_potential(const int N, const double fac,
                              const struct space *s,
                              int local_0_start, int local_n0,
                              double *potential_slice,
                              hashmap_t *potential_map) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Determine rank, number of MPI ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Create a hashmap to store the indexes of the cells we need */
  hashmap_t map;
  hashmap_init(&map);

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell += 1) {
    struct cell *cell = &(s->cells_top[local_cells[icell]]);
    if(cell->grav.count > 0) {

      /* Determine range of FFT mesh cells we need for particles in this top
         level cell. The 5 point stencil used for accelerations requires
         2 neighbouring FFT mesh cells in each direction and for CIC
         evaluation of the accelerations we need one extra FFT mesh cell
         in the +ve direction.

	 We also have to add a small buffer to avoid problems with rounding
	 (e.g. if we decide a cell isn't needed here but then try to look
	 up its value in the hashmap later).

	 TODO: can we calculate exactly how big the rounding error can be?
	       Will just assume that 1% of a mesh cell is enough for now.
      */
      int ixmin[3];
      int ixmax[3];
      for(int idim=0;idim<3;idim+=1) {
        const double xmin = cell->loc[idim] - 2.01/fac;
        const double xmax = cell->loc[idim] + cell->width[idim] + 3.01/fac;
        ixmin[idim] = (int) floor(xmin*fac);
        ixmax[idim] = (int) floor(xmax*fac);
      }

      /* Add the required cells to the map */
      for(int i=ixmin[0]; i<=ixmax[0]; i+=1) {
        for(int j=ixmin[1]; j<=ixmax[1]; j+=1) {
          for(int k=ixmin[2]; k<=ixmax[2]; k+=1) {
            const size_t index = row_major_id_periodic_size_t_padded(i, j, k, N);
            /* We don't have a value associated with the entry yet */
	    hashmap_value_t *value = hashmap_get(&map, (hashmap_key_t) index);
	    value->value_dbl = 0.0;
          }
        }
      }
    }
  }

  /* Make an array with the cell IDs we need to request from other ranks */
  size_t nr_send_tot = hashmap_size(&map);
  struct mesh_key_value *send_cells;
  if (swift_memalign("send_cells", (void **) &send_cells, 32,
                     nr_send_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate array for cells to request!");
  nr_send_tot = hashmap_to_sorted_array(&map, send_cells, nr_send_tot);
  hashmap_free(&map);

  /* Get width of the mesh slice on each rank */
  int *slice_width = malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&local_n0, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine first mesh x coordinate stored on each rank */
  int *slice_offset = malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Count how many mesh cells we need to request from each MPI rank */
  int dest_rank = 0;
  size_t *nr_send = malloc(sizeof(size_t)*nr_nodes);
  for(int i=0; i<nr_nodes; i+=1) {
    nr_send[i] = 0;
  }
  for(size_t i=0; i<nr_send_tot; i+=1) {
    while(get_xcoord_from_padded_row_major_id(send_cells[i].key, N) >= (slice_offset[dest_rank]+slice_width[dest_rank]) || slice_width[dest_rank] == 0) {
      dest_rank += 1;
    }
#ifdef SWIFT_DEBUG_CHECKS
    if(dest_rank >= nr_nodes || dest_rank < 0)error("Destination rank out of range");
#endif
    nr_send[dest_rank] += 1;
  }

  /* Determine how many requests we'll receive from each MPI rank */
  size_t *nr_recv = malloc(sizeof(size_t)*nr_nodes);
  MPI_Alltoall(nr_send, sizeof(size_t), MPI_BYTE, nr_recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t nr_recv_tot = 0;
  for(int i=0; i<nr_nodes; i+=1) {
    nr_recv_tot += nr_recv[i];
  }

  /* Allocate buffer to receive requests */
  struct mesh_key_value *recv_cells;
  if (swift_memalign("recv_cells", (void **) &recv_cells, 32,
                     nr_recv_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate array for mesh receive buffer!");

  /* Send requests for cells to other ranks */
  exchange_structs(nr_send, (char *) send_cells,
                   nr_recv, (char *) recv_cells,
                   sizeof(struct mesh_key_value));
  
  /* Look up potential in the requested cells */
  for(size_t i=0; i<nr_recv_tot; i+=1) {
#ifdef SWIFT_DEBUG_CHECKS
    const size_t cells_in_slab  = ((size_t) N)*(2*(N/2+1));
    const size_t first_local_id = local_0_start*cells_in_slab;
    const size_t num_local_ids  = local_n0*cells_in_slab;
    if(recv_cells[i].key < first_local_id || recv_cells[i].key >= first_local_id+num_local_ids) {
      error("Requested potential mesh cell ID is out of range");
    }
#endif
    size_t local_id = get_index_in_local_slice(recv_cells[i].key, N, local_0_start);
#ifdef SWIFT_DEBUG_CHECKS
    const size_t Ns = N;
    if(local_id >= Ns*(2*(Ns/2+1))*local_n0)error("Local potential mesh cell ID is out of range");
#endif
    recv_cells[i].value = potential_slice[local_id];
  }

  /* Return the results */
  exchange_structs(nr_recv, (char *) recv_cells,
                   nr_send, (char *) send_cells,
                   sizeof(struct mesh_key_value));

  /* Store the results in the hashmap */
  for(size_t i=0; i<nr_send_tot; i+=1) {
    int created = 0;
    hashmap_value_t *value = hashmap_get_new(potential_map, send_cells[i].key,
					     &created);
    value->value_dbl = send_cells[i].value;
#ifdef SWIFT_DEBUG_CHECKS
    if(!created)error("Received duplicate potential hash map value");
    const size_t Ns = N;
    if(send_cells[i].key >= Ns*Ns*(2*(Ns/2+1)))error("Received potential mesh cell ID out of range");
#endif
  }

  /* Tidy up */
  swift_free("send_cells", send_cells);
  swift_free("recv_cells", recv_cells);
  free(slice_width);
  free(slice_offset);
  free(nr_send);
  free(nr_recv);
#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}


/**
 * @brief Computes the potential on a gpart from a given mesh using the CIC
 * method.
 *
 * @param gp The #gpart.
 * @param patch The local mesh patch
 */
#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)
void mesh_patch_to_gparts_CIC(struct gpart *gp, const struct pm_mesh_patch *patch) {

  const double fac = patch->fac;

  /* Box wrap the gpart's position to the copy nearest the cell centre */
  const double pos_x = box_wrap(gp->x[0], patch->wrap_min[0], patch->wrap_max[0]);
  const double pos_y = box_wrap(gp->x[1], patch->wrap_min[1], patch->wrap_max[1]);
  const double pos_z = box_wrap(gp->x[2], patch->wrap_min[2], patch->wrap_max[2]);

  /* Workout the CIC coefficients */
  int i = (int) floor(fac * pos_x);
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int) floor(fac * pos_y);
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;
    
  int k = (int) floor(fac * pos_z);
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (gp->a_grav_mesh[0] != 0.) error("Particle with non-initalised stuff");
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  if (gp->potential_mesh != 0.) error("Particle with non-initalised stuff");
#endif
#endif

  /* Some local accumulators */
  double p = 0.;
  double a[3] = {0.};

  /* Get coordinates within the mesh patch */
  const int ii = i - patch->mesh_min[0];
  const int jj = j - patch->mesh_min[1];
  const int kk = k - patch->mesh_min[2];

  /* Simple CIC for the potential itself */
  p += pm_mesh_patch_CIC_get(patch, ii, jj, kk, tx, ty, tz, dx, dy, dz);

  /* 5-point stencil along each axis for the accelerations */
  a[0] += (1. / 12.) * pm_mesh_patch_CIC_get(patch, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (2. / 3.)  * pm_mesh_patch_CIC_get(patch, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] += (2. / 3.)  * pm_mesh_patch_CIC_get(patch, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  a[0] -= (1. / 12.) * pm_mesh_patch_CIC_get(patch, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);
  
  a[1] += (1. / 12.) * pm_mesh_patch_CIC_get(patch, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (2. / 3.)  * pm_mesh_patch_CIC_get(patch, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] += (2. / 3.)  * pm_mesh_patch_CIC_get(patch, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  a[1] -= (1. / 12.) * pm_mesh_patch_CIC_get(patch, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);
  
  a[2] += (1. / 12.) * pm_mesh_patch_CIC_get(patch, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  a[2] -= (2. / 3.)  * pm_mesh_patch_CIC_get(patch, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  a[2] += (2. / 3.)  * pm_mesh_patch_CIC_get(patch, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  a[2] -= (1. / 12.) * pm_mesh_patch_CIC_get(patch, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);

  /* Store things back */
  gp->a_grav_mesh[0] = fac * a[0];
  gp->a_grav_mesh[1] = fac * a[1];
  gp->a_grav_mesh[2] = fac * a[2];
  gravity_add_comoving_mesh_potential(gp, p);
}
#endif


/**
 * @brief Interpolate the forces and potential from the mesh to the #gpart.
 *
 * This is for the case where the mesh is distributed between MPI ranks
 * and stored in the form of a hashmap. This function updates the particles
 * in one #cell.
 *
 * @param c The #cell containing the #gpart to update
 * @param potential Hashmap containing the potential to interpolate from.
 * @param N Size of the full mesh
 * @param fac Inverse of the FFT mesh cell size
 * @param const_G Gravitional constant
 * @param dim Dimensions of the #space
 */
void cell_distributed_mesh_to_gpart_CIC(const struct cell *c, hashmap_t *potential, 
					const int N, const double fac, const float const_G,
					const double dim[3]) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  const int gcount = c->grav.count;
  struct gpart* gparts = c->grav.parts;

  /* Check for empty cell as this would cause problems finding the extent */
  if(gcount==0)return;

  /* Allocate the local mesh patch */
  struct pm_mesh_patch patch;
  pm_mesh_patch_init(&patch, c, N, fac, dim, /*boundary_size=*/2);
  
  /* Populate the mesh patch with values from the potential hashmap */
  pm_mesh_patch_set_values_from_hashmap(&patch, potential);

  /* Get the potential from the mesh patch to the active gparts using CIC */
  for (int i = 0; i < gcount; ++i) {
    struct gpart* gp = &gparts[i];

    if (gp->time_bin == time_bin_inhibited) continue;

    gp->a_grav_mesh[0] = 0.f;
    gp->a_grav_mesh[1] = 0.f;
    gp->a_grav_mesh[2] = 0.f;
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
    gp->potential_mesh = 0.f;
#endif

    mesh_patch_to_gparts_CIC(gp, &patch);

    gp->a_grav_mesh[0] *= const_G;
    gp->a_grav_mesh[1] *= const_G;
    gp->a_grav_mesh[2] *= const_G;
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
    gp->potential_mesh *= const_G;
#endif
  }

  /* Discard the temporary mesh patch */
  pm_mesh_patch_clean(&patch);

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}


/**
 * @brief Shared information about the mesh to be used by all the threads in the
 * pool.
 */
struct distributed_cic_mapper_data {
  const struct cell* cells;
  hashmap_t* potential;
  int N;
  double fac;
  double dim[3];
  float const_G;
};

/**
 * @brief Threadpool mapper function for the mesh CIC assignment of a cell.
 *
 * @param map_data A chunk of the list of local cells.
 * @param num The number of cells in the chunk.
 * @param extra The information about the mesh and cells.
 */
void cell_distributed_mesh_to_gpart_CIC_mapper(void* map_data, int num, void* extra) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

  /* Unpack the shared information */
  const struct distributed_cic_mapper_data* data = (struct distributed_cic_mapper_data*)extra;
  const struct cell* cells = data->cells;
  hashmap_t* potential = data->potential;
  const int N = data->N;
  const double fac = data->fac;
  const double dim[3] = {data->dim[0], data->dim[1], data->dim[2]};
  const float const_G = data->const_G;

  /* Pointer to the chunk to be processed */
  int* local_cells = (int*)map_data;

  /* Loop over the elements assigned to this thread */
  for (int i = 0; i < num; ++i) {

    /* Pointer to local cell */
    const struct cell* c = &cells[local_cells[i]];

    /* Update acceleration and potential for gparts in this cell */
    cell_distributed_mesh_to_gpart_CIC(c, potential, N, fac, const_G, dim);
  }

#else
  error("FFTW MPI not found - unable to use distributed mesh");
#endif
}



 void mpi_mesh_update_gparts(struct pm_mesh* mesh, const struct space* s,
			     struct threadpool* tp, const int N, 
			     const double cell_fac) {

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

   const int* local_cells = s->local_cells_top;
   const int nr_local_cells = s->nr_local_cells;

   /* Gather the mesh shared information to be used by the threads */
   struct distributed_cic_mapper_data data;
   data.cells = s->cells_top;
   data.potential = mesh->potential_local;
   data.N = N;
   data.fac = cell_fac;
   data.dim[0] = s->dim[0];
   data.dim[1] = s->dim[1];
   data.dim[2] = s->dim[2];
   data.const_G = s->e->physical_constants->const_newton_G;

   if (nr_local_cells == 0) {
     error("Distributed mesh not implemented without cells");
   } else {
     /* Evaluate acceleration and potential for each gpart */
     threadpool_map(tp, cell_distributed_mesh_to_gpart_CIC_mapper, 
		    (void*)local_cells, nr_local_cells, sizeof(int),
		    threadpool_auto_chunk_size, (void*)&data);
   }
#else
   error("FFTW MPI not found - unable to use distributed mesh");
#endif
 }
