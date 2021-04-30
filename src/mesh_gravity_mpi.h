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

#ifndef SWIFT_MESH_GRAVITY_MPI_H
#define SWIFT_MESH_GRAVITY_MPI_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "hashmap.h"

/* Forward declarations */
struct space;
struct gpart;
struct threadpool;
struct engine;
struct cell;
struct pm_mesh;

void mpi_mesh_accumulate_gparts_to_hashmap(struct threadpool* tp,
                                           const int N, const double fac,
                                           const struct space *s, hashmap_t *map);

void mpi_mesh_hashmaps_to_slices(const int N, const int Nslice, hashmap_t *map,
                                 double *mesh);

void mpi_mesh_fetch_potential(const int N, const double fac,
                              const struct space *s,
                              int local_0_start, int local_n0,
                              double *potential_slice,
                              hashmap_t *potential_map);

void mpi_mesh_update_gparts(struct pm_mesh* mesh, const struct space* s,
			    struct threadpool* tp, const int N, 
			    const double cell_fac);
#endif
