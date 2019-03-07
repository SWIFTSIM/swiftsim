/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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

#ifndef SWIFT_FOF_H
#define SWIFT_FOF_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "align.h"
#include "parser.h"

/* Avoid cyclic inclusions */
struct space;

struct fof_CoM {

  double x, y, z;

} SWIFT_STRUCT_ALIGN;

/* MPI message required for FOF. */
struct fof_mpi {

  /* The local particle's root ID.*/
  size_t group_i;

  /* The local group's size.*/
  size_t group_i_size;

  /* The local group's mass.*/
  double group_i_mass;

  /* The local group's CoM.*/
  struct fof_CoM group_i_CoM;

  /* The foreign particle's root ID.*/
  size_t group_j;

  /* The local group's size.*/
  size_t group_j_size;

  /* The local group's mass.*/
  double group_j_mass;

  /* The local group's CoM.*/
  struct fof_CoM group_j_CoM;

} SWIFT_STRUCT_ALIGN;

/* MPI message used for linking FoF fragments between nodes */
struct fof_mpi_links {
  /* Index of the local group */
  size_t group_i;
  /* Index of the remote group */
  size_t group_j;
  /* Minimum group index associated with this link */
  size_t group_index_min;
} SWIFT_STRUCT_ALIGN;

/* MPI message used for updating FoF sizes etc between nodes */
struct fof_mpi_sizes {
  /* Index of the local group */
  size_t group_i;
  /* Index of the remote group */
  size_t group_j;
  /* Number of particles in the local group */
  size_t size;
} SWIFT_STRUCT_ALIGN;


struct fof {

  size_t *group_index;
  size_t *group_size;
  double *group_mass;
  struct fof_CoM *group_CoM;
  int num_groups;
  size_t min_group_size;
  size_t group_id_default;
  size_t group_id_offset;
  int group_links_size_default;

  int group_link_count;
  struct fof_mpi *group_links;
  int group_links_size;

  char base_name[PARSER_MAX_LINE_SIZE];

} SWIFT_STRUCT_ALIGN;

/* Store group size and offset into array. */
struct group_length {

  size_t index, size;

} SWIFT_STRUCT_ALIGN;

/* Store local and foreign cell indices that touch. */
struct cell_pair_indices {

  struct cell *local, *foreign;

} SWIFT_STRUCT_ALIGN;

/* Function prototypes. */
void fof_init(struct space *s);
void fof_search_cell(struct space *s, struct cell *c);
void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj);
void fof_search_pair_cells_foreign(struct space *s, struct cell *ci,
                                   struct cell *cj, int *link_count,
                                   struct fof_mpi **group_links,
                                   int *group_links_size);
void fof_search_tree(struct space *s);
void fof_dump_group_data(char *out_file, struct space *s,
                         struct group_length *group_sizes);
void rec_fof_search_self(struct cell *c, struct space *s, const double dim[3],
                         const double search_r2);
void rec_fof_search_pair(struct cell *restrict ci, struct cell *restrict cj,
                         struct space *s, const double dim[3],
                         const double search_r2);

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype fof_mpi_type;
extern MPI_Datatype group_length_mpi_type;
#endif

#endif /* SWIFT_FOF_H */
