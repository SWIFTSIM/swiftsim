/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#include <float.h>
#include <limits.h>
#include <math.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "space_unique_id.h"

/* Local headers. */
#include "engine.h"
#include "lock.h"
#include "space.h"

/**
 * @brief Update the unique id structure by requesting a
 * new slice if required.
 *
 * @param s The #space.
 */
void space_update_unique_id(struct space *s) {
  const int require_new_slice = s->unique_id.next.current == 0;

#ifdef WITH_MPI
  const struct engine *e = s->e;

  /* Check if the other ranks need a slices */
  int *all_requires = (int *)malloc(sizeof(int) * e->nr_nodes);

  /* Do the communication */
  MPI_Allgather(&require_new_slice, 1, MPI_INT, all_requires, 1, MPI_INT,
                MPI_COMM_WORLD);

  /* Compute the position of this rank slice and the position of
     the next free slice. */
  int local_index = 0;
  int total_shift = 0;
  for (int i = 0; i < e->nr_nodes; i++) {
    total_shift += 1;
    if (i < engine_rank) {
      local_index += 1;
    }
  }

  /* Free the allocated resources. */
  free(all_requires);

#else

  int local_index = 0;
  int total_shift = require_new_slice;

#endif  // WITH_MPI

  /* Compute the size of the each slice. */
  const long long slice_size = (space_extra_parts + space_extra_sparts +
                                space_extra_gparts + space_extra_bparts) *
                               s->nr_cells;

  /* Get a new slice. */
  if (require_new_slice) {
    /* First check against an overflow. */
    const long long local_shift = local_index * slice_size;
    if (s->unique_id.global_next_id > LLONG_MAX - (local_shift + slice_size)) {
      error("Overflow for the unique IDs.");
    }
    /* Now assign it. */
    s->unique_id.next.current = s->unique_id.global_next_id + local_shift;
    s->unique_id.next.max =
        s->unique_id.global_next_id + local_shift + slice_size;
  }

  /* Shift the position of the next available slice. */
  const long long shift = total_shift * slice_size;
  if (s->unique_id.global_next_id > LLONG_MAX - shift) {
    error("Overflow for the unique IDs.");
  }
  s->unique_id.global_next_id += shift;
}

/**
 * @brief Get a new unique ID.
 *
 * @param s the #space.
 *
 * @return The new unique ID
 */
long long space_get_new_unique_id(struct space *s) {
  /* Get the lock. */
  lock_lock(&s->unique_id.lock);

  /* Get the current available id. */
  const long long id = s->unique_id.current.current;

  /* Update the counter. */
  s->unique_id.current.current++;

  /* Check if everything is fine */
  if (s->unique_id.current.current > s->unique_id.current.max) {
    error("Failed to get a new ID");
  }

  /* Check if need to move to the next slice. */
  if (s->unique_id.current.current == s->unique_id.current.max) {

    /* Check if the next slice is already used */
    if (s->unique_id.next.current == 0) {
      error("Failed to obtain a new unique ID.");
    }
    s->unique_id.current = s->unique_id.next;

    /* Reset the next slice. */
    s->unique_id.next.current = 0;
    s->unique_id.next.max = 0;
  }

  /* Release the lock. */
  if (lock_unlock(&s->unique_id.lock) != 0) {
    error("Impossible to unlock the unique id.");
  }

  return id;
}
