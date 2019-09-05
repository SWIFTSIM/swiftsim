/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "memswap.h"

/* Load the profiler header, if needed. */
#ifdef WITH_PROFILER
#include <gperftools/profiler.h>
#endif

/**
 * @brief Unskip any hydro tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_hydro(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->hydro.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_hydro(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_hydro_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any stars tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param with_star_formation Are we running with star formation switched on?
 */
static void engine_do_unskip_stars(struct cell *c, struct engine *e,
                                   const int with_star_formation) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  const int non_empty =
      c->stars.count > 0 || (with_star_formation && c->hydro.count > 0);

  /* Ignore empty cells. */
  if (!non_empty) return;

  const int ci_active = cell_is_active_stars(c, e) ||
                        (with_star_formation && cell_is_active_hydro(c, e));

  /* Skip inactive cells. */
  if (!ci_active) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_stars(cp, e, with_star_formation);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild =
      cell_unskip_stars_tasks(c, &e->sched, with_star_formation);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any black hole tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_black_holes(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->black_holes.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_black_holes(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_black_holes(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_black_holes_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any gravity tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_gravity(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->grav.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse */
  if (c->split && ((c->maxdepth - c->depth) >= space_subdepth_diff_grav)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_gravity(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  cell_unskip_gravity_tasks(c, &e->sched);
}

/**
 * @brief Mapper function to unskip active tasks.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  const int with_star_formation = e->policy & engine_policy_star_formation;
  const int nodeID = e->nodeID;
  struct space *s = e->s;
  int *local_cells = (int *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];
    if (c != NULL) {

      /* Hydro tasks */
      if (e->policy & engine_policy_hydro) engine_do_unskip_hydro(c, e);

      /* All gravity tasks */
      if ((e->policy & engine_policy_self_gravity) ||
          ((e->policy & engine_policy_external_gravity) && c->nodeID == nodeID))
        engine_do_unskip_gravity(c, e);

      /* Stars tasks */
      if (e->policy & engine_policy_stars)
        engine_do_unskip_stars(c, e, with_star_formation);

      /* Black hole tasks */
      if (e->policy & engine_policy_black_holes)
        engine_do_unskip_black_holes(c, e);
    }
  }
}

/**
 * @brief Unskip all the tasks that act on active cells at this time.
 *
 * @param e The #engine.
 */
void engine_unskip(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  const int with_hydro = e->policy & engine_policy_hydro;
  const int with_self_grav = e->policy & engine_policy_self_gravity;
  const int with_ext_grav = e->policy & engine_policy_external_gravity;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_feedback = e->policy & engine_policy_feedback;
  const int with_black_holes = e->policy & engine_policy_black_holes;

#ifdef WITH_PROFILER
  static int count = 0;
  char filename[100];
  sprintf(filename, "/tmp/swift_engine_do_usnkip_mapper_%06i.prof", count++);
  ProfilerStart(filename);
#endif  // WITH_PROFILER

  /* Move the active local cells to the top of the list. */
  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if ((with_hydro && cell_is_active_hydro(c, e)) ||
        (with_self_grav && cell_is_active_gravity(c, e)) ||
        (with_ext_grav && c->nodeID == nodeID &&
         cell_is_active_gravity(c, e)) ||
        (with_feedback && cell_is_active_stars(c, e)) ||
        (with_stars && c->nodeID == nodeID && cell_is_active_stars(c, e)) ||
        (with_black_holes && cell_is_active_black_holes(c, e))) {

      if (num_active_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_cells], sizeof(int));
      num_active_cells += 1;
    }
  }

  /* Should we duplicate the list of active cells to better parallelise the
     unskip over the threads ? */
  int multiplier = 0;
  multiplier += (with_hydro > 0);
  multiplier += (with_self_grav > 0 || with_ext_grav > 0);
  multiplier += (with_feedback > 0 || with_stars > 0);
  multiplier += (with_black_holes > 0);

  int *local_active_cells;
  if (multiplier > 1) {

    /* Make space for copies of the list */
    local_active_cells = (int *)malloc(multiplier * num_active_cells);
    if (local_active_cells == NULL)
      error(
          "Couldn't allocate memory for duplicated list of local active "
          "cells.");

    /* Make blind copies of the list */
    for (int m = 0; m < multiplier; m++) {
      memcpy(local_active_cells + m * num_active_cells, local_cells,
             num_active_cells * sizeof(int));
    }
  } else {
    local_active_cells = local_cells;
  }

  /* We now have a list of local active cells duplicated as many times as
   * we have broad task types. We can now release all the threads on the list */

  /* Activate all the regular tasks */
  threadpool_map(&e->threadpool, engine_do_unskip_mapper, local_active_cells,
                 num_active_cells * multiplier, sizeof(int), 0, e);

#ifdef WITH_PROFILER
  ProfilerStop();
#endif  // WITH_PROFILER

  /* Free stuff? */
  if (multiplier > 1) {
    free(local_active_cells);
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
