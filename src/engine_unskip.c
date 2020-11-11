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
 * @brief Broad categories of tasks.
 *
 * Each category is unskipped independently
 * of the others.
 */
enum task_broad_types {
  task_broad_types_hydro = 1,
  task_broad_types_gravity,
  task_broad_types_stars,
  task_broad_types_sinks,
  task_broad_types_black_holes,
  task_broad_types_rt,
  task_broad_types_count,
};

/**
 * @brief Meta-data for the unskipping
 */
struct unskip_data {

  /*! The #engine */
  struct engine *e;

  /*! Pointer to the start of the list of cells to unskip */
  struct cell **list_base;

  /*! Number of times the list has been duplicated */
  int multiplier;

  /*! The number of active cells (without dulication) */
  int num_active_cells;

  /*! The #task_broad_types corresponding to each copy of the list */
  enum task_broad_types task_types[task_broad_types_count];
};

enum recursion_type {
    recursion_type_above,
    recursion_type_below,
    recursion_type_both,
};

/**
 * @brief Unskip any hydro tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param recursion The recursion type (only do above, below or both of them).
 */
static void engine_do_unskip_hydro(struct cell *c, struct engine *e,
                                   enum recursion_type recursion) {

  /* Unskip the parents of the reference cell. */
  if (recursion == recursion_type_both && c->parent != NULL) {
    struct cell *parent = c->parent;
    /* Do the cells above only if c is the first progeny. */
    for(int i = 0; i < 8; i++) {
      if (parent->progeny[i] == c) {
        engine_do_unskip_hydro(parent, e, recursion_type_above);
        break;
      }
      if (parent->progeny[i] != NULL) {
        break;
      }
    }
  }

  /* We only need to climb the three to the top level cell. */
  if (recursion == recursion_type_above) {
    if (c->parent != NULL) {
      engine_do_unskip_hydro(c->parent, e, recursion_type_above);
    }
  }

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->hydro.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_hydro(c, e)) return;

  /* Do the recursion below the current level. */
  if (c->split && (recursion == recursion_type_below ||
                   recursion == recursion_type_both)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_hydro(cp, e, recursion_type_below);
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
 * @param recursion The recursion type (only do above, below or both of them).
 */
static void engine_do_unskip_stars(struct cell *c, struct engine *e,
                                   const int with_star_formation, enum recursion_type recursion) {

  /* Unskip the parents of the reference cell. */
  if (recursion == recursion_type_both && c->parent != NULL) {
    struct cell *parent = c->parent;
    /* Do the cells above only if c is the first progeny. */
    for(int i = 0; i < 8; i++) {
      if (parent->progeny[i] == c) {
        engine_do_unskip_stars(parent, e, with_star_formation, recursion_type_above);
        break;
      }
      if (parent->progeny[i] != NULL) {
        break;
      }
    }
  }

  /* We only need to climb the three to the top level cell. */
  if (recursion == recursion_type_above) {
    if (c->parent != NULL) {
      engine_do_unskip_stars(c->parent, e, with_star_formation, recursion_type_above);
    }
  }

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

  /* Do the recursion below the current level. */
  if (c->split && (recursion == recursion_type_below ||
                   recursion == recursion_type_both)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_stars(cp, e, with_star_formation, recursion_type_below);
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
 * @param recursion The recursion type (only do above, below or both of them).
 */
static void engine_do_unskip_black_holes(struct cell *c, struct engine *e,
                                         enum recursion_type recursion) {

  /* Unskip the parents of the reference cell. */
  if (recursion == recursion_type_both && c->parent != NULL) {
    struct cell *parent = c->parent;
    /* Do the cells above only if c is the first progeny. */
    for(int i = 0; i < 8; i++) {
      if (parent->progeny[i] == c) {
        engine_do_unskip_black_holes(parent, e, recursion_type_above);
        break;
      }
      if (parent->progeny[i] != NULL) {
        break;
      }
    }
  }

  /* We only need to climb the three to the top level cell. */
  if (recursion == recursion_type_above) {
    if (c->parent != NULL) {
      engine_do_unskip_black_holes(c->parent, e, recursion_type_above);
    }
  }

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->black_holes.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_black_holes(c, e)) return;

  /* Do the recursion below the current level. */
  if (c->split && (recursion == recursion_type_below ||
                   recursion == recursion_type_both)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_black_holes(cp, e, recursion_type_below);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_black_holes_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any sink tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param recursion The recursion type (only do above, below or both of them).
 */
static void engine_do_unskip_sinks(struct cell *c, struct engine *e,
                                   enum recursion_type recursion) {

  /* Unskip the parents of the reference cell. */
  if (recursion == recursion_type_both && c->parent != NULL) {
    struct cell *parent = c->parent;
    /* Do the cells above only if c is the first progeny. */
    for(int i = 0; i < 8; i++) {
      if (parent->progeny[i] == c) {
        engine_do_unskip_sinks(parent, e, recursion_type_above);
        break;
      }
      if (parent->progeny[i] != NULL) {
        break;
      }
    }
  }

  /* We only need to climb the three to the top level cell. */
  if (recursion == recursion_type_above) {
    if (c->parent != NULL) {
      engine_do_unskip_sinks(c->parent, e, recursion_type_above);
    }
  }

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->sinks.count == 0 && c->hydro.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_sinks(c, e) && !cell_is_active_hydro(c, e)) return;

  /* Do the recursion below the current level. */
  if (c->split && (recursion == recursion_type_below ||
                   recursion == recursion_type_both)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_sinks(cp, e, recursion_type_below);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_sinks_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any gravity tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param recursion The recursion type (only do above, below or both of them).
 */
static void engine_do_unskip_gravity(struct cell *c, struct engine *e,
                                     enum recursion_type recursion) {

  /* Unskip the parents of the reference cell. */
  if (recursion == recursion_type_both && c->parent != NULL) {
    struct cell *parent = c->parent;
    /* Do the cells above only if c is the first progeny. */
    for(int i = 0; i < 8; i++) {
      if (parent->progeny[i] == c) {
        engine_do_unskip_gravity(parent, e, recursion_type_above);
        break;
      }
      if (parent->progeny[i] != NULL) {
        break;
      }
    }
  }

  /* We only need to climb the three to the top level cell. */
  if (recursion == recursion_type_above) {
    if (c->parent != NULL) {
      engine_do_unskip_gravity(c->parent, e, recursion_type_above);
    }
  }

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->grav.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_gravity(c, e)) return;

  /* Do the recursion below the current level. */
  if (c->split && ((c->maxdepth - c->depth) >= space_subdepth_diff_grav)
      && (recursion == recursion_type_below ||
          recursion == recursion_type_both)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_gravity(cp, e, recursion_type_below);
      }
    }
  }

  /* Unskip any active tasks. */
  cell_unskip_gravity_tasks(c, &e->sched);
}

/**
 * @brief Unskip any radiative transfer tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param recursion The recursion type (only do above, below or both of them).
 */
static void engine_do_unskip_rt(struct cell *c, struct engine *e,
                                enum recursion_type recursion) {

  /* Unskip the parents of the reference cell. */
  if (recursion == recursion_type_both && c->parent != NULL) {
    struct cell *parent = c->parent;
    /* Do the cells above only if c is the first progeny. */
    for(int i = 0; i < 8; i++) {
      if (parent->progeny[i] == c) {
        engine_do_unskip_rt(parent, e, recursion_type_above);
        break;
      }
      if (parent->progeny[i] != NULL) {
        break;
      }
    }
  }

  /* We only need to climb the three to the top level cell. */
  if (recursion == recursion_type_above) {
    if (c->parent != NULL) {
      engine_do_unskip_rt(c->parent, e, recursion_type_above);
    }
  }

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->hydro.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_hydro(c, e)) return;

  /* Do the recursion below the current level. */
  if (c->split
      && (recursion == recursion_type_below ||
          recursion == recursion_type_both)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_rt(cp, e, recursion_type_below);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_rt_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Mapper function to unskip active tasks.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an unskip_data structure.
 */
void engine_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  /* Unpack the meta data */
  struct unskip_data *data = (struct unskip_data *)extra_data;
  const int num_active_cells = data->num_active_cells;
  const enum task_broad_types *const task_types = data->task_types;
  struct cell **list_base = data->list_base;
  struct engine *e = data->e;

  /* What policies are we running? */
  const int with_star_formation = e->policy & engine_policy_star_formation;

  /* The current chunk of active cells */
  struct cell **const local_cells = (struct cell **)map_data;

  /* Loop over this thread's chunk of cells to unskip */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Handle on the cell */
    struct cell *const c = local_cells[ind];

    /* In what copy of the global list are we?
     * This gives us the broad type of task we are working on. */
    const ptrdiff_t delta = &local_cells[ind] - list_base;
    const int type = delta / num_active_cells;

#ifdef SWIFT_DEBUG_CHECKS
    if (type >= data->multiplier) error("Invalid broad task type!");
    if (c == NULL) error("Got an invalid cell index!");
#endif

    /* What broad type of tasks are we unskipping? */
    switch (task_types[type]) {
      case task_broad_types_hydro:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_hydro))
          error("Trying to unskip hydro tasks in a non-hydro run!");
#endif
        engine_do_unskip_hydro(c, e, recursion_type_both);
        break;
      case task_broad_types_gravity:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_self_gravity) &&
            !(e->policy & engine_policy_external_gravity))
          error("Trying to unskip gravity tasks in a non-gravity run!");
#endif
        engine_do_unskip_gravity(c, e, recursion_type_both);
        break;
      case task_broad_types_stars:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_stars))
          error("Trying to unskip star tasks in a non-stars run!");
#endif
        engine_do_unskip_stars(c, e, with_star_formation, recursion_type_both);
        break;
      case task_broad_types_sinks:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_sinks))
          error("Trying to unskip sink tasks in a non-sinks run!");
#endif
        engine_do_unskip_sinks(c, e, recursion_type_both);
        break;
      case task_broad_types_black_holes:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_black_holes))
          error("Trying to unskip black holes tasks in a non-BH run!");
#endif
        engine_do_unskip_black_holes(c, e, recursion_type_both);
        break;
      case task_broad_types_rt:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_rt))
          error("Trying to unskip radiative transfer tasks in a non-rt run!");
#endif
        engine_do_unskip_rt(c, e, recursion_type_both);
        break;
      default:
#ifdef SWIFT_DEBUG_CHECKS
        error("Invalid broad task type!");
#endif
        continue;
    }
  }
}

/**
 * @brief Copy the address of the cells that requires an unskip.
 *
 * @param c The current #cell.
 * @param local_active_cells The output array containing the address of the cells.
 * @param level Level at which the unskip is run (0 for top level).
 * @param e The #engine.
 */
struct cell **engine_unskip_copy_cell(struct cell *c, struct cell **local_active_cells, const int level, const struct engine *e) {

  /* Should we go lower? */
  if (c->split && c->depth != level) {
    for(int i = 0; i < 8; i++) {
      if (c->progeny[i] != NULL) {
        local_active_cells = engine_unskip_copy_cell(c->progeny[i], local_active_cells, level, e);
      }
    }
  }
  /* Correct cell found. */
  else {
    local_active_cells[0] = c;
    local_active_cells++;
  }

  return local_active_cells;
}

/**
 * @brief Unskip all the tasks that act on active cells at this time.
 *
 * @param e The #engine.
 */
void engine_unskip(struct engine *e) {

  const int unskip_level = e->unskip_level;

  const ticks tic = getticks();
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  const int with_hydro = e->policy & engine_policy_hydro;
  const int with_self_grav = e->policy & engine_policy_self_gravity;
  const int with_ext_grav = e->policy & engine_policy_external_gravity;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_sinks = e->policy & engine_policy_sinks;
  const int with_feedback = e->policy & engine_policy_feedback;
  const int with_black_holes = e->policy & engine_policy_black_holes;
  const int with_rt = e->policy & engine_policy_rt;

#ifdef WITH_PROFILER
  static int count = 0;
  char filename[100];
  sprintf(filename, "/tmp/swift_engine_do_usnkip_mapper_%06i.prof", count++);
  ProfilerStart(filename);
#endif  // WITH_PROFILER

  /* Move the active local cells to the top of the list. */
  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_top_cells = 0;
  int total_size = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if ((with_hydro && cell_is_active_hydro(c, e)) ||
        (with_self_grav && cell_is_active_gravity(c, e)) ||
        (with_ext_grav && c->nodeID == nodeID &&
         cell_is_active_gravity(c, e)) ||
        (with_feedback && cell_is_active_stars(c, e)) ||
        (with_stars && c->nodeID == nodeID && cell_is_active_stars(c, e)) ||
        (with_sinks && cell_is_active_sinks(c, e)) ||
        (with_black_holes && cell_is_active_black_holes(c, e))) {

      if (num_active_top_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_top_cells], sizeof(int));
      num_active_top_cells += 1;
      const int depth = min(c->maxdepth, unskip_level);
      total_size += 1 << (3 * depth);
    }
  }

  /* What kind of tasks do we have? */
  struct unskip_data data;
  bzero(&data, sizeof(struct unskip_data));
  int multiplier = 0;
  if (with_hydro) {
    data.task_types[multiplier] = task_broad_types_hydro;
    multiplier++;
  }
  if (with_self_grav || with_ext_grav) {
    data.task_types[multiplier] = task_broad_types_gravity;
    multiplier++;
  }
  if (with_feedback || with_stars) {
    data.task_types[multiplier] = task_broad_types_stars;
    multiplier++;
  }
  if (with_sinks) {
    data.task_types[multiplier] = task_broad_types_sinks;
    multiplier++;
  }
  if (with_black_holes) {
    data.task_types[multiplier] = task_broad_types_black_holes;
    multiplier++;
  }
  if (with_rt) {
    data.task_types[multiplier] = task_broad_types_rt;
    multiplier++;
  }

  /* Create the list of active cells for the unskip. */
  struct cell **local_active_cells;
  local_active_cells = (struct cell **) malloc(
      multiplier * total_size * sizeof(struct cell*));
  if (local_active_cells == NULL) {
    error("Failed to allocate the cell array.");
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Reset to 0 in order to ensure a segfault. */
  bzero(local_active_cells, total_size);
#endif

  /* Construct the cell array */
  struct cell **current = local_active_cells;
  for(int i = 0; i < num_active_top_cells; i++) {
    struct cell *c = &s->cells_top[local_cells[i]];
    current = engine_unskip_copy_cell(c, current, unskip_level, e);
  }
  const int number_cells = current - local_active_cells;

  /* Make blind copies of the list */
  for (int m = 1; m < multiplier; m++) {
    memcpy(local_active_cells + m * number_cells, local_active_cells,
           number_cells * sizeof(struct cell *));
  }


  /* We now have a list of local active cells duplicated as many times as
   * we have broad task types. We can now release all the threads on the list */

  data.e = e;
  data.list_base = local_active_cells;
  data.num_active_cells = number_cells;
  data.multiplier = multiplier;

  /* Activate all the regular tasks */
  threadpool_map(&e->threadpool, engine_do_unskip_mapper, local_active_cells,
                 number_cells * multiplier, sizeof(struct cell*), /*chunk=*/1,
                 &data);

#ifdef WITH_PROFILER
  ProfilerStop();
#endif  // WITH_PROFILER

  /* Free stuff? */
  free(local_active_cells);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void engine_unskip_timestep_communications_mapper(void *map_data,
                                                  int num_elements,
                                                  void *extra_data) {
  /* Unpack the data */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *const tasks = (struct task *)map_data;
  const int nr_tasks = num_elements;

  /* Unskip the tasks in this part of the array */
  for (int i = 0; i < nr_tasks; ++i) {

    struct task *const t = &tasks[i];

    if (t->type == task_type_send || t->type == task_type_recv) {

      if (t->subtype == task_subtype_tend_part ||
          t->subtype == task_subtype_tend_gpart) {

        scheduler_activate(s, t);
      }
    }
  }
}

/**
 * @brief Blindly unskips all the tend communications for #part and #gpart.
 *
 * This function is only necessary when running with the time-step limiter
 * or the time-step synchronization policy as the time-steps of inactive
 * sections of the tree might have been changed by these tasks.
 *
 * @param e The #engine.
 */
void engine_unskip_timestep_communications(struct engine *e) {

#ifdef WITH_MPI

  const ticks tic = getticks();

  struct scheduler *s = &e->sched;
  struct task *tasks = e->sched.tasks;
  const int nr_tasks = e->sched.nr_tasks;

  /* Activate all the part and gpart ti_end tasks */
  threadpool_map(&e->threadpool, engine_unskip_timestep_communications_mapper,
                 tasks, nr_tasks, sizeof(struct task),
                 threadpool_auto_chunk_size, s);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
