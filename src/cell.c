/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Switch off timers. */
#ifdef TIMER
#undef TIMER
#endif

/* This object's header. */
#include "cell.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "drift.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "memswap.h"
#include "minmax.h"
#include "scheduler.h"
#include "space.h"
#include "timers.h"

/* Global variables. */
int cell_next_tag = 0;

/**
 * @brief Get the size of the cell subtree.
 *
 * @param c The #cell.
 */
int cell_getsize(struct cell *c) {

  /* Number of cells in this subtree. */
  int count = 1;

  /* Sum up the progeny if split. */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) count += cell_getsize(c->progeny[k]);

  /* Return the final count. */
  return count;
}

/**
 * @brief Unpack the data of a given cell and its sub-cells.
 *
 * @param pc An array of packed #pcell.
 * @param c The #cell in which to unpack the #pcell.
 * @param s The #space in which the cells are created.
 *
 * @return The number of cells created.
 */
int cell_unpack(struct pcell *pc, struct cell *c, struct space *s) {

#ifdef WITH_MPI

  /* Unpack the current pcell. */
  c->h_max = pc->h_max;
  c->ti_end_min = pc->ti_end_min;
  c->ti_end_max = pc->ti_end_max;
  c->ti_old = pc->ti_old;
  c->count = pc->count;
  c->gcount = pc->gcount;
  c->tag = pc->tag;

  /* Number of new cells created. */
  int count = 1;

  /* Fill the progeny recursively, depth-first. */
  for (int k = 0; k < 8; k++)
    if (pc->progeny[k] >= 0) {
      struct cell *temp;
      space_getcells(s, 1, &temp);
      temp->count = 0;
      temp->gcount = 0;
      temp->loc[0] = c->loc[0];
      temp->loc[1] = c->loc[1];
      temp->loc[2] = c->loc[2];
      temp->width[0] = c->width[0] / 2;
      temp->width[1] = c->width[1] / 2;
      temp->width[2] = c->width[2] / 2;
      temp->dmin = c->dmin / 2;
      if (k & 4) temp->loc[0] += temp->width[0];
      if (k & 2) temp->loc[1] += temp->width[1];
      if (k & 1) temp->loc[2] += temp->width[2];
      temp->depth = c->depth + 1;
      temp->split = 0;
      temp->dx_max = 0.f;
      temp->nodeID = c->nodeID;
      temp->parent = c;
      temp->ti_old = c->ti_old;
      c->progeny[k] = temp;
      c->split = 1;
      count += cell_unpack(&pc[pc->progeny[k]], temp, s);
    }

  /* Return the total number of unpacked cells. */
  c->pcell_size = count;
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Link the cells recursively to the given #part array.
 *
 * @param c The #cell.
 * @param parts The #part array.
 *
 * @return The number of particles linked.
 */
int cell_link_parts(struct cell *c, struct part *parts) {

  c->parts = parts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_parts(c->progeny[k], &parts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->count;
}

/**
 * @brief Link the cells recursively to the given #gpart array.
 *
 * @param c The #cell.
 * @param gparts The #gpart array.
 *
 * @return The number of particles linked.
 */
int cell_link_gparts(struct cell *c, struct gpart *gparts) {

  c->gparts = gparts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_gparts(c->progeny[k], &gparts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->gcount;
}

/**
 * @brief Pack the data of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param pc Pointer to an array of packed cells in which the
 *      cells will be packed.
 *
 * @return The number of packed cells.
 */
int cell_pack(struct cell *c, struct pcell *pc) {

#ifdef WITH_MPI

  /* Start by packing the data of the current cell. */
  pc->h_max = c->h_max;
  pc->ti_end_min = c->ti_end_min;
  pc->ti_end_max = c->ti_end_max;
  pc->ti_old = c->ti_old;
  pc->count = c->count;
  pc->gcount = c->gcount;
  c->tag = pc->tag = atomic_inc(&cell_next_tag) % cell_max_tag;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      pc->progeny[k] = count;
      count += cell_pack(c->progeny[k], &pc[count]);
    } else
      pc->progeny[k] = -1;

  /* Return the number of packed cells used. */
  c->pcell_size = count;
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the time information of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param ti_ends (output) The time information we pack into
 *
 * @return The number of packed cells.
 */
int cell_pack_ti_ends(struct cell *c, int *ti_ends) {

#ifdef WITH_MPI

  /* Pack this cell's data. */
  ti_ends[0] = c->ti_end_min;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_ti_ends(c->progeny[k], &ti_ends[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the time information of a given cell and its sub-cells.
 *
 * @param c The #cell
 * @param ti_ends The time information to unpack
 *
 * @return The number of cells created.
 */
int cell_unpack_ti_ends(struct cell *c, int *ti_ends) {

#ifdef WITH_MPI

  /* Unpack this cell's data. */
  c->ti_end_min = ti_ends[0];

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_ti_ends(c->progeny[k], &ti_ends[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Lock a cell for access to its array of #part and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_locktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to lock this cell. */
  if (c->hold || lock_trylock(&c->lock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->hold) {

    /* Unlock this cell. */
    if (lock_unlock(&c->lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {

    /* Lock this cell. */
    if (lock_trylock(&finger->lock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->hold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->lock) != 0) error("Failed to unlock cell.");
  }

  /* If we reached the top of the tree, we're done. */
  if (finger == NULL) {
    TIMER_TOC(timer_locktree);
    return 0;
  }

  /* Otherwise, we hit a snag. */
  else {

    /* Undo the holds up to finger. */
    for (struct cell *finger2 = c->parent; finger2 != finger;
         finger2 = finger2->parent)
      atomic_dec(&finger2->hold);

    /* Unlock this cell. */
    if (lock_unlock(&c->lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Lock a cell for access to its array of #gpart and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_glocktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to lock this cell. */
  if (c->ghold || lock_trylock(&c->glock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->ghold) {

    /* Unlock this cell. */
    if (lock_unlock(&c->glock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {

    /* Lock this cell. */
    if (lock_trylock(&finger->glock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->ghold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->glock) != 0) error("Failed to unlock cell.");
  }

  /* If we reached the top of the tree, we're done. */
  if (finger == NULL) {
    TIMER_TOC(timer_locktree);
    return 0;
  }

  /* Otherwise, we hit a snag. */
  else {

    /* Undo the holds up to finger. */
    for (struct cell *finger2 = c->parent; finger2 != finger;
         finger2 = finger2->parent)
      atomic_dec(&finger2->ghold);

    /* Unlock this cell. */
    if (lock_unlock(&c->glock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Unlock a cell's parents for access to #part array.
 *
 * @param c The #cell.
 */
void cell_unlocktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->hold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Unlock a cell's parents for access to #gpart array.
 *
 * @param c The #cell.
 */
void cell_gunlocktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->glock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->ghold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Sort the parts into eight bins along the given pivots.
 *
 * @param c The #cell array to be sorted.
 * @param buff A buffer with at least max(c->count, c->gcount) entries,
 *        used for sorting indices.
 * @param gbuff A buffer with at least max(c->count, c->gcount) entries,
 *        used for sorting indices for the gparts.
 */
void cell_split(struct cell *c, struct cell_buff *buff,
                struct cell_buff *gbuff) {

  const int count = c->count, gcount = c->gcount;
  const double pivot[3] = {c->loc[0] + c->width[0] / 2,
                           c->loc[1] + c->width[1] / 2,
                           c->loc[2] + c->width[2] / 2};
  int bucket_count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int bucket_offset[9];

  /* Fill the buffer with the indices. */
  for (int k = 0; k < count; k++) {
    const int bid = (buff[k].x[0] > pivot[0]) * 4 +
                    (buff[k].x[1] > pivot[1]) * 2 + (buff[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    buff[k].ind = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = buff[k].ind;
      if (bid != bucket) {
        struct cell_buff temp_buff = buff[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (buff[j].ind == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&buff[j], &temp_buff, sizeof(struct cell_buff));
          bid = temp_buff.ind;
        }
        buff[k] = temp_buff;
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->count = bucket_count[k];
    c->progeny[k]->parts = &c->parts[bucket_offset[k]];
    c->progeny[k]->xparts = &c->xparts[bucket_offset[k]];
  }

  /* Now do the same song and dance for the gparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < gcount; k++) {
    const int bid = (gbuff[k].x[0] > pivot[0]) * 4 +
                    (gbuff[k].x[1] > pivot[1]) * 2 + (gbuff[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    gbuff[k].ind = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = gbuff[k].ind;
      if (bid != bucket) {
        struct cell_buff temp_buff = gbuff[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (gbuff[j].ind == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&gbuff[j], &temp_buff, sizeof(struct cell_buff));
          bid = temp_buff.ind;
        }
        gbuff[k] = temp_buff;
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->gcount = bucket_count[k];
    c->progeny[k]->gparts = &c->gparts[bucket_offset[k]];
  }
}

/**
 * @brief Re-arrange the #part array accourding to the indices stored in the
 * buff array.
 *
 * @param parts The array of #part to re-arrange.
 * @param xparts The corresponding array of #xpart.
 * @param count The size of the parts array.
 * @param buff An array of #cell_buff containing the indices of parts.
 */
void cell_reorder_parts(struct part *parts, struct xpart *xparts, int count,
                        struct cell_buff *buff) {

  for (int i = 0; i < count; ++i) {

    // If this particle is in the right place already, do nothing.
    if (buff[i].offset == i) continue;
    
    // Swap this particle to its correct location.
    struct part temp;
    struct xpart xtemp;
    vector_memcpy(&temp, &parts[i], sizeof(struct part));
    vector_memcpy(&xtemp, &xparts[i], sizeof(struct xpart));
    int j = i;
    while (1) {
      const int k = buff[j].offset;
      buff[j].offset = j;
      if (k == i) break;
      vector_memcpy(&parts[j], &parts[k], sizeof(struct part));
      vector_memcpy(&xparts[j], &xparts[k], sizeof(struct xpart));
      j = k;
    }
    vector_memcpy(&parts[j], &temp, sizeof(struct part));
    vector_memcpy(&xparts[j], &xtemp, sizeof(struct xpart));
  }
}

/**
 * @brief Re-arrange the #gpart array accourding to the indices stored in the
 * buff array.
 *
 * @param gparts The array of #gpart to re-arrange.
 * @param gcount The size of the gparts array.
 * @param buff An array of #cell_buff containing the indices of gparts.
 */
void cell_reorder_gparts(struct gpart *gparts, int gcount,
                         struct cell_buff *buff) {

  for (int i = 0; i < gcount; ++i) {

    // If this particle is in the right place already, do nothing.
    if (buff[i].offset == i) continue;
    
    // Swap this particle to its correct location.
    struct gpart temp;
    vector_memcpy(&temp, &gparts[i], sizeof(struct gpart));
    int j = i;
    while (1) {
      const int k = buff[j].offset;
      buff[j].offset = j;
      if (k == i) break;
      vector_memcpy(&gparts[j], &gparts[k], sizeof(struct gpart));
      j = k;
    }
    vector_memcpy(&gparts[j], &temp, sizeof(struct gpart));
  }
}

/**
 * @brief Recursively update the h_max and time-step information of a cell and
 * all its progenitors.
 *
 * @param c The #cell to act on.
 */
void cell_init_counters(struct cell *c) {

  const int count = c->count;
  const int gcount = c->gcount;
  float h_max = 0.0f;
  int ti_end_min = max_nr_timesteps, ti_end_max = 0;
  struct part *parts = c->parts;
  struct xpart *xparts = c->xparts;
  struct gpart *gparts = c->gparts;

  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        /* Recurse */
        cell_init_counters(c->progeny[k]);

        /* Update */
        h_max = max(h_max, c->progeny[k]->h_max);
        ti_end_min = min(ti_end_min, c->progeny[k]->ti_end_min);
        ti_end_max = max(ti_end_max, c->progeny[k]->ti_end_max);
      }
    }

  } else {

    /* Get dt_min/dt_max. */
    for (int k = 0; k < count; k++) {
      const struct part *p = &parts[k];
      struct xpart *xp = &xparts[k];
      const float h = p->h;
      const int ti_end = p->ti_end;
      xp->x_diff[0] = 0.f;
      xp->x_diff[1] = 0.f;
      xp->x_diff[2] = 0.f;
      if (h > h_max) h_max = h;
      if (ti_end < ti_end_min) ti_end_min = ti_end;
      if (ti_end > ti_end_max) ti_end_max = ti_end;
    }
    for (int k = 0; k < gcount; k++) {
      struct gpart *gp = &gparts[k];
      const int ti_end = gp->ti_end;
      gp->x_diff[0] = 0.f;
      gp->x_diff[1] = 0.f;
      gp->x_diff[2] = 0.f;
      if (ti_end < ti_end_min) ti_end_min = ti_end;
      if (ti_end > ti_end_max) ti_end_max = ti_end;
    }
  }

  /* Set the values for this cell. */
  c->h_max = h_max;
  c->ti_end_min = ti_end_min;
  c->ti_end_max = ti_end_max;
}

/**
 * @brief Sanitizes the smoothing length values of cells by setting large
 * outliers to more sensible values.
 *
 * We compute the mean and standard deviation of the smoothing lengths in
 * logarithmic space and limit values to mean + 4 sigma.
 *
 * @param c The cell.
 */
void cell_sanitize(struct cell *c) {

  const int count = c->count;
  struct part *parts = c->parts;

  /* First collect some statistics */
  float h_mean = 0.f, h_mean2 = 0.f;
  float h_min = FLT_MAX, h_max = 0.f;
  for (int i = 0; i < count; ++i) {

    const float h = logf(parts[i].h);
    h_mean += h;
    h_mean2 += h * h;
    h_max = max(h_max, h);
    h_min = min(h_min, h);
  }
  h_mean /= count;
  h_mean2 /= count;
  const float h_var = h_mean2 - h_mean * h_mean;
  const float h_std = (h_var > 0.f) ? sqrtf(h_var) : 0.1f * h_mean;

  /* Choose a cut */
  const float h_limit = expf(h_mean + 4.f * h_std);

  /* Be verbose this is not innocuous */
  message("Cell properties: h_min= %f h_max= %f geometric mean= %f.",
          expf(h_min), expf(h_max), expf(h_mean));

  if (c->h_max > h_limit) {

    message("Smoothing lengths will be limited to (mean + 4sigma)= %f.",
            h_limit);

    /* Apply the cut */
    for (int i = 0; i < count; ++i) parts->h = min(parts[i].h, h_limit);

    c->h_max = h_limit;

  } else {

    message("Smoothing lengths will not be limited.");
  }
}

/**
 * @brief Converts hydro quantities to a valid state after the initial density
 * calculation
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_convert_hydro(struct cell *c, void *data) {

  struct part *p = c->parts;

  for (int i = 0; i < c->count; ++i) {
    hydro_convert_quantities(&p[i]);
  }
}

/**
 * @brief Cleans the links in a given cell.
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_clean_links(struct cell *c, void *data) {
  c->density = NULL;
  c->gradient = NULL;
  c->force = NULL;
  c->grav = NULL;
}

/**
 * @brief Checks that a cell is at the current point in time
 *
 * Calls error() if the cell is not at the current time.
 *
 * @param c Cell to act upon
 * @param data The current time on the integer time-line
 */
void cell_check_drift_point(struct cell *c, void *data) {

  const int ti_current = *(int *)data;

  if (c->ti_old != ti_current && c->nodeID == engine_rank)
    error("Cell in an incorrect time-zone! c->ti_old=%d ti_current=%d",
          c->ti_old, ti_current);
}

/**
 * @brief Checks whether the cells are direct neighbours ot not. Both cells have
 * to be of the same size
 *
 * @param ci First #cell.
 * @param cj Second #cell.
 *
 * @todo Deal with periodicity.
 */
int cell_are_neighbours(const struct cell *restrict ci,
                        const struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->width[0] != cj->width[0]) error("Cells of different size !");
#endif

  /* Maximum allowed distance */
  const double min_dist =
      1.2 * ci->width[0]; /* 1.2 accounts for rounding errors */

  /* (Manhattan) Distance between the cells */
  for (int k = 0; k < 3; k++) {
    const double center_i = ci->loc[k];
    const double center_j = cj->loc[k];
    if (fabs(center_i - center_j) > min_dist) return 0;
  }

  return 1;
}

/**
 * @brief Computes the multi-pole brutally and compare to the
 * recursively computed one.
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_check_multipole(struct cell *c, void *data) {

  struct multipole ma;

  if (c->gcount > 0) {

    /* Brute-force calculation */
    multipole_init(&ma, c->gparts, c->gcount);

    /* Compare with recursive one */
    struct multipole mb = c->multipole;

    if (fabsf(ma.mass - mb.mass) / fabsf(ma.mass + mb.mass) > 1e-5)
      error("Multipole masses are different (%12.15e vs. %12.15e)", ma.mass,
            mb.mass);

    for (int k = 0; k < 3; ++k)
      if (fabs(ma.CoM[k] - mb.CoM[k]) / fabs(ma.CoM[k] + mb.CoM[k]) > 1e-5)
        error("Multipole CoM are different (%12.15e vs. %12.15e", ma.CoM[k],
              mb.CoM[k]);

#if const_gravity_multipole_order >= 2
    if (fabsf(ma.I_xx - mb.I_xx) / fabsf(ma.I_xx + mb.I_xx) > 1e-5 &&
        ma.I_xx > 1e-9)
      error("Multipole I_xx are different (%12.15e vs. %12.15e)", ma.I_xx,
            mb.I_xx);
    if (fabsf(ma.I_yy - mb.I_yy) / fabsf(ma.I_yy + mb.I_yy) > 1e-5 &&
        ma.I_yy > 1e-9)
      error("Multipole I_yy are different (%12.15e vs. %12.15e)", ma.I_yy,
            mb.I_yy);
    if (fabsf(ma.I_zz - mb.I_zz) / fabsf(ma.I_zz + mb.I_zz) > 1e-5 &&
        ma.I_zz > 1e-9)
      error("Multipole I_zz are different (%12.15e vs. %12.15e)", ma.I_zz,
            mb.I_zz);
    if (fabsf(ma.I_xy - mb.I_xy) / fabsf(ma.I_xy + mb.I_xy) > 1e-5 &&
        ma.I_xy > 1e-9)
      error("Multipole I_xy are different (%12.15e vs. %12.15e)", ma.I_xy,
            mb.I_xy);
    if (fabsf(ma.I_xz - mb.I_xz) / fabsf(ma.I_xz + mb.I_xz) > 1e-5 &&
        ma.I_xz > 1e-9)
      error("Multipole I_xz are different (%12.15e vs. %12.15e)", ma.I_xz,
            mb.I_xz);
    if (fabsf(ma.I_yz - mb.I_yz) / fabsf(ma.I_yz + mb.I_yz) > 1e-5 &&
        ma.I_yz > 1e-9)
      error("Multipole I_yz are different (%12.15e vs. %12.15e)", ma.I_yz,
            mb.I_yz);
#endif
  }
}

/**
 * @brief Frees up the memory allocated for this #cell.
 *
 * @param c The #cell.
 */
void cell_clean(struct cell *c) {

  free(c->sort);

  /* Recurse */
  for (int k = 0; k < 8; k++)
    if (c->progeny[k]) cell_clean(c->progeny[k]);
}

/**
 * @brief Checks whether a given cell needs drifting or not.
 *
 * @param c the #cell.
 * @param e The #engine (holding current time information).
 *
 * @return 1 If the cell needs drifting, 0 otherwise.
 */
int cell_is_drift_needed(struct cell *c, const struct engine *e) {

  /* Do we have at least one active particle in the cell ?*/
  if (cell_is_active(c, e)) return 1;

  /* Loop over the pair tasks that involve this cell */
  for (struct link *l = c->density; l != NULL; l = l->next) {

    if (l->t->type != task_type_pair && l->t->type != task_type_sub_pair)
      continue;

    /* Is the other cell in the pair active ? */
    if ((l->t->ci == c && cell_is_active(l->t->cj, e)) ||
        (l->t->cj == c && cell_is_active(l->t->ci, e)))
      return 1;
  }

  /* No neighbouring cell has active particles. Drift not necessary */
  return 0;
}

/**
 * @brief Un-skips all the tasks associated with a given cell and checks
 * if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_tasks(struct cell *c, struct scheduler *s) {

  /* Un-skip the density tasks involved with this cell. */
  for (struct link *l = c->density; l != NULL; l = l->next) {
    struct task *t = l->t;
    const struct cell *ci = t->ci;
    const struct cell *cj = t->cj;
    scheduler_activate(s, t);

    /* Set the correct sorting flags */
    if (t->type == task_type_pair) {
      if (!(ci->sorted & (1 << t->flags))) {
        atomic_or(&ci->sorts->flags, (1 << t->flags));
        scheduler_activate(s, ci->sorts);
      }
      if (!(cj->sorted & (1 << t->flags))) {
        atomic_or(&cj->sorts->flags, (1 << t->flags));
        scheduler_activate(s, cj->sorts);
      }
    }

    /* Check whether there was too much particle motion */
    if (t->type == task_type_pair || t->type == task_type_sub_pair) {
      if (t->tight &&
          (max(ci->h_max, cj->h_max) + ci->dx_max + cj->dx_max > cj->dmin ||
           ci->dx_max > space_maxreldx * ci->h_max ||
           cj->dx_max > space_maxreldx * cj->h_max))
        return 1;

#ifdef WITH_MPI
      /* Activate the send/recv flags. */
      if (ci->nodeID != engine_rank) {

        /* Activate the tasks to recv foreign cell ci's data. */
        scheduler_activate(s, ci->recv_xv);
        scheduler_activate(s, ci->recv_rho);
        scheduler_activate(s, ci->recv_ti);

        /* Look for the local cell cj's send tasks. */
        struct link *l = NULL;
        for (l = cj->send_xv; l != NULL && l->t->cj->nodeID != ci->nodeID;
             l = l->next)
          ;
        if (l == NULL) error("Missing link to send_xv task.");
        scheduler_activate(s, l->t);

        if (cj->super->drift)
          scheduler_activate(s, cj->super->drift);
        else
          error("Drift task missing !");

        for (l = cj->send_rho; l != NULL && l->t->cj->nodeID != ci->nodeID;
             l = l->next)
          ;
        if (l == NULL) error("Missing link to send_rho task.");
        scheduler_activate(s, l->t);

        for (l = cj->send_ti; l != NULL && l->t->cj->nodeID != ci->nodeID;
             l = l->next)
          ;
        if (l == NULL) error("Missing link to send_ti task.");
        scheduler_activate(s, l->t);

      } else if (cj->nodeID != engine_rank) {

        /* Activate the tasks to recv foreign cell cj's data. */
        scheduler_activate(s, cj->recv_xv);
        scheduler_activate(s, cj->recv_rho);
        scheduler_activate(s, cj->recv_ti);

        /* Look for the local cell ci's send tasks. */
        struct link *l = NULL;
        for (l = ci->send_xv; l != NULL && l->t->cj->nodeID != cj->nodeID;
             l = l->next)
          ;
        if (l == NULL) error("Missing link to send_xv task.");
        scheduler_activate(s, l->t);

        if (ci->super->drift)
          scheduler_activate(s, ci->super->drift);
        else
          error("Drift task missing !");

        for (l = ci->send_rho; l != NULL && l->t->cj->nodeID != cj->nodeID;
             l = l->next)
          ;
        if (l == NULL) error("Missing link to send_rho task.");
        scheduler_activate(s, l->t);

        for (l = ci->send_ti; l != NULL && l->t->cj->nodeID != cj->nodeID;
             l = l->next)
          ;
        if (l == NULL) error("Missing link to send_ti task.");
        scheduler_activate(s, l->t);
      }
#endif
    }
  }

  /* Unskip all the other task types. */
  for (struct link *l = c->gradient; l != NULL; l = l->next)
    scheduler_activate(s, l->t);
  for (struct link *l = c->force; l != NULL; l = l->next)
    scheduler_activate(s, l->t);
  for (struct link *l = c->grav; l != NULL; l = l->next)
    scheduler_activate(s, l->t);
  if (c->extra_ghost != NULL) scheduler_activate(s, c->extra_ghost);
  if (c->ghost != NULL) scheduler_activate(s, c->ghost);
  if (c->init != NULL) scheduler_activate(s, c->init);
  if (c->drift != NULL) scheduler_activate(s, c->drift);
  if (c->kick != NULL) scheduler_activate(s, c->kick);
  if (c->cooling != NULL) scheduler_activate(s, c->cooling);
  if (c->sourceterms != NULL) scheduler_activate(s, c->sourceterms);

  return 0;
}

/**
 * @brief Set the super-cell pointers for all cells in a hierarchy.
 *
 * @param c The top-level #cell to play with.
 * @param super Pointer to the deepest cell with tasks in this part of the tree.
 */
void cell_set_super(struct cell *c, struct cell *super) {

  /* Are we in a cell with some kind of self/pair task ? */
  if (super == NULL && c->nr_tasks > 0) super = c;

  /* Set the super-cell */
  c->super = super;

  /* Recurse */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_set_super(c->progeny[k], super);
}

/**
 * @brief Recursively drifts all particles and g-particles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 */
void cell_drift(struct cell *c, const struct engine *e) {

  const double timeBase = e->timeBase;
  const int ti_old = c->ti_old;
  const int ti_current = e->ti_current;
  struct part *const parts = c->parts;
  struct xpart *const xparts = c->xparts;
  struct gpart *const gparts = c->gparts;

  /* Drift from the last time the cell was drifted to the current time */
  const double dt = (ti_current - ti_old) * timeBase;
  float dx_max = 0.f, dx2_max = 0.f, h_max = 0.f;

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old) error("Attempt to drift to the past");

  /* Are we not in a leaf ? */
  if (c->split) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        cell_drift(cp, e);
        dx_max = max(dx_max, cp->dx_max);
        h_max = max(h_max, cp->h_max);
      }

  } else if (ti_current > ti_old) {

    /* Loop over all the g-particles in the cell */
    const size_t nr_gparts = c->gcount;
    for (size_t k = 0; k < nr_gparts; k++) {

      /* Get a handle on the gpart. */
      struct gpart *const gp = &gparts[k];

      /* Drift... */
      drift_gpart(gp, dt, timeBase, ti_old, ti_current);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = gp->x_diff[0] * gp->x_diff[0] +
                        gp->x_diff[1] * gp->x_diff[1] +
                        gp->x_diff[2] * gp->x_diff[2];
      dx2_max = (dx2_max > dx2) ? dx2_max : dx2;
    }

    /* Loop over all the particles in the cell */
    const size_t nr_parts = c->count;
    for (size_t k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* Drift... */
      drift_part(p, xp, dt, timeBase, ti_old, ti_current);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = xp->x_diff[0] * xp->x_diff[0] +
                        xp->x_diff[1] * xp->x_diff[1] +
                        xp->x_diff[2] * xp->x_diff[2];
      dx2_max = (dx2_max > dx2) ? dx2_max : dx2;

      /* Maximal smoothing length */
      h_max = (h_max > p->h) ? h_max : p->h;
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);

  } else {

    h_max = c->h_max;
    dx_max = c->dx_max;
  }

  /* Store the values */
  c->h_max = h_max;
  c->dx_max = dx_max;

  /* Update the time of the last drift */
  c->ti_old = ti_current;
}

/**
 * @brief Recursively check that the cell properties match the particle content.
 *
 * This is a debugging function. Don't call in actual code.
 *
 * @param c The #cell to check.
 */
void cell_check_content(struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS

  const int count = c->count;
  const int gcount = c->gcount;

  const double pos_min[3] = {c->loc[0], c->loc[1], c->loc[2]};
  const double pos_max[3] = {c->loc[0] + c->width[0], c->loc[1] + c->width[1],
                             c->loc[2] + c->width[2]};

  const double h_max = c->h_max;
  const int ti_end_min = c->ti_end_min;
  const int ti_end_max = c->ti_end_max;

  /* Check all the parts at this level */
  for (int i = 0; i < count; ++i) {

    /* Handle on the particle */
    const struct part *p = &c->parts[i];

    /* Check that particle is in the cell */
    for (int k = 0; k < 3; ++k) {
      if (p->x[k] < pos_min[k] || p->x[k] > pos_max[k])
        error(
            "Particle not in the cell ! c->loc=[%e %e %e] c->width=[%e %e %e] "
            "p->x=[%e %e %e]",
            pos_min[0], pos_min[1], pos_min[2], c->width[0], c->width[1],
            c->width[2], p->x[0], p->x[1], p->x[2]);
    }

    /* Check time-zone */
    if (p->ti_end < ti_end_min)
      error("Particle in wrong time-zone ! p->ti_end=%d c->ti_end_min=%d",
            p->ti_end, ti_end_min);

    if (p->ti_end > ti_end_max)
      error("Particle in wrong time-zone ! p->ti_end=%d c->ti_end_max=%d",
            p->ti_end, ti_end_max);

    /* Check smoothing length */
    if (p->h > h_max)
      error("Particle with wrong h ! p->h=%e c->h_max=%e", p->h, h_max);
  }

  /* Check all the parts at this level */
  for (int i = 0; i < gcount; ++i) {

    /* Handle on the particle */
    const struct gpart *gp = &c->gparts[i];

    /* Check that particle is in the cell */
    for (int k = 0; k < 3; ++k) {
      if (gp->x[k] < pos_min[k] || gp->x[k] > pos_max[k])
        error(
            "Particle not in the cell ! c->loc=[%e %e %e] c->width=[%e %e %e] "
            "gp->x=[%e %e %e]",
            pos_min[0], pos_min[1], pos_min[2], c->width[0], c->width[1],
            c->width[2], gp->x[0], gp->x[1], gp->x[2]);
    }

    /* Check time-zone */
    if (gp->ti_end < ti_end_min)
      error("Particle in wrong time-zone ! gp->ti_end=%d c->ti_end_min=%d",
            gp->ti_end, ti_end_min);

    if (gp->ti_end > ti_end_max)
      error("Particle in wrong time-zone ! gp->ti_end=%d c->ti_end_max=%d",
            gp->ti_end, ti_end_max);
  }

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) cell_check_content(c->progeny[k]);
    }
  }

#else

  error("Debugging function called in actual code.");

#endif
}
