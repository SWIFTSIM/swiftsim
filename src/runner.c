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
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "approx_math.h"
#include "atomic.h"
#include "const.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro_properties.h"
#include "hydro.h"
#include "minmax.h"
#include "scheduler.h"
#include "space.h"
#include "task.h"
#include "timers.h"

/* Orientation of the cell pairs */
const float runner_shift[13 * 3] = {
    5.773502691896258e-01,
    5.773502691896258e-01,
    5.773502691896258e-01,
    7.071067811865475e-01,
    7.071067811865475e-01,
    0.0,
    5.773502691896258e-01,
    5.773502691896258e-01,
    -5.773502691896258e-01,
    7.071067811865475e-01,
    0.0,
    7.071067811865475e-01,
    1.0,
    0.0,
    0.0,
    7.071067811865475e-01,
    0.0,
    -7.071067811865475e-01,
    5.773502691896258e-01,
    -5.773502691896258e-01,
    5.773502691896258e-01,
    7.071067811865475e-01,
    -7.071067811865475e-01,
    0.0,
    5.773502691896258e-01,
    -5.773502691896258e-01,
    -5.773502691896258e-01,
    0.0,
    7.071067811865475e-01,
    7.071067811865475e-01,
    0.0,
    1.0,
    0.0,
    0.0,
    7.071067811865475e-01,
    -7.071067811865475e-01,
    0.0,
    0.0,
    1.0,
};

/* Does the axis need flipping ? */
const char runner_flip[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* #define MY_CELL 428428428 */
/* #define DX 0.1 */
/* #define NX 1000 */
/* #define CELL_ID                                                  \ */
/*   ((int)(c->loc[0] / DX) * NX *NX + (int)(c->loc[1] / DX) * NX + \ */
/*    (int)(c->loc[2] / DX)) */
/* #define OUT \ */
/*   if (CELL_ID == MY_CELL) { \ */
/*     message(" cell= %d gcount=%d time=%f \n", CELL_ID, c->gcount,
 * r->e->time); \ */
/*   } */
/* #define OUT  message(" cell %d %d %f \n",CELL_ID, c->count, r->e->time); */
/* #define OUT  if(CELL_ID == MY_CELL) message("\n cell %f %f %f %d %d */
/* %f\n",c->loc[0],c->loc[1],c->loc[2], CELL_ID, c->count, r->e->time); */

/* Import the density loop functions. */
#define FUNCTION density
#include "runner_doiact.h"

/* Import the force loop functions. */
#undef FUNCTION
#define FUNCTION force
#include "runner_doiact.h"

/* Import the gravity loop functions. */
#include "runner_doiact_grav.h"

/**
 * @brief Calculate gravity acceleration from external potential
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_grav_external(struct runner *r, struct cell *c, int timer) {

  struct gpart *g, *gparts = c->gparts;
  int i, k, gcount = c->gcount;
  const int ti_current = r->e->ti_current;
  const struct external_potential *potential = r->e->external_potential;
  const struct phys_const *constants = r->e->physical_constants;

  TIMER_TIC;

  /* Recurse? */
  if (c->split) {
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_grav_external(r, c->progeny[k], 0);
    return;
  }

#ifdef TASK_VERBOSE
  OUT;
#endif

  /* Loop over the parts in this cell. */
  for (i = 0; i < gcount; i++) {

    /* Get a direct pointer on the part. */
    g = &gparts[i];

    /* Is this part within the time step? */
    if (g->ti_end <= ti_current) {

      external_gravity(potential, constants, g);

      /* /\* check for energy and angular momentum conservation - begin by */
      /*  * synchronizing velocity*\/ */
      /* const float dx = g->x[0] - r->e->potential->point_mass.x; */
      /* const float dy = g->x[1] - r->e->potential->point_mass.y; */
      /* const float dz = g->x[2] - r->e->potential->point_mass.z; */
      /* const float dr = sqrtf((dx * dx) + (dy * dy) + (dz * dz)); */
      /* const float rinv = 1.f / sqrtf(dx * dx + dy * dy + dz * dz); */

      /* const int current_dti = g->ti_end - g->ti_begin; */
      /* const float dt = 0.5f * current_dti * r->e->timeBase; */
      /* const float vx = g->v_full[0] + dt * g->a_grav[0]; */
      /* const float vy = g->v_full[1] + dt * g->a_grav[1]; */
      /* const float vz = g->v_full[2] + dt * g->a_grav[2]; */

      /* /\* E/L *\/ */
      /* float L[3], E; */
      /* E = 0.5 * ((vx * vx) + (vy * vy) + (vz * vz)) - */
      /*     r->e->physical_constants->newton_gravity * */
      /*         r->e->potential->point_mass.mass * rinv; */
      /* L[0] = dy * vz - dz * vy; */
      /* L[1] = dz * vx - dx * vz; */
      /* L[2] = dx * vy - dy * vx; */
      /* if (abs(g->id) == 1) { */
      /*   float v2 = vx * vx + vy * vy + vz * vz; */
      /*   float fg = r->e->physical_constants->newton_gravity * */
      /*              r->e->potential->point_mass.mass * rinv; */
      /*   float fga = sqrtf((g->a_grav[0] * g->a_grav[0]) + */
      /*                     (g->a_grav[1] * g->a_grav[1]) + */
      /*                     (g->a_grav[2] * g->a_grav[2])) * */
      /*               dr; */
      /*   // message("grav_external time= %f\t V_c^2= %f GM/r= %f E= %f L[2]=
       * %f */
      /*   // x= %f y= %f vx= %f vy= %f\n", r->e->time, v2, fg, E, L[2],
       * g->x[0], */
      /*   // g->x[1], vx, vy); */
      /*   message("%f\t %f %f %f %f %f %f %f %f %f %f %f %f %f\n", r->e->time,
       */
      /*           g->tx, g->tv, dt, v2, fg, fga, dr, E, L[2], g->x[0], g->x[1],
       */
      /*           vx, vy); */
      /*   /\* message(" G=%e M=%e\n", r->e->physical_constants->newton_gravity,
       */
      /*    * r->e->potential->point_mass.mass); *\/ */
      /*   /\* exit(-1); *\/ */
      /* } */
    }
  }
  if (timer) TIMER_TOC(timer_dograv_external);
}

/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */

void runner_do_sort_ascending(struct entry *sort, int N) {

  struct {
    short int lo, hi;
  } qstack[10];
  int qpos, i, j, lo, hi, imin;
  struct entry temp;
  float pivot;

  /* Sort parts in cell_i in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    if (hi - lo < 15) {
      for (i = lo; i < hi; i++) {
        imin = i;
        for (j = i + 1; j <= hi; j++)
          if (sort[j].d < sort[imin].d) imin = j;
        if (imin != i) {
          temp = sort[imin];
          sort[imin] = sort[i];
          sort[i] = temp;
        }
      }
    } else {
      pivot = sort[(lo + hi) / 2].d;
      i = lo;
      j = hi;
      while (i <= j) {
        while (sort[i].d < pivot) i++;
        while (sort[j].d > pivot) j--;
        if (i <= j) {
          if (i < j) {
            temp = sort[i];
            sort[i] = sort[j];
            sort[j] = temp;
          }
          i += 1;
          j -= 1;
        }
      }
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }
}

/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */

void runner_do_sort(struct runner *r, struct cell *c, int flags, int clock) {

  struct entry *finger;
  struct entry *fingers[8];
  struct part *parts = c->parts;
  struct entry *sort;
  int j, k, count = c->count;
  int i, ind, off[8], inds[8], temp_i, missing;
  // float shift[3];
  float buff[8], px[3];

  TIMER_TIC

  /* Clean-up the flags, i.e. filter out what's already been sorted. */
  flags &= ~c->sorted;
  if (flags == 0) return;

  /* start by allocating the entry arrays. */
  if (c->sort == NULL || c->sortsize < count) {
    if (c->sort != NULL) free(c->sort);
    c->sortsize = count * 1.1;
    if ((c->sort = (struct entry *)malloc(sizeof(struct entry) *
                                          (c->sortsize + 1) * 13)) == NULL)
      error("Failed to allocate sort memory.");
  }
  sort = c->sort;

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    for (k = 0; k < 8; k++) {
      if (c->progeny[k] == NULL) continue;
      missing = flags & ~c->progeny[k]->sorted;
      if (missing) runner_do_sort(r, c->progeny[k], missing, 0);
    }

    /* Loop over the 13 different sort arrays. */
    for (j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      for (off[0] = 0, k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      for (k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->count > 0) {
          fingers[k] = &c->progeny[k]->sort[j * (c->progeny[k]->count + 1)];
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (i = 0; i < 7; i++)
        for (k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      finger = &sort[j * (count + 1)];
      for (ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      sort[j * (count + 1) + count].d = FLT_MAX;
      sort[j * (count + 1) + count].i = 0;

      /* Mark as sorted. */
      c->sorted |= (1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Fill the sort array. */
    for (k = 0; k < count; k++) {
      px[0] = parts[k].x[0];
      px[1] = parts[k].x[1];
      px[2] = parts[k].x[2];
      for (j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          sort[j * (count + 1) + k].i = k;
          sort[j * (count + 1) + k].d = px[0] * runner_shift[3 * j + 0] +
                                        px[1] * runner_shift[3 * j + 1] +
                                        px[2] * runner_shift[3 * j + 2];
        }
    }

    /* Add the sentinel and sort. */
    for (j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        sort[j * (count + 1) + count].d = FLT_MAX;
        sort[j * (count + 1) + count].i = 0;
        runner_do_sort_ascending(&sort[j * (count + 1)], count);
        c->sorted |= (1 << j);
      }
  }

  /* Verify the sorting. */
  /* for ( j = 0 ; j < 13 ; j++ ) {
      if ( !( flags & (1 << j) ) )
          continue;
      finger = &sort[ j*(count + 1) ];
      for ( k = 1 ; k < count ; k++ ) {
          if ( finger[k].d < finger[k-1].d )
              error( "Sorting failed, ascending array." );
          if ( finger[k].i >= count )
              error( "Sorting failed, indices borked." );
          }
      } */

  if (clock) TIMER_TOC(timer_dosort);
}

void runner_do_gsort(struct runner *r, struct cell *c, int flags, int clock) {

  struct entry *finger;
  struct entry *fingers[8];
  struct gpart *gparts = c->gparts;
  struct entry *gsort;
  int j, k, count = c->gcount;
  int i, ind, off[8], inds[8], temp_i, missing;
  // float shift[3];
  float buff[8], px[3];

  TIMER_TIC

  /* Clean-up the flags, i.e. filter out what's already been sorted. */
  flags &= ~c->gsorted;
  if (flags == 0) return;

  /* start by allocating the entry arrays. */
  if (c->gsort == NULL || c->gsortsize < count) {
    if (c->gsort != NULL) free(c->gsort);
    c->gsortsize = count * 1.1;
    if ((c->gsort = (struct entry *)malloc(sizeof(struct entry) *
                                           (c->gsortsize + 1) * 13)) == NULL)
      error("Failed to allocate sort memory.");
  }
  gsort = c->gsort;

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    for (k = 0; k < 8; k++) {
      if (c->progeny[k] == NULL) continue;
      missing = flags & ~c->progeny[k]->gsorted;
      if (missing) runner_do_gsort(r, c->progeny[k], missing, 0);
    }

    /* Loop over the 13 different sort arrays. */
    for (j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      for (off[0] = 0, k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->gcount;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      for (k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->gcount > 0) {
          fingers[k] = &c->progeny[k]->gsort[j * (c->progeny[k]->gcount + 1)];
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (i = 0; i < 7; i++)
        for (k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      finger = &gsort[j * (count + 1)];
      for (ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      gsort[j * (count + 1) + count].d = FLT_MAX;
      gsort[j * (count + 1) + count].i = 0;

      /* Mark as sorted. */
      c->gsorted |= (1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Fill the sort array. */
    for (k = 0; k < count; k++) {
      px[0] = gparts[k].x[0];
      px[1] = gparts[k].x[1];
      px[2] = gparts[k].x[2];
      for (j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          gsort[j * (count + 1) + k].i = k;
          gsort[j * (count + 1) + k].d = px[0] * runner_shift[3 * j + 0] +
                                         px[1] * runner_shift[3 * j + 1] +
                                         px[2] * runner_shift[3 * j + 2];
        }
    }

    /* Add the sentinel and sort. */
    for (j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        gsort[j * (count + 1) + count].d = FLT_MAX;
        gsort[j * (count + 1) + count].i = 0;
        runner_do_sort_ascending(&gsort[j * (count + 1)], count);
        c->gsorted |= (1 << j);
      }
  }

  /* Verify the sorting. */
  /* for ( j = 0 ; j < 13 ; j++ ) {
      if ( !( flags & (1 << j) ) )
          continue;
      finger = &c->gsort[ j*(count + 1) ];
      for ( k = 1 ; k < count ; k++ ) {
          if ( finger[k].d < finger[k-1].d )
              error( "Sorting failed, ascending array." );
          if ( finger[k].i < 0 || finger[k].i >= count )
              error( "Sorting failed, indices borked." );
          }
      } */

  if (clock) TIMER_TOC(timer_dosort);
}

/**
 * @brief Initialize the particles before the density calculation
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_init(struct runner *r, struct cell *c, int timer) {

  struct part *const parts = c->parts;
  struct gpart *const gparts = c->gparts;
  const int count = c->count;
  const int gcount = c->gcount;
  const int ti_current = r->e->ti_current;

  TIMER_TIC;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_init(r, c->progeny[k], 0);
    return;
  } else {

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      struct part *const p = &parts[i];

      if (p->ti_end <= ti_current) {

        /* Get ready for a density calculation */
        hydro_init_part(p);
      }
    }

    /* Loop over the gparts in this cell. */
    for (int i = 0; i < gcount; i++) {

      /* Get a direct pointer on the part. */
      struct gpart *const gp = &gparts[i];

      if (gp->ti_end <= ti_current) {

        /* Get ready for a density calculation */
        gravity_init_part(gp);
      }
    }
  }

  if (timer) TIMER_TOC(timer_init);
}

/**
 * @brief Intermediate task between density and force
 *
 * @param r The runner thread.
 * @param c The cell.
 */

void runner_do_ghost(struct runner *r, struct cell *c) {

  struct part *p, *parts = c->parts;
  struct xpart *xp, *xparts = c->xparts;
  struct cell *finger;
  int redo, count = c->count;
  int *pid;
  float h_corr;
  const int ti_current = r->e->ti_current;
  const double timeBase = r->e->timeBase;
  const float target_wcount = r->e->hydro_properties->target_neighbours;
  const float max_wcount =
      target_wcount + r->e->hydro_properties->delta_neighbours;
  const float min_wcount =
      target_wcount - r->e->hydro_properties->delta_neighbours;
  const int max_smoothing_iter =
      r->e->hydro_properties->max_smoothing_iterations;

  TIMER_TIC;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_ghost(r, c->progeny[k]);
    return;
  }

  /* Init the IDs that have to be updated. */
  if ((pid = (int *)alloca(sizeof(int) * count)) == NULL)
    error("Call to alloca failed.");
  for (int k = 0; k < count; k++) pid[k] = k;

  /* While there are particles that need to be updated... */
  for (int num_reruns = 0; count > 0 && num_reruns < max_smoothing_iter;
       num_reruns++) {

    /* Reset the redo-count. */
    redo = 0;

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      p = &parts[pid[i]];
      xp = &xparts[pid[i]];

      /* Is this part within the timestep? */
      if (p->ti_end <= ti_current) {

        /* Finish the density calculation */
        hydro_end_density(p, ti_current);

        /* If no derivative, double the smoothing length. */
        if (p->density.wcount_dh == 0.0f) h_corr = p->h;

        /* Otherwise, compute the smoothing length update (Newton step). */
        else {
          h_corr = (target_wcount - p->density.wcount) / p->density.wcount_dh;

          /* Truncate to the range [ -p->h/2 , p->h ]. */
          h_corr = fminf(h_corr, p->h);
          h_corr = fmaxf(h_corr, -p->h * 0.5f);
        }

        /* Did we get the right number density? */
        if (p->density.wcount > max_wcount || p->density.wcount < min_wcount) {

          /* Ok, correct then */
          p->h += h_corr;

          /* Flag for another round of fun */
          pid[redo] = pid[i];
          redo += 1;

          /* Re-initialise everything */
          hydro_init_part(p);

          /* Off we go ! */
          continue;
        }

        /* We now have a particle whose smoothing length has converged */

        /* As of here, particle force variables will be set. */

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, ti_current, timeBase);

        /* The particle force values are now set.  Do _NOT_
           try to read any particle density variables! */

        /* Prepare the particle for the force loop over neighbours */
        hydro_reset_acceleration(p);
      }
    }

    /* We now need to treat the particles whose smoothing length had not
     * converged again */

    /* Re-set the counter for the next loop (potentially). */
    count = redo;
    if (count > 0) {

      /* Climb up the cell hierarchy. */
      for (finger = c; finger != NULL; finger = finger->parent) {

        /* Run through this cell's density interactions. */
        for (struct link *l = finger->density; l != NULL; l = l->next) {

          /* Self-interaction? */
          if (l->t->type == task_type_self)
            runner_doself_subset_density(r, finger, parts, pid, count);

          /* Otherwise, pair interaction? */
          else if (l->t->type == task_type_pair) {

            /* Left or right? */
            if (l->t->ci == finger)
              runner_dopair_subset_density(r, finger, parts, pid, count,
                                           l->t->cj);
            else
              runner_dopair_subset_density(r, finger, parts, pid, count,
                                           l->t->ci);

          }

          /* Otherwise, sub interaction? */
          else if (l->t->type == task_type_sub) {

            /* Left or right? */
            if (l->t->ci == finger)
              runner_dosub_subset_density(r, finger, parts, pid, count,
                                          l->t->cj, -1, 1);
            else
              runner_dosub_subset_density(r, finger, parts, pid, count,
                                          l->t->ci, -1, 1);
          }
        }
      }
    }
  }

  if (count)
    message("Smoothing length failed to converge on %i particles.", count);

  TIMER_TOC(timer_do_ghost);
}

/**
 * @brief Drift particles and g-particles forward in time
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift(struct runner *r, struct cell *c, int timer) {

  const double timeBase = r->e->timeBase;
  const double dt = (r->e->ti_current - r->e->ti_old) * timeBase;
  const int ti_old = r->e->ti_old;
  const int ti_current = r->e->ti_current;
  struct part *const parts = c->parts;
  struct xpart *const xparts = c->xparts;
  struct gpart *const gparts = c->gparts;
  float dx_max = 0.f, dx2_max = 0.f, h_max = 0.f;

  TIMER_TIC
#ifdef TASK_VERBOSE
  OUT;
#endif

  /* No children? */
  if (!c->split) {

    /* Loop over all the g-particles in the cell */
    const int nr_gparts = c->gcount;
    for (size_t k = 0; k < nr_gparts; k++) {

      /* Get a handle on the gpart. */
      struct gpart *const gp = &gparts[k];

      /* Drift... */
      gp->x[0] += gp->v_full[0] * dt;
      gp->x[1] += gp->v_full[1] * dt;
      gp->x[2] += gp->v_full[2] * dt;

      /* Compute offset since last cell construction */
      gp->x_diff[0] -= gp->v_full[0] * dt;
      gp->x_diff[1] -= gp->v_full[1] * dt;
      gp->x_diff[2] -= gp->v_full[2] * dt;

      /* Compute (square of) motion since last cell construction */
      const float dx2 = gp->x_diff[0] * gp->x_diff[0] +
                        gp->x_diff[1] * gp->x_diff[1] +
                        gp->x_diff[2] * gp->x_diff[2];
      dx2_max = fmaxf(dx2_max, dx2);
    }

    /* Loop over all the particles in the cell (more work for these !) */
    const size_t nr_parts = c->count;
    for (size_t k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* Useful quantity */
      const float h_inv = 1.0f / p->h;

      /* Drift... */
      p->x[0] += xp->v_full[0] * dt;
      p->x[1] += xp->v_full[1] * dt;
      p->x[2] += xp->v_full[2] * dt;

      /* Predict velocities (for hydro terms) */
      p->v[0] += p->a_hydro[0] * dt;
      p->v[1] += p->a_hydro[1] * dt;
      p->v[2] += p->a_hydro[2] * dt;

      /* Predict smoothing length */
      const float w1 = p->h_dt * h_inv * dt;
      if (fabsf(w1) < 0.2f)
        p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
      else
        p->h *= expf(w1);

      /* Predict density */
      const float w2 = -3.0f * p->h_dt * h_inv * dt;
      if (fabsf(w2) < 0.2f)
        p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
      else
        p->rho *= expf(w2);

      /* Predict the values of the extra fields */
      hydro_predict_extra(p, xp, ti_old, ti_current, timeBase);

      /* Compute offset since last cell construction */
      xp->x_diff[0] -= xp->v_full[0] * dt;
      xp->x_diff[1] -= xp->v_full[1] * dt;
      xp->x_diff[2] -= xp->v_full[2] * dt;

      /* Compute (square of) motion since last cell construction */
      const float dx2 = xp->x_diff[0] * xp->x_diff[0] +
                        xp->x_diff[1] * xp->x_diff[1] +
                        xp->x_diff[2] * xp->x_diff[2];
      dx2_max = fmaxf(dx2_max, dx2);

      /* Maximal smoothing length */
      h_max = fmaxf(p->h, h_max);
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);
  }

  /* Otherwise, aggregate data from children. */
  else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_do_drift(r, cp, 0);

        dx_max = fmaxf(dx_max, cp->dx_max);
        h_max = fmaxf(h_max, cp->h_max);
      }
  }

  /* Store the values */
  c->h_max = h_max;
  c->dx_max = dx_max;

  if (timer) TIMER_TOC(timer_drift);
}

/**
 * @brief Kick particles in momentum space and collect statistics
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_kick(struct runner *r, struct cell *c, int timer) {

  const float global_dt_min = r->e->dt_min;
  const float global_dt_max = r->e->dt_max;
  const int ti_current = r->e->ti_current;
  const double timeBase = r->e->timeBase;
  const double timeBase_inv = 1.0 / r->e->timeBase;
  const int count = c->count;
  const int gcount = c->gcount;
  struct part *const parts = c->parts;
  struct xpart *const xparts = c->xparts;
  struct gpart *const gparts = c->gparts;
  const struct external_potential *potential = r->e->external_potential;
  const struct hydro_props *hydro_properties = r->e->hydro_properties;
  const struct phys_const *constants = r->e->physical_constants;
  const float ln_max_h_change = hydro_properties->log_max_h_change;
  const int is_fixdt =
      (r->e->policy & engine_policy_fixdt) == engine_policy_fixdt;

  int updated = 0, g_updated = 0;
  int ti_end_min = max_nr_timesteps, ti_end_max = 0;
  double e_kin = 0.0, e_int = 0.0, e_pot = 0.0, mass = 0.0;
  float mom[3] = {0.0f, 0.0f, 0.0f};
  float ang[3] = {0.0f, 0.0f, 0.0f};

  TIMER_TIC

#ifdef TASK_VERBOSE
  OUT;
#endif

  /* No children? */
  if (!c->split) {

    /* Loop over the g-particles and kick the active ones. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *const gp = &gparts[k];

      /* If the g-particle has no counterpart and needs to be kicked */
      if (gp->id < 0) {

        if (is_fixdt || gp->ti_end <= ti_current) {

          /* First, finish the force calculation */
          gravity_end_force(gp);

          /* Now we are ready to compute the next time-step size */
          int new_dti;

          if (is_fixdt) {

            /* Now we have a time step, proceed with the kick */
            new_dti = global_dt_max * timeBase_inv;

          } else {

            /* Compute the next timestep (gravity condition) */
            const float new_dt_external =
                gravity_compute_timestep_external(potential, constants, gp);
            const float new_dt_self =
                gravity_compute_timestep_self(constants, gp);

            float new_dt = fminf(new_dt_external, new_dt_self);

            /* Limit timestep within the allowed range */
            new_dt = fminf(new_dt, global_dt_max);
            new_dt = fmaxf(new_dt, global_dt_min);

            /* Convert to integer time */
            new_dti = new_dt * timeBase_inv;

            /* Recover the current timestep */
            const int current_dti = gp->ti_end - gp->ti_begin;

            /* Limit timestep increase */
            if (current_dti > 0) new_dti = min(new_dti, 2 * current_dti);

            /* Put this timestep on the time line */
            int dti_timeline = max_nr_timesteps;
            while (new_dti < dti_timeline) dti_timeline /= 2;

            new_dti = dti_timeline;

            /* Make sure we are allowed to increase the timestep size */
            if (new_dti > current_dti) {
              if ((max_nr_timesteps - gp->ti_end) % new_dti > 0)
                new_dti = current_dti;
            }

            /* Now we have a time step, proceed with the kick */
          }

          /* Compute the time step for this kick */
          const int ti_start = (gp->ti_begin + gp->ti_end) / 2;
          const int ti_end = gp->ti_end + new_dti / 2;
          const double dt = (ti_end - ti_start) * timeBase;
          const double half_dt = (ti_end - gp->ti_end) * timeBase;

          /* Move particle forward in time */
          gp->ti_begin = gp->ti_end;
          gp->ti_end = gp->ti_begin + new_dti;

          /* Kick particles in momentum space */
          gp->v_full[0] += gp->a_grav[0] * dt;
          gp->v_full[1] += gp->a_grav[1] * dt;
          gp->v_full[2] += gp->a_grav[2] * dt;

          /* Extra kick work */
          gravity_kick_extra(gp, dt, half_dt);

          /* Number of updated g-particles */
          g_updated++;
        }

        /* Minimal time for next end of time-step */
        ti_end_min = min(gp->ti_end, ti_end_min);
        ti_end_max = max(gp->ti_end, ti_end_max);
      }
    }

    /* Now do the hydro ones... */

    /* Loop over the particles and kick the active ones. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* If particle needs to be kicked */
      if (is_fixdt || p->ti_end <= ti_current) {

        /* First, finish the force loop */
        p->h_dt *= p->h * 0.333333333f;

        /* And do the same of the extra variable */
        hydro_end_force(p);
        if (p->gpart != NULL) gravity_end_force(p->gpart);

        /* Now we are ready to compute the next time-step size */
        int new_dti;

        if (is_fixdt) {

          /* Now we have a time step, proceed with the kick */
          new_dti = global_dt_max * timeBase_inv;

        } else {

          /* Compute the next timestep (hydro condition) */
          const float new_dt_hydro =
              hydro_compute_timestep(p, xp, hydro_properties);

          /* Compute the next timestep (gravity condition) */
          float new_dt_grav = FLT_MAX;
          if (p->gpart != NULL) {

            const float new_dt_external = gravity_compute_timestep_external(
                potential, constants, p->gpart);
            const float new_dt_self =
                gravity_compute_timestep_self(constants, p->gpart);

            new_dt_grav = fminf(new_dt_external, new_dt_self);
          }

          /* Final time-step is minimum of hydro and gravity */
          float new_dt = fminf(new_dt_hydro, new_dt_grav);

          /* Limit change in h */
          const float dt_h_change =
              (p->h_dt != 0.0f) ? fabsf(ln_max_h_change * p->h / p->h_dt)
                                : FLT_MAX;

          new_dt = fminf(new_dt, dt_h_change);

          /* Limit timestep within the allowed range */
          new_dt = fminf(new_dt, global_dt_max);
          new_dt = fmaxf(new_dt, global_dt_min);

          /* Convert to integer time */
          new_dti = new_dt * timeBase_inv;

          /* Recover the current timestep */
          const int current_dti = p->ti_end - p->ti_begin;

          /* Limit timestep increase */
          if (current_dti > 0) new_dti = min(new_dti, 2 * current_dti);

          /* Put this timestep on the time line */
          int dti_timeline = max_nr_timesteps;
          while (new_dti < dti_timeline) dti_timeline /= 2;

          new_dti = dti_timeline;

          /* Make sure we are allowed to increase the timestep size */
          if (new_dti > current_dti) {
            if ((max_nr_timesteps - p->ti_end) % new_dti > 0)
              new_dti = current_dti;
          }

          /* Now we have a time step, proceed with the kick */
        }

        /* Compute the time step for this kick */
        const int ti_start = (p->ti_begin + p->ti_end) / 2;
        const int ti_end = p->ti_end + new_dti / 2;
        const double dt = (ti_end - ti_start) * timeBase;
        const double half_dt = (ti_end - p->ti_end) * timeBase;

        /* Move particle forward in time */
        p->ti_begin = p->ti_end;
        p->ti_end = p->ti_begin + new_dti;
        if (p->gpart != NULL) {
          p->gpart->ti_begin = p->ti_begin;
          p->gpart->ti_end = p->ti_end;
        }

        /* Get the acceleration */
        float a_tot[3] = {p->a_hydro[0], p->a_hydro[1], p->a_hydro[2]};
        if (p->gpart != NULL) {
          a_tot[0] += p->gpart->a_grav[0];
          a_tot[1] += p->gpart->a_grav[1];
          a_tot[1] += p->gpart->a_grav[2];
        }

        /* Kick particles in momentum space */
        xp->v_full[0] += a_tot[0] * dt;
        xp->v_full[1] += a_tot[1] * dt;
        xp->v_full[2] += a_tot[2] * dt;

        if (p->gpart != NULL) {
          p->gpart->v_full[0] = xp->v_full[0];
          p->gpart->v_full[1] = xp->v_full[1];
          p->gpart->v_full[2] = xp->v_full[2];
        }

        /* Go back by half-step for the hydro velocity */
        p->v[0] = xp->v_full[0] - half_dt * a_tot[0];
        p->v[1] = xp->v_full[1] - half_dt * a_tot[1];
        p->v[2] = xp->v_full[2] - half_dt * a_tot[2];

        /* Extra kick work */
        hydro_kick_extra(p, xp, dt, half_dt);
        if (p->gpart != NULL) gravity_kick_extra(p->gpart, dt, half_dt);

        /* Number of updated particles */
        updated++;
        if (p->gpart != NULL) g_updated++;
      }

      /* Now collect quantities for statistics */

      const double x[3] = {p->x[0], p->x[1], p->x[2]};
      const float v_full[3] = {xp->v_full[0], xp->v_full[1], xp->v_full[2]};
      const float m = p->mass;

      /* Collect mass */
      mass += m;

      /* Collect momentum */
      mom[0] += m * v_full[0];
      mom[1] += m * v_full[1];
      mom[2] += m * v_full[2];

      /* Collect angular momentum */
      ang[0] += m * (x[1] * v_full[2] - x[2] * v_full[1]);
      ang[1] += m * (x[2] * v_full[0] - x[0] * v_full[2]);
      ang[2] += m * (x[0] * v_full[1] - x[1] * v_full[0]);

      /* Collect total energy. */
      e_kin += 0.5 * m * (v_full[0] * v_full[0] + v_full[1] * v_full[1] +
                          v_full[2] * v_full[2]);
      e_pot += 0.f; /* No gravitational potential thus far */
      e_int += hydro_get_internal_energy(p);

      /* Minimal time for next end of time-step */
      ti_end_min = min(p->ti_end, ti_end_min);
      ti_end_max = max(p->ti_end, ti_end_max);
    }
  }

  /* Otherwise, aggregate data from children. */
  else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *const cp = c->progeny[k];

        /* Recurse */
        runner_do_kick(r, cp, 0);

        /* And aggregate */
        updated += cp->updated;
        g_updated += cp->g_updated;
        e_kin += cp->e_kin;
        e_int += cp->e_int;
        e_pot += cp->e_pot;
        mass += cp->mass;
        mom[0] += cp->mom[0];
        mom[1] += cp->mom[1];
        mom[2] += cp->mom[2];
        ang[0] += cp->ang[0];
        ang[1] += cp->ang[1];
        ang[2] += cp->ang[2];
        ti_end_min = min(cp->ti_end_min, ti_end_min);
        ti_end_max = max(cp->ti_end_max, ti_end_max);
      }
  }

  /* Store the values. */
  c->updated = updated;
  c->g_updated = g_updated;
  c->e_kin = e_kin;
  c->e_int = e_int;
  c->e_pot = e_pot;
  c->mass = mass;
  c->mom[0] = mom[0];
  c->mom[1] = mom[1];
  c->mom[2] = mom[2];
  c->ang[0] = ang[0];
  c->ang[1] = ang[1];
  c->ang[2] = ang[2];
  c->ti_end_min = ti_end_min;
  c->ti_end_max = ti_end_max;

  if (timer) TIMER_TOC(timer_kick);
}

/**
 * @brief Construct the cell properties from the received particles
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_recv_cell(struct runner *r, struct cell *c, int timer) {

  const struct part *const parts = c->parts;
  const struct gpart *const gparts = c->gparts;
  const size_t nr_parts = c->count;
  const size_t nr_gparts = c->gcount;
  // const int ti_current = r->e->ti_current;

  TIMER_TIC;

  int ti_end_min = max_nr_timesteps;
  int ti_end_max = 0;
  float h_max = 0.f;

  /* Collect everything... */
  for (size_t k = 0; k < nr_parts; k++) {
    const int ti_end = parts[k].ti_end;
    // if(ti_end < ti_current) error("Received invalid particle !");
    ti_end_min = min(ti_end_min, ti_end);
    ti_end_max = max(ti_end_max, ti_end);
    h_max = fmaxf(h_max, parts[k].h);
  }
  for (size_t k = 0; k < nr_gparts; k++) {
    const int ti_end = gparts[k].ti_end;
    // if(ti_end < ti_current) error("Received invalid particle !");
    ti_end_min = min(ti_end_min, ti_end);
    ti_end_max = max(ti_end_max, ti_end);
  }

  /* ... and store. */
  c->ti_end_min = ti_end_min;
  c->ti_end_max = ti_end_max;
  c->h_max = h_max;

  if (timer) TIMER_TOC(timer_dorecv_cell);
}

/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 */
void *runner_main(void *data) {

  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;

  /* Main loop. */
  while (1) {

    /* Wait at the barrier. */
    engine_barrier(e, r->id);

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;

    /* Loop while there are tasks... */
    while (1) {

      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, prev);
        TIMER_TOC(timer_gettask);

        /* Did I get anything? */
        if (t == NULL) break;
      }

      /* Get the cells. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;
      t->rid = r->cpuid;
      t->last_rid = r->cpuid;

      /* Different types of tasks... */
      switch (t->type) {
        case task_type_self:
          if (t->subtype == task_subtype_density)
            runner_doself1_density(r, ci);
          else if (t->subtype == task_subtype_force)
            runner_doself2_force(r, ci);
          else
            error("Unknown task subtype.");
          break;
        case task_type_pair:
          if (t->subtype == task_subtype_density)
            runner_dopair1_density(r, ci, cj);
          else if (t->subtype == task_subtype_force)
            runner_dopair2_force(r, ci, cj);
          else
            error("Unknown task subtype.");
          break;
        case task_type_sort:
          runner_do_sort(r, ci, t->flags, 1);
          break;
        case task_type_sub:
          if (t->subtype == task_subtype_density)
            runner_dosub1_density(r, ci, cj, t->flags, 1);
          else if (t->subtype == task_subtype_force)
            runner_dosub2_force(r, ci, cj, t->flags, 1);
          else if (t->subtype == task_subtype_grav)
            runner_dosub_grav(r, ci, cj, 1);
          else
            error("Unknown task subtype.");
          break;
        case task_type_init:
          runner_do_init(r, ci, 1);
          break;
        case task_type_ghost:
          runner_do_ghost(r, ci);
          break;
        case task_type_drift:
          runner_do_drift(r, ci, 1);
          break;
        case task_type_kick:
          runner_do_kick(r, ci, 1);
          break;
        case task_type_send:
          break;
        case task_type_recv:
          runner_do_recv_cell(r, ci, 1);
          break;
        case task_type_grav_pp:
          if (t->cj == NULL)
            runner_doself_grav(r, t->ci);
          else
            runner_dopair_grav(r, t->ci, t->cj);
          break;
        case task_type_grav_mm:
          runner_dograv_mm(r, t->ci, t->cj);
          break;
        case task_type_grav_up:
          runner_dograv_up(r, t->ci);
          break;
        case task_type_grav_down:
          runner_dograv_down(r, t->ci);
          break;
        case task_type_grav_external:
          runner_do_grav_external(r, t->ci, 1);
          break;
        case task_type_part_sort:
          space_do_parts_sort();
          break;
        case task_type_gpart_sort:
          space_do_gparts_sort();
          break;
        case task_type_split_cell:
          space_do_split(e->s, t->ci);
          break;
        case task_type_rewait:
          scheduler_do_rewait((struct task *)t->ci, (struct task *)t->cj,
                              t->flags, t->rank);
          break;
        default:
          error("Unknown task type.");
      }

      /* We're done with this task, see if we get a next one. */
      prev = t;
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
