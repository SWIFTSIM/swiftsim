/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_RUNNER_DOIACT_STARS_H
#define SWIFT_RUNNER_DOIACT_STARS_H

#include "swift.h"

/**
 * @brief Recursively checks that the flags are consistent in a cell hierarchy.
 *
 * Debugging function.
 *
 * @param c The #cell to check.
 * @param flags The sorting flags to check.
 */
void runner_check_stars_sorts(struct cell *c, int flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (flags & ~c->stars.sorted) error("Inconsistent sort flags (downward)!");
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0)
        runner_check_stars_sorts(c->progeny[k], c->stars.sorted);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Sort the stars particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param cleanup If true, re-build the sorts for the selected flags instead
 *        of just adding them.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */
void runner_do_stars_sort(struct runner *r, struct cell *c, int flags, int cleanup,
                    int clock) {

  struct entry *fingers[8];
  const int count = c->stars.count;
  struct spart *parts = c->stars.parts;
  float buff[8];

  TIMER_TIC;

  /* We need to do the local sorts plus whatever was requested further up. */
  flags |= c->stars.do_sort;
  if (cleanup) {
    c->stars.sorted = 0;
  } else {
    flags &= ~c->stars.sorted;
  }
  if (flags == 0 && !c->stars.do_sub_sort) return;

  /* Check that the particles have been moved to the current time */
  if (flags && !cell_are_gpart_drifted(c, r->e))
    error("Sorting un-drifted cell c->nodeID=%d", c->nodeID);

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure the sort flags are consistent (downward). */
  runner_check_stars_sorts(c, c->stars.sorted);

  /* Make sure the sort flags are consistent (upard). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->stars.sorted & ~c->stars.sorted)
      error("Inconsistent sort flags (upward).");
  }

  /* Update the sort timer which represents the last time the sorts
     were re-set. */
  if (c->stars.sorted == 0) c->stars.ti_sort = r->e->ti_current;
#endif

  /* start by allocating the entry arrays in the requested dimensions. */
  for (int j = 0; j < 13; j++) {
    if ((flags & (1 << j)) && c->stars.sort[j] == NULL) {
      if ((c->stars.sort[j] = (struct entry *)malloc(sizeof(struct entry) *
                                                     (count + 1))) == NULL)
        error("Failed to allocate sort memory.");
    }
  }

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    float dx_max_sort = 0.0f;
    float dx_max_sort_old = 0.0f;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0) {
        /* Only propagate cleanup if the progeny is stale. */
        runner_do_stars_sort(r, c->progeny[k], flags,
                       cleanup && (c->progeny[k]->stars.dx_max_sort >
                                   space_maxreldx * c->progeny[k]->dmin),
                       0);
        dx_max_sort = max(dx_max_sort, c->progeny[k]->stars.dx_max_sort);
        dx_max_sort_old =
            max(dx_max_sort_old, c->progeny[k]->stars.dx_max_sort_old);
      }
    }
    c->stars.dx_max_sort = dx_max_sort;
    c->stars.dx_max_sort_old = dx_max_sort_old;

    /* Loop over the 13 different sort arrays. */
    for (int j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      int off[8];
      off[0] = 0;
      for (int k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->stars.count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      int inds[8];
      for (int k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0) {
          fingers[k] = c->progeny[k]->stars.sort[j];
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (int i = 0; i < 7; i++)
        for (int k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            int temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      struct entry *finger = c->stars.sort[j];
      for (int ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (int k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          int temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      c->stars.sort[j][count].d = FLT_MAX;
      c->stars.sort[j][count].i = 0;

      /* Mark as sorted. */
      atomic_or(&c->stars.sorted, 1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Reset the sort distance */
    if (c->stars.sorted == 0) {

      /* And the individual sort distances if we are a local cell */
      for (int k = 0; k < count; k++) {
	parts[k].x_diff_sort[0] = 0.0f;
	parts[k].x_diff_sort[1] = 0.0f;
        parts[k].x_diff_sort[2] = 0.0f;
      }
      c->stars.dx_max_sort_old = 0.f;
      c->stars.dx_max_sort = 0.f;
    }

    /* Fill the sort array. */
    for (int k = 0; k < count; k++) {
      const double px[3] = {parts[k].x[0], parts[k].x[1], parts[k].x[2]};
      for (int j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          c->stars.sort[j][k].i = k;
          c->stars.sort[j][k].d = px[0] * runner_shift[j][0] +
                                  px[1] * runner_shift[j][1] +
                                  px[2] * runner_shift[j][2];
        }
    }

    /* Add the sentinel and sort. */
    for (int j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        c->stars.sort[j][count].d = FLT_MAX;
        c->stars.sort[j][count].i = 0;
        runner_do_sort_ascending(c->stars.sort[j], count);
        atomic_or(&c->stars.sorted, 1 << j);
      }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify the sorting. */
  for (int j = 0; j < 13; j++) {
    if (!(flags & (1 << j))) continue;
    struct entry *finger = c->stars.sort[j];
    for (int k = 1; k < count; k++) {
      if (finger[k].d < finger[k - 1].d)
        error("Sorting failed, ascending array.");
      if (finger[k].i >= count) error("Sorting failed, indices borked.");
    }
  }

  /* Make sure the sort flags are consistent (downward). */
  runner_check_stars_sorts(c, flags);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->stars.sorted & ~c->stars.sorted)
      error("Inconsistent sort flags.");
  }
#endif

  /* Clear the cell's sort flags. */
  c->stars.do_sort = 0;
  c->stars.do_sub_sort = 0;
  c->stars.requires_sorts = 0;

  if (clock) TIMER_TOC(timer_dostars_sort);
}


/**
 * @brief Calculate the number density of #part around the #spart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_doself_stars_density(struct runner *r, struct cell *c, int timer) {
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  const int count = c->hydro.count;
  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts[sid];
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      if (r2 > 0.f && r2 < hig2) {
        runner_iact_nonsym_stars_density(r2, dx, hi, hj, si, pj, a, H);
      }
    } /* loop over the parts in ci. */
  }   /* loop over the sparts in ci. */

  TIMER_TOC(timer_doself_stars_density);
}

/**
 * @brief Calculate the number density of cj #part around the ci #spart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void runner_dosubpair_stars_density(struct runner *r, struct cell *restrict ci,
                                    struct cell *restrict cj) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (!cell_is_active_stars(ci, e) && !cell_is_active_stars(cj, e)) return;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts_i[sid];
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      if (r2 < hig2)
        runner_iact_nonsym_stars_density(r2, dx, hi, hj, si, pj, a, H);

    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

void runner_dopair_stars_density(struct runner *r, struct cell *restrict ci,
                                 struct cell *restrict cj, int timer) {

  TIMER_TIC;

  runner_dosubpair_stars_density(r, ci, cj);
  runner_dosubpair_stars_density(r, cj, ci);

  if (timer) TIMER_TOC(timer_dopair_stars_density);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void runner_dopair_subset_stars_density(struct runner *r,
                                        struct cell *restrict ci,
                                        struct spart *restrict sparts_i,
                                        int *restrict ind, int scount,
                                        struct cell *restrict cj,
                                        const double *shift) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over the parts_i. */
  for (int pid = 0; pid < scount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *restrict spi = &sparts_i[ind[pid]];
    double spix[3];
    for (int k = 0; k < 3; k++) spix[k] = spi->x[k] - shift[k];
    const float hi = spi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(spi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = spix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        runner_iact_nonsym_stars_density(r2, dx, hi, pj->h, spi, pj, a, H);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(timer_dopair_subset_naive);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts The #spart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void runner_doself_subset_stars_density(struct runner *r,
                                        struct cell *restrict ci,
                                        struct spart *restrict sparts,
                                        int *restrict ind, int scount) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;

  /* Loop over the parts in ci. */
  for (int spid = 0; spid < scount; spid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *spi = &sparts[ind[spid]];
    const float spix[3] = {(float)(spi->x[0] - ci->loc[0]),
                           (float)(spi->x[1] - ci->loc[1]),
                           (float)(spi->x[2] - ci->loc[2])};
    const float hi = spi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(spi, e))
      error("Inactive particle in subset function!");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                            (float)(pj->x[1] - ci->loc[1]),
                            (float)(pj->x[2] - ci->loc[2])};
      float dx[3] = {spix[0] - pjx[0], spix[1] - pjx[1], spix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 > 0.f && r2 < hig2) {
        runner_iact_nonsym_stars_density(r2, dx, hi, hj, spi, pj, a, H);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(timer_doself_subset_stars_density);
}

/**
* @brief Determine which version of DOSELF_SUBSET needs to be called depending
* on the optimisation level.

* @param r The #runner.
* @param ci The first #cell.
* @param sparts The #spart to interact.
* @param ind The list of indices of particles in @c ci to interact with.
* @param scount The number of particles in @c ind.
*/
void runner_doself_subset_branch_stars_density(struct runner *r,
                                               struct cell *restrict ci,
                                               struct spart *restrict sparts,
                                               int *restrict ind, int scount) {

  runner_doself_subset_stars_density(r, ci, sparts, ind, scount);
}

/**
 * @brief Determine which version of DOPAIR_SUBSET needs to be called depending
 * on the
 * orientation of the cells or whether DOPAIR_SUBSET needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #spart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void runner_dopair_subset_branch_stars_density(struct runner *r,
                                               struct cell *restrict ci,
                                               struct spart *restrict sparts_i,
                                               int *restrict ind, int scount,
                                               struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  runner_dopair_subset_stars_density(r, ci, sparts_i, ind, scount, cj, shift);
}

void runner_dosub_subset_stars_density(struct runner *r, struct cell *ci,
                                       struct spart *sparts, int *ind,
                                       int scount, struct cell *cj, int sid,
                                       int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (!cell_is_active_stars(ci, e) &&
      (cj == NULL || !cell_is_active_stars(cj, e)))
    return;
  if (ci->stars.count == 0 || (cj != NULL && cj->stars.count == 0)) return;

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&sparts[ind[0]] >= &ci->progeny[k]->stars.parts[0] &&
            &sparts[ind[0]] <
                &ci->progeny[k]->stars.parts[ci->progeny[k]->stars.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_stars_task(ci)) {

      /* Loop over all progeny. */
      runner_dosub_subset_stars_density(r, sub, sparts, ind, scount, NULL, -1,
                                        0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          runner_dosub_subset_stars_density(r, sub, sparts, ind, scount,
                                            ci->progeny[j], -1, 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      runner_doself_subset_branch_stars_density(r, ci, sparts, ind, scount);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_stars_task(ci) &&
        cell_can_recurse_in_pair_stars_task(cj)) {

      /* Get the type of pair if not specified explicitly. */
      double shift[3] = {0.0, 0.0, 0.0};
      sid = space_getsid(s, &ci, &cj, shift);

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[4], sparts, ind,
                                              scount, cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[3], sparts, ind,
                                              scount, ci->progeny[4], -1, 0);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[2], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[2], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[2], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[5] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[2], sparts, ind,
                                              scount, cj->progeny[5], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[5], sparts, ind,
                                              scount, ci->progeny[2], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[5] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[5], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[5] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[5], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[5] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[5], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[5] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[5], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[5] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[5], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[2], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[5] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[2], sparts, ind,
                                              scount, cj->progeny[5], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[5], sparts, ind,
                                              scount, ci->progeny[2], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[1], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[5] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[6], sparts, ind,
                                              scount, cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[5], sparts, ind,
                                              scount, ci->progeny[6], -1, 0);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[1], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[1], -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[1], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[1], -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[1], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[1], -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[6] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[1], sparts, ind,
                                              scount, cj->progeny[6], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[6] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[6], sparts, ind,
                                              scount, ci->progeny[1], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[6] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[3], sparts, ind,
                                              scount, cj->progeny[6], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[6] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[6], sparts, ind,
                                              scount, ci->progeny[3], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[6] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[5], sparts, ind,
                                              scount, cj->progeny[6], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[6] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[6], sparts, ind,
                                              scount, ci->progeny[5], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[0], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[2], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[4], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[6] != NULL)
            runner_dosub_subset_stars_density(r, ci->progeny[7], sparts, ind,
                                              scount, cj->progeny[6], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[6] == sub)
            runner_dosub_subset_stars_density(r, cj->progeny[6], sparts, ind,
                                              scount, ci->progeny[7], -1, 0);
          break;
      }

    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_stars(ci, e) || cell_is_active_stars(cj, e)) {

      /* Do any of the cells need to be drifted first? */
      if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");

      runner_dopair_subset_branch_stars_density(r, ci, sparts, ind, scount, cj);
    }

  } /* otherwise, pair interaction. */

  if (gettimer) TIMER_TOC(timer_dosub_subset);
}

/**
 * @brief Determine which version of runner_doself needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void runner_doself_branch_stars_density(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->stars.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  runner_doself_stars_density(r, c, 1);
}

/**
 * @brief Determine which version of DOPAIR1 needs to be called depending on the
 * orientation of the cells or whether DOPAIR1 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void runner_dopair_branch_stars_density(struct runner *r, struct cell *ci,
                                        struct cell *cj) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (!cell_is_active_stars(ci, e) && !cell_is_active_stars(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_i = ci->hydro.sort[sid];
  const struct entry *restrict sort_j = cj->hydro.sort[sid];

  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->hydro.count; pid++) {
    const struct part *p = &ci->hydro.parts[sort_i[pid].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->hydro.dx_max_sort >
            1.0e-4 * max(fabsf(d), ci->hydro.dx_max_sort_old) &&
        fabsf(d - sort_i[pid].d) - ci->hydro.dx_max_sort >
            ci->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->hydro.dx_max_sort=%e "
          "ci->hydro.dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->hydro.dx_max_sort,
          ci->hydro.dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->hydro.count; pjd++) {
    const struct part *p = &cj->hydro.parts[sort_j[pjd].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if ((fabsf(d - sort_j[pjd].d) - cj->hydro.dx_max_sort) >
            1.0e-4 * max(fabsf(d), cj->hydro.dx_max_sort_old) &&
        (fabsf(d - sort_j[pjd].d) - cj->hydro.dx_max_sort) >
            cj->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->hydro.dx_max_sort=%e "
          "cj->hydro.dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->hydro.dx_max_sort,
          cj->hydro.dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

  runner_dopair_stars_density(r, ci, cj, 1);
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction linking the cells
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void runner_dosub_pair_stars_density(struct runner *r, struct cell *ci,
                                     struct cell *cj, int sid, int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  TIMER_TIC;

  /* Should we even bother? */
  if (!cell_is_active_stars(ci, e) && !cell_is_active_stars(cj, e)) return;
  if (ci->stars.count == 0 || cj->stars.count == 0) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_stars_task(ci) &&
      cell_can_recurse_in_pair_stars_task(cj)) {

    /* Different types of flags. */
    switch (sid) {

      /* Regular sub-cell interactions of a single cell. */
      case 0: /* (  1 ,  1 ,  1 ) */
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        break;

      case 1: /* (  1 ,  1 ,  0 ) */
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[1], -1,
                                          0);
        break;

      case 2: /* (  1 ,  1 , -1 ) */
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[1], -1,
                                          0);
        break;

      case 3: /* (  1 ,  0 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[2], -1,
                                          0);
        break;

      case 4: /* (  1 ,  0 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[3], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[3], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[3], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[3], -1,
                                          0);
        break;

      case 5: /* (  1 ,  0 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[3], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[3], -1,
                                          0);
        break;

      case 6: /* (  1 , -1 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[2], -1,
                                          0);
        break;

      case 7: /* (  1 , -1 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[3], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[3], -1,
                                          0);
        break;

      case 8: /* (  1 , -1 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[4], cj->progeny[3], -1,
                                          0);
        break;

      case 9: /* (  0 ,  1 ,  1 ) */
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[4], -1,
                                          0);
        break;

      case 10: /* (  0 ,  1 ,  0 ) */
        if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[2], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[2], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[2], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[2], cj->progeny[5], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[5], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[5], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[5], -1,
                                          0);
        break;

      case 11: /* (  0 ,  1 , -1 ) */
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[2], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[2], cj->progeny[5], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[1], -1,
                                          0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[6], cj->progeny[5], -1,
                                          0);
        break;

      case 12: /* (  0 ,  0 ,  1 ) */
        if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[1], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[1], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[1], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[1], cj->progeny[6], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[3], cj->progeny[6], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[5], cj->progeny[6], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[0], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[2], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[4], -1,
                                          0);
        if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
          runner_dosub_pair_stars_density(r, ci->progeny[7], cj->progeny[6], -1,
                                          0);
        break;
    }

  }

  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_stars(ci, e) || cell_is_active_stars(cj, e)) {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
      error("Interacting undrifted cells.");

    /* Do any of the cells need to be sorted first? */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx)
      error("Interacting unsorted cell.");
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx)
      error("Interacting unsorted cell.");

    /* Compute the interactions. */
    runner_dopair_branch_stars_density(r, ci, cj);
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void runner_dosub_self_stars_density(struct runner *r, struct cell *ci,
                                     int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  if (ci->stars.count == 0 || !cell_is_active_stars(ci, r->e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_stars_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        runner_dosub_self_stars_density(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            runner_dosub_pair_stars_density(r, ci->progeny[k], ci->progeny[j],
                                            -1, 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    runner_doself_branch_stars_density(r, ci);
  }

  if (gettimer) TIMER_TOC(timer_dosub_self_stars_density);
}

#endif // SWIFT_RUNNER_DOIACT_STARS_H
