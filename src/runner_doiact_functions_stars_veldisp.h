/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2020 Joel Pfeffer (j.l.pfeffer@ljmu.ac.uk)
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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_stars_veldisp.h"

/**
 * @brief Calculate the number density of #spart around the #spart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_SVD(struct runner *r, struct cell *c, int timer) {
#ifdef STARS_MOSAICS

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_star_formation = (e->policy & engine_policy_star_formation);

  /* Anything to do here? */
  if (c->stars.count == 0) return;
  if (!(with_star_formation && cell_is_active_hydro(c, e))) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  struct spart *restrict sparts = c->stars.parts;

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts[sid];

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    /* Only want the star particles formed in this timestep */
    if (!si->new_star) continue;

    const float hi = si->hbirth;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the sparts in cj. */
    for (int sjd = 0; sjd < scount; sjd++) {

      /* Get a pointer to the jth particle. */
      struct spart *restrict sj = &sparts[sjd];

      /* Early abort? */
      if (spart_is_inhibited(sj, e)) continue;

      /* This particle was already included as gas in the hydro loop */
      if (sj->new_star) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - c->loc[0]),
                            (float)(sj->x[1] - c->loc[1]),
                            (float)(sj->x[2] - c->loc[2])};
      float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      if (r2 < hig2) {
        runner_iact_nonsym_star_veldisp(r2, dx, hi, si, sj, a, H);
      }
    } /* loop over the sparts in ci. */
  }   /* loop over the sparts in ci. */

  TIMER_TOC(TIMER_DOSELF_SVD);
#endif
}

/**
 * @brief Calculate the number density of cj #spart around the ci #spart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_SVD_NAIVE(struct runner *r, struct cell *restrict ci,
                                 struct cell *restrict cj) {
#ifdef STARS_MOSAICS

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_star_formation = (e->policy & engine_policy_star_formation);

  /* Anything to do here? */
  if (ci->stars.count == 0 || cj->stars.count == 0) return;
  if (!(with_star_formation && cell_is_active_hydro(ci, e))) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->stars.count;
  const int scount_j = cj->stars.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct spart *restrict sparts_j = cj->stars.parts;

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

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    /* Only want the star particles formed in this timestep */
    if (!si->new_star) continue;

    const float hi = si->hbirth;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the sparts in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct spart *restrict sj = &sparts_j[sjd];

      /* Skip inhibited particles. */
      if (spart_is_inhibited(sj, e)) continue;

      /* This particle was already included as gas in the hydro loop */
      if (sj->new_star) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - cj->loc[0]),
                            (float)(sj->x[1] - cj->loc[1]),
                            (float)(sj->x[2] - cj->loc[2])};
      float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      if (r2 < hig2) {
        runner_iact_nonsym_star_veldisp(r2, dx, hi, si, sj, a, H);
      }
    } /* loop over the sparts in cj. */
  }   /* loop over the sparts in ci. */
#endif
}

void DOPAIR1_SVD_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct cell *restrict cj, int timer) {
#ifdef STARS_MOSAICS

  TIMER_TIC;

  const int do_ci_stars = ci->nodeID == r->e->nodeID;
  const int do_cj_stars = cj->nodeID == r->e->nodeID;
  if (do_ci_stars && ci->stars.count != 0 && cj->stars.count != 0)
    DO_NONSYM_PAIR1_SVD_NAIVE(r, ci, cj);
  if (do_cj_stars && cj->stars.count != 0 && ci->stars.count != 0)
    DO_NONSYM_PAIR1_SVD_NAIVE(r, cj, ci);

  TIMER_TOC(TIMER_DOPAIR_SVD);
#endif
}

/**
 * @brief Determine which version of DOSELF1_SVD needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_SVD(struct runner *r, struct cell *c) {
#ifdef STARS_MOSAICS

  const struct engine *restrict e = r->e;
  const int with_star_formation = (e->policy & engine_policy_star_formation);

  /* Anything to do here? */
  if (c->stars.count == 0) return;

  /* Anything to do here? */
  if (!(with_star_formation && cell_is_active_hydro(c, e))) return;

  /* Did we mess up the recursion? */
  if (c->hydro.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  DOSELF1_SVD(r, c, 1);
#endif
}

/**
 * @brief Determine which version of DOPAIR1_SVD needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_SVD needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_SVD(struct runner *r, struct cell *ci, struct cell *cj) {
#ifdef STARS_MOSAICS

  const struct engine *restrict e = r->e;
  const int with_star_formation = (e->policy & engine_policy_star_formation);

  const int ci_active = with_star_formation && cell_is_active_hydro(ci, e);
  const int cj_active = with_star_formation && cell_is_active_hydro(cj, e);
  const int do_ci_stars = ci->nodeID == e->nodeID;
  const int do_cj_stars = cj->nodeID == e->nodeID;
  const int do_ci = (ci->stars.count != 0 && cj->stars.count != 0 &&
                     ci_active && do_ci_stars);
  const int do_cj = (cj->stars.count != 0 && ci->stars.count != 0 &&
                     cj_active && do_cj_stars);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci && 
      (!cell_are_spart_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_spart_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  DOPAIR1_SVD_NAIVE(r, ci, cj, 1);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1_SVD(struct runner *r, struct cell *ci, struct cell *cj,
                       int gettimer) {
#ifdef STARS_MOSAICS

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;
  const int with_star_formation = (e->policy & engine_policy_star_formation);

  /* Should we even bother? */
  const int should_do_ci = ci->stars.count != 0 && cj->stars.count != 0 &&
      (with_star_formation && cell_is_active_hydro(ci, e));
  const int should_do_cj = cj->stars.count != 0 && ci->stars.count != 0 &&
      (with_star_formation && cell_is_active_hydro(cj, e));
  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_hydro_task(ci) &&
      cell_can_recurse_in_pair_hydro_task(cj)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_SVD(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci_stars = ci->nodeID == e->nodeID;
    const int do_cj_stars = cj->nodeID == e->nodeID;
    const int do_ci = should_do_ci && do_ci_stars;
    const int do_cj = should_do_cj && do_cj_stars;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (!cell_are_spart_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (!cell_are_spart_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_SVD(r, ci, cj);
  }

  TIMER_TOC(TIMER_DOSUB_PAIR_SVD);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_SVD(struct runner *r, struct cell *ci, int gettimer) {
#ifdef STARS_MOSAICS

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  const struct engine *e = r->e;
  const int with_star_formation = (e->policy & engine_policy_star_formation);

  /* Should we even bother? */
  if (ci->stars.count == 0 || 
      !(with_star_formation && cell_is_active_hydro(ci, r->e)))
    return;

  /* Recurse? */
  if (cell_can_recurse_in_self_hydro_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_SVD(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_SVD(r, ci->progeny[k], ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_spart_drifted(ci, r->e))
        error("Interacting undrifted cell.");

    DOSELF1_BRANCH_SVD(r, ci);
  }

  TIMER_TOC(TIMER_DOSUB_SELF_SVD);
#endif
}
