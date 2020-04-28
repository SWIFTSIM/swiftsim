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

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (c->stars.count == 0) return;
  if (!cell_is_active_stars(c, e) && !cell_is_active_hydro(c, e)) return;

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

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (ci->stars.count == 0 || cj->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e) && !cell_is_active_hydro(ci, e)) return;

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
}

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DO_SYM_PAIR1_SVD(struct runner *r, struct cell *ci, struct cell *cj,
                        const int sid, const double *shift) {

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  const int do_ci_stars = (ci->nodeID == e->nodeID) && (ci->stars.count != 0) &&
                          (cj->stars.count != 0);
  const int do_cj_stars = (cj->nodeID == e->nodeID) && (cj->stars.count != 0) &&
                          (ci->stars.count != 0);

  if (do_ci_stars) {

    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_j = cell_get_stars_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->stars.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->stars.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->stars.dx_max_part, cj->stars.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const int count_i = ci->stars.count;
    const int count_j = cj->stars.count;
    struct spart *restrict sparts_i = ci->stars.parts;
    struct spart *restrict sparts_j = cj->stars.parts;
    const double dj_min = sort_j[0].d;
    const float dx_max_rshift = cj->stars.dx_max_sort - rshift;

    /* Loop over the sparts in ci. */
    /* Naive loop is good enough here since the new stars are few (and we don't
     * have hi_max) */
    for (int pid = 0; pid < count_i; pid++) {

      /* Get a hold of the ith part in ci. */
      struct spart *restrict spi = &sparts_i[pid];
      const float hi = spi->hbirth;

      /* Skip inactive particles */
      if (!spart_is_active(spi, e)) continue;

      /* Only want the star particles formed in this timestep */
      if (!spi->new_star) continue;

      /* Compute distance from the other cell. */
      const double px[3] = {spi->x[0], spi->x[1], spi->x[2]};
      float dist = px[0] * runner_shift[sid][0] + px[1] * runner_shift[sid][1] +
                   px[2] * runner_shift[sid][2];

      /* Is there anything we need to interact with ? */
      const double di = dist + hi * kernel_gamma + dx_max_rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float hig2 = hi * hi * kernel_gamma2;
      const float pix = spi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = spi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = spi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct spart *spj = &sparts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (spart_is_inhibited(spj, e)) continue;

        /* This particle was already included as gas in the hydro loop */
        if (spj->new_star) continue;

        const float pjx = spj->x[0] - cj->loc[0];
        const float pjy = spj->x[1] - cj->loc[1];
        const float pjz = spj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle spi not drifted to current time");
        if (spj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          runner_iact_nonsym_star_veldisp(r2, dx, hi, spi, spj, a, H);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* do_ci_stars */

  if (do_cj_stars) {
    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = cell_get_stars_sorts(ci, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->stars.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->stars.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->stars.dx_max_part, cj->stars.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const int count_i = ci->stars.count;
    const int count_j = cj->stars.count;
    struct spart *restrict sparts_i = ci->stars.parts;
    struct spart *restrict sparts_j = cj->stars.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const float dx_max_rshift = ci->stars.dx_max_sort - rshift;

    /* Loop over the parts in cj. */
    /* Naive loop is good enough here since the new stars are few and we don't
     * have hj_max */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a hold of the jth part in cj. */
      struct spart *spj = &sparts_j[pjd];
      const float hj = spj->hbirth;

      /* Skip inactive particles */
      if (!spart_is_active(spj, e)) continue;

      /* Only want the star particles formed in this timestep */
      if (!spj->new_star) continue;

      /* Compute distance from the other cell. */
      const double px[3] = {spj->x[0], spj->x[1], spj->x[2]};
      float dist = px[0] * runner_shift[sid][0] + px[1] * runner_shift[sid][1] +
                   px[2] * runner_shift[sid][2];

      /* Is there anything we need to interact with ? */
      const double dj = dist - hj * kernel_gamma - dx_max_rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float hjg2 = hj * hj * kernel_gamma2;
      const float pjx = spj->x[0] - cj->loc[0];
      const float pjy = spj->x[1] - cj->loc[1];
      const float pjz = spj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct spart *spi = &sparts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (spart_is_inhibited(spi, e)) continue;

        /* This particle was already included as gas in the hydro loop */
        if (spi->new_star) continue;

        const float pix = spi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = spi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = spi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (spj->ti_drift != e->ti_current)
          error("Particle spj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {
          runner_iact_nonsym_star_veldisp(r2, dx, hj, spj, spi, a, H);
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR_SVD);
}

void DOPAIR1_SVD_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct cell *restrict cj, int timer) {

  TIMER_TIC;

  const int do_ci_stars = ci->nodeID == r->e->nodeID;
  const int do_cj_stars = cj->nodeID == r->e->nodeID;
  if (do_ci_stars && ci->stars.count != 0 && cj->stars.count != 0)
    DO_NONSYM_PAIR1_SVD_NAIVE(r, ci, cj);
  if (do_cj_stars && cj->stars.count != 0 && ci->stars.count != 0)
    DO_NONSYM_PAIR1_SVD_NAIVE(r, cj, ci);

  TIMER_TOC(TIMER_DOPAIR_SVD);
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

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->stars.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e) && !cell_is_active_hydro(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->stars.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  DOSELF1_SVD(r, c, 1);
}

#define RUNNER_CHECK_SORT(TYPE, PART, cj, ci, sid)                          \
  ({                                                                        \
    const struct sort_entry *restrict sort_j =                              \
        cell_get_##TYPE##_sorts(cj, sid);                                   \
                                                                            \
    for (int pjd = 0; pjd < cj->TYPE.count; pjd++) {                        \
      const struct PART *p = &cj->TYPE.parts[sort_j[pjd].i];                \
      if (PART##_is_inhibited(p, e)) continue;                              \
                                                                            \
      const float d = p->x[0] * runner_shift[sid][0] +                      \
                      p->x[1] * runner_shift[sid][1] +                      \
                      p->x[2] * runner_shift[sid][2];                       \
      if ((fabsf(d - sort_j[pjd].d) - cj->TYPE.dx_max_sort) >               \
              1.0e-4 * max(fabsf(d), cj->TYPE.dx_max_sort_old) &&           \
          (fabsf(d - sort_j[pjd].d) - cj->TYPE.dx_max_sort) >               \
              cj->width[0] * 1.0e-10)                                       \
        error(                                                              \
            "particle shift diff exceeds dx_max_sort in cell cj. "          \
            "cj->nodeID=%d "                                                \
            "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->" #TYPE                \
            ".dx_max_sort=%e "                                              \
            "cj->" #TYPE                                                    \
            ".dx_max_sort_old=%e, cellID=%i super->cellID=%i"               \
            "cj->depth=%d cj->maxdepth=%d",                                 \
            cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->TYPE.dx_max_sort, \
            cj->TYPE.dx_max_sort_old, cj->cellID, cj->hydro.super->cellID,  \
            cj->depth, cj->maxdepth);                                       \
    }                                                                       \
  })

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

  const struct engine *restrict e = r->e;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  const int ci_active = 
      cell_is_active_stars(ci, e) || cell_is_active_hydro(ci, e);
  const int cj_active = 
      cell_is_active_stars(cj, e) || cell_is_active_hydro(cj, e);
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
      ((!cell_are_spart_drifted(ci, e) && !cell_are_spart_drifted(ci, e)) || 
      !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_ci && (!(ci->stars.sorted & (1 << sid)) ||
                ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

  if (do_ci && (!(cj->stars.sorted & (1 << sid)) ||
                cj->stars.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

  if (do_cj &&
      ((!cell_are_spart_drifted(ci, e) && !cell_are_spart_drifted(ci, e)) || 
      !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_cj && (!(ci->stars.sorted & (1 << sid)) ||
                ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

  if (do_cj && (!(cj->stars.sorted & (1 << sid)) ||
                cj->stars.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
  if (do_ci) {
    // MATTHIEU: This test is faulty. To be fixed...
    // RUNNER_CHECK_SORT(hydro, part, cj, ci, sid);
    RUNNER_CHECK_SORT(stars, spart, cj, ci, sid);
    RUNNER_CHECK_SORT(stars, spart, ci, cj, sid);
  }

  if (do_cj) {
    // MATTHIEU: This test is faulty. To be fixed...
    // RUNNER_CHECK_SORT(hydro, part, ci, cj, sid);
    RUNNER_CHECK_SORT(stars, spart, ci, cj, sid);
    RUNNER_CHECK_SORT(stars, spart, cj, ci, sid);
  }
#endif /* SWIFT_DEBUG_CHECKS */

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_SVD
  DOPAIR1_SVD_NAIVE(r, ci, cj, 1);
#else
  DO_SYM_PAIR1_SVD(r, ci, cj, sid, shift);
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

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  const int should_do_ci = ci->stars.count != 0 && cj->stars.count != 0 &&
      (cell_is_active_stars(ci, e) || cell_is_active_hydro(ci, e));
  const int should_do_cj = cj->stars.count != 0 && ci->stars.count != 0 &&
      (cell_is_active_stars(cj, e) || cell_is_active_hydro(cj, e));
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
    const int do_ci = ci->stars.count != 0 && cj->stars.count != 0 &&
        (cell_is_active_stars(ci, e) || cell_is_active_hydro(ci, e)) && 
        do_ci_stars;
    const int do_cj = cj->stars.count != 0 && ci->stars.count != 0 &&
        (cell_is_active_stars(cj, e) || cell_is_active_hydro(cj, e)) && 
        do_cj_stars;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e) && !cell_are_part_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (!cell_are_spart_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");

      /* Do any of the cells need to be sorted first? */
      if (!(ci->stars.sorted & (1 << sid)) ||
          ci->stars.dx_max_sort_old > ci->dmin * space_maxreldx) {
        error("Interacting unsorted cell (sparts).");
      }

      if (!(cj->stars.sorted & (1 << sid)) ||
          cj->stars.dx_max_sort_old > cj->dmin * space_maxreldx)
        error("Interacting unsorted cell (sparts). %i", cj->nodeID);
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (!cell_are_spart_drifted(cj, e) && !cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");

      /* Do any of the cells need to be sorted first? */
      if (!(ci->stars.sorted & (1 << sid)) ||
          ci->stars.dx_max_sort_old > ci->dmin * space_maxreldx) {
        error("Interacting unsorted cell (sparts).");
      }

      if (!(cj->stars.sorted & (1 << sid)) ||
          cj->stars.dx_max_sort_old > cj->dmin * space_maxreldx) {
        error("Interacting unsorted cell (sparts).");
      }
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_SVD(r, ci, cj);
  }

  TIMER_TOC(TIMER_DOSUB_PAIR_SVD);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_SVD(struct runner *r, struct cell *ci, int gettimer) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (ci->stars.count == 0 ||
      (!cell_is_active_stars(ci, r->e) && !cell_is_active_hydro(ci, r->e)))
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
    if (!cell_are_spart_drifted(ci, r->e) && !cell_are_spart_drifted(ci, r->e))
        error("Interacting undrifted cell.");

    DOSELF1_BRANCH_SVD(r, ci);
  }

  TIMER_TOC(TIMER_DOSUB_SELF_SVD);
}
