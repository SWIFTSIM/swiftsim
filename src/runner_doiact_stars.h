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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_doself_FUNCTION and runner_dosub_FUNCTION
   calling the pairwise interaction function runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOSELF1_STARS(f) PASTE(runner_doself_stars, f)
#define DOSELF1_STARS _DOSELF1_STARS(FUNCTION)

#define _DO_SYM_PAIR1_STARS(f) PASTE(runner_do_sym_pair_stars, f)
#define DO_SYM_PAIR1_STARS _DO_SYM_PAIR1_STARS(FUNCTION)

#define _DO_NONSYM_PAIR1_STARS_NAIVE(f) \
  PASTE(runner_do_nonsym_pair_stars_naive, f)
#define DO_NONSYM_PAIR1_STARS_NAIVE _DO_NONSYM_PAIR1_STARS_NAIVE(FUNCTION)

#define _DOPAIR1_STARS_NAIVE(f) PASTE(runner_dopair_stars_naive, f)
#define DOPAIR1_STARS_NAIVE _DOPAIR1_STARS_NAIVE(FUNCTION)

#define _DOPAIR1_SUBSET_STARS(f) PASTE(runner_dopair_subset_stars, f)
#define DOPAIR1_SUBSET_STARS _DOPAIR1_SUBSET_STARS(FUNCTION)

#define _DOPAIR1_SUBSET_STARS_NAIVE(f) \
  PASTE(runner_dopair_subset_stars_naive, f)
#define DOPAIR1_SUBSET_STARS_NAIVE _DOPAIR1_SUBSET_STARS_NAIVE(FUNCTION)

#define _DOSELF1_SUBSET_STARS(f) PASTE(runner_doself_subset_stars, f)
#define DOSELF1_SUBSET_STARS _DOSELF1_SUBSET_STARS(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_STARS(f) \
  PASTE(runner_doself_subset_branch_stars, f)
#define DOSELF1_SUBSET_BRANCH_STARS _DOSELF1_SUBSET_BRANCH_STARS(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_STARS(f) \
  PASTE(runner_dopair_subset_branch_stars, f)
#define DOPAIR1_SUBSET_BRANCH_STARS _DOPAIR1_SUBSET_BRANCH_STARS(FUNCTION)

#define _DOSUB_SUBSET_STARS(f) PASTE(runner_dosub_subset_stars, f)
#define DOSUB_SUBSET_STARS _DOSUB_SUBSET_STARS(FUNCTION)

#define _DOSELF1_BRANCH_STARS(f) PASTE(runner_doself_branch_stars, f)
#define DOSELF1_BRANCH_STARS _DOSELF1_BRANCH_STARS(FUNCTION)

#define _DOPAIR1_BRANCH_STARS(f) PASTE(runner_dopair_branch_stars, f)
#define DOPAIR1_BRANCH_STARS _DOPAIR1_BRANCH_STARS(FUNCTION)

#define _DOSUB_PAIR1_STARS(f) PASTE(runner_dosub_pair_stars, f)
#define DOSUB_PAIR1_STARS _DOSUB_PAIR1_STARS(FUNCTION)

#define _DOSUB_SELF1_STARS(f) PASTE(runner_dosub_self_stars, f)
#define DOSUB_SELF1_STARS _DOSUB_SELF1_STARS(FUNCTION)

#define _TIMER_DOSELF_STARS(f) PASTE(timer_doself_stars, f)
#define TIMER_DOSELF_STARS _TIMER_DOSELF_STARS(FUNCTION)

#define _TIMER_DOPAIR_STARS(f) PASTE(timer_dopair_stars, f)
#define TIMER_DOPAIR_STARS _TIMER_DOPAIR_STARS(FUNCTION)

#define _TIMER_DOSUB_SELF_STARS(f) PASTE(timer_dosub_self_stars, f)
#define TIMER_DOSUB_SELF_STARS _TIMER_DOSUB_SELF_STARS(FUNCTION)

#define _TIMER_DOSUB_PAIR_STARS(f) PASTE(timer_dosub_pair_stars, f)
#define TIMER_DOSUB_PAIR_STARS _TIMER_DOSUB_PAIR_STARS(FUNCTION)

#define _IACT_STARS(f) PASTE(runner_iact_nonsym_stars, f)
#define IACT_STARS _IACT_STARS(FUNCTION)

<<<<<<< HEAD
/**
 * @brief Calculate the number density of #part around the #spart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_STARS(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->stars.count == 0) return;
  if (!cell_is_active_stars(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  const int count = c->hydro.count;
  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts[sid];

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    /* Skip inactive particles */
    if (!feedback_is_active(si, e->time, cosmo, with_cosmology)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];
      struct xpart *restrict xpj = &xparts[pjd];
      const float hj = pj->h;

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

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

      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, si, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                            ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                          ti_current, e->time, e->step);
#endif
      }
    } /* loop over the parts in ci. */
  }   /* loop over the sparts in ci. */

  TIMER_TOC(TIMER_DOSELF_STARS);
}

/**
 * @brief Calculate the number density of cj #part around the ci #spart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_STARS_NAIVE(struct runner *r, struct cell *restrict ci,
                                 struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#else
  if (cj->nodeID != engine_rank) error("Should be run on a different node");
#endif
#endif

  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (cj->hydro.count == 0 || ci->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;

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

    /* Skip inactive particles */
    if (!feedback_is_active(si, e->time, cosmo, with_cosmology)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = &xparts_j[pjd];
      const float hj = pj->h;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

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

      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, si, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                            ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                          ti_current, e->time, e->step);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
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
void DO_SYM_PAIR1_STARS(struct runner *r, struct cell *ci, struct cell *cj,
                        const int sid, const double *shift) {

  TIMER_TIC;

  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_stars = (ci->nodeID == e->nodeID) && (ci->stars.count != 0) &&
                          (cj->hydro.count != 0) && cell_is_active_stars(ci, e);
  const int do_cj_stars = (cj->nodeID == e->nodeID) && (cj->stars.count != 0) &&
                          (ci->hydro.count != 0) && cell_is_active_stars(cj, e);
#else
  /* here we are updating the hydro -> switch ci, cj for local */
  const int do_ci_stars = (cj->nodeID == e->nodeID) && (ci->stars.count != 0) &&
                          (cj->hydro.count != 0) && cell_is_active_stars(ci, e);
  const int do_cj_stars = (ci->nodeID == e->nodeID) && (cj->stars.count != 0) &&
                          (ci->hydro.count != 0) && cell_is_active_stars(cj, e);
#endif

  if (do_ci_stars) {

    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_j = cj->hydro.sort[sid];
    const struct sort_entry *restrict sort_i = ci->stars.sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->stars.dx_max_part, cj->hydro.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->stars.dx_max_part, cj->hydro.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->stars.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const double hi_max = ci->stars.h_max * kernel_gamma - rshift;
    const int count_i = ci->stars.count;
    const int count_j = cj->hydro.count;
    struct spart *restrict sparts_i = ci->stars.parts;
    struct part *restrict parts_j = cj->hydro.parts;
    struct xpart *restrict xparts_j = cj->hydro.xparts;
    const double dj_min = sort_j[0].d;
    const float dx_max_rshift =
        (ci->stars.dx_max_sort + cj->hydro.dx_max_sort) - rshift;
    const float dx_max = (ci->stars.dx_max_sort + cj->hydro.dx_max_sort);

    /* Loop over the sparts in ci. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct spart *restrict spi = &sparts_i[sort_i[pid].i];
      const float hi = spi->h;

      /* Skip inactive particles */
      if (!spart_is_active(spi, e)) continue;

      /* Skip inactive particles */
      if (!feedback_is_active(spi, e->time, cosmo, with_cosmology)) continue;

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
        struct part *pj = &parts_j[sort_j[pjd].i];
        struct xpart *xpj = &xparts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

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
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, xpj,
                                              cosmo, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj, cosmo,
                                            ti_current, e->time, e->step);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* do_ci_stars */

  if (do_cj_stars) {
    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = ci->hydro.sort[sid];
    const struct sort_entry *restrict sort_j = cj->stars.sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->hydro.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->hydro.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->hydro.dx_max_part, cj->stars.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const double hj_max = cj->hydro.h_max * kernel_gamma;
    const int count_i = ci->hydro.count;
    const int count_j = cj->stars.count;
    struct part *restrict parts_i = ci->hydro.parts;
    struct xpart *restrict xparts_i = ci->hydro.xparts;
    struct spart *restrict sparts_j = cj->stars.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const float dx_max_rshift =
        (ci->hydro.dx_max_sort + cj->stars.dx_max_sort) + rshift;
    const float dx_max = (ci->hydro.dx_max_sort + cj->stars.dx_max_sort);

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct spart *spj = &sparts_j[sort_j[pjd].i];
      const float hj = spj->h;

      /* Skip inactive particles */
      if (!spart_is_active(spj, e)) continue;

      /* Skip inactive particles */
      if (!feedback_is_active(spj, e->time, cosmo, with_cosmology)) continue;

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
        struct part *pi = &parts_i[sort_i[pid].i];
        struct xpart *xpi = &xparts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

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
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (spj->ti_drift != e->ti_current)
          error("Particle spj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {

          IACT_STARS(r2, dx, hj, hi, spj, pi, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hj, hi, spj, pi, xpi,
                                              cosmo, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          runner_iact_nonsym_feedback_apply(r2, dx, hj, hi, spj, pi, xpi, cosmo,
                                            ti_current, e->time, e->step);
#endif
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR_STARS);
}

void DOPAIR1_STARS_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct cell *restrict cj, int timer) {

  TIMER_TIC;

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_stars = ci->nodeID == r->e->nodeID;
  const int do_cj_stars = cj->nodeID == r->e->nodeID;
#else
  /* here we are updating the hydro -> switch ci, cj */
  const int do_ci_stars = cj->nodeID == r->e->nodeID;
  const int do_cj_stars = ci->nodeID == r->e->nodeID;
#endif
  if (do_ci_stars && ci->stars.count != 0 && cj->hydro.count != 0)
    DO_NONSYM_PAIR1_STARS_NAIVE(r, ci, cj);
  if (do_cj_stars && cj->stars.count != 0 && ci->hydro.count != 0)
    DO_NONSYM_PAIR1_STARS_NAIVE(r, cj, ci);

  TIMER_TOC(TIMER_DOPAIR_STARS);
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
 * @param sid The direction of the pair.
 * @param flipped Flag to check whether the cells have been flipped or not.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_STARS(struct runner *r, struct cell *restrict ci,
                          struct spart *restrict sparts_i, int *restrict ind,
                          int scount, struct cell *restrict cj, const int sid,
                          const int flipped, const double *shift) {

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_j = cj->hydro.sort[sid];
  const float dxj = cj->hydro.dx_max_sort;

  /* Sparts are on the left? */
  if (!flipped) {

    /* Loop over the sparts_i. */
    for (int pid = 0; pid < scount; pid++) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict spi = &sparts_i[ind[pid]];
      const double pix = spi->x[0] - (shift[0]);
      const double piy = spi->x[1] - (shift[1]);
      const double piz = spi->x[2] - (shift[2]);
      const float hi = spi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = hi * kernel_gamma + dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];
        struct xpart *restrict xpj = &xparts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];
        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, xpj,
                                              cosmo, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj, cosmo,
                                            ti_current, e->time, e->step);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the sparts in ci. */
  }

  /* Sparts are on the right. */
  else {

    /* Loop over the sparts_i. */
    for (int pid = 0; pid < scount; pid++) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict spi = &sparts_i[ind[pid]];
      const double pix = spi->x[0] - (shift[0]);
      const double piy = spi->x[1] - (shift[1]);
      const double piz = spi->x[2] - (shift[2]);
      const float hi = spi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = -hi * kernel_gamma - dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];
        struct xpart *restrict xpj = &xparts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];
        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, xpj,
                                              cosmo, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj, cosmo,
                                            ti_current, e->time, e->step);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the sparts in ci. */
  }
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
void DOPAIR1_SUBSET_STARS_NAIVE(struct runner *r, struct cell *restrict ci,
                                struct spart *restrict sparts_i,
                                int *restrict ind, int scount,
                                struct cell *restrict cj, const double *shift) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the parts_i. */
  for (int pid = 0; pid < scount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *restrict spi = &sparts_i[ind[pid]];

    const double pix = spi->x[0] - (shift[0]);
    const double piy = spi->x[1] - (shift[1]);
    const double piz = spi->x[2] - (shift[2]);
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
      struct xpart *restrict xpj = &xparts_j[pjd];

      /* Skip inhibited particles */
      if (part_is_inhibited(pj, e)) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                     (float)(piz - pjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, xpj, cosmo,
                                            ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj, cosmo,
                                          ti_current, e->time, e->step);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
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
void DOSELF1_SUBSET_STARS(struct runner *r, struct cell *restrict ci,
                          struct spart *restrict sparts, int *restrict ind,
                          int scount) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;
  struct xpart *restrict xparts_j = ci->hydro.xparts;

  /* Early abort? */
  if (count_i == 0) return;

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
      struct xpart *restrict xpj = &xparts_j[pjd];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

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
      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, pj->h, spi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, pj->h, spi, pj, xpj,
                                            cosmo, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        runner_iact_nonsym_feedback_apply(r2, dx, hi, pj->h, spi, pj, xpj,
                                          cosmo, ti_current, e->time, e->step);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}
=======
void DOSELF1_BRANCH_STARS(struct runner *r, struct cell *c);
void DOPAIR1_BRANCH_STARS(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1_STARS(struct runner *r, struct cell *ci, int gettimer);
void DOSUB_PAIR1_STARS(struct runner *r, struct cell *ci, struct cell *cj,
                       int gettimer);
>>>>>>> upstream/master

void DOSELF1_SUBSET_BRANCH_STARS(struct runner *r, struct cell *restrict ci,
                                 struct spart *restrict sparts,
                                 int *restrict ind, int scount);

void DOPAIR1_SUBSET_BRANCH_STARS(struct runner *r, struct cell *restrict ci,
                                 struct spart *restrict sparts_i,
                                 int *restrict ind, int scount,
                                 struct cell *restrict cj);

void DOSUB_SUBSET_STARS(struct runner *r, struct cell *ci, struct spart *sparts,
                        int *ind, int scount, struct cell *cj, int gettimer);
