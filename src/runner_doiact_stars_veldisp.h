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

#define _DOSELF1_SVD(f) PASTE(runner_doself_stars, f)
#define DOSELF1_SVD _DOSELF1_SVD(FUNCTION)

#define _DO_SYM_PAIR1_SVD(f) PASTE(runner_do_sym_pair_stars, f)
#define DO_SYM_PAIR1_SVD _DO_SYM_PAIR1_SVD(FUNCTION)

#define _DO_NONSYM_PAIR1_SVD_NAIVE(f) \
  PASTE(runner_do_nonsym_pair_stars_naive, f)
#define DO_NONSYM_PAIR1_SVD_NAIVE _DO_NONSYM_PAIR1_SVD_NAIVE(FUNCTION)

#define _DOPAIR1_SVD_NAIVE(f) PASTE(runner_dopair_stars_naive, f)
#define DOPAIR1_SVD_NAIVE _DOPAIR1_SVD_NAIVE(FUNCTION)

#define _DOSELF1_BRANCH_SVD(f) PASTE(runner_doself_branch_stars, f)
#define DOSELF1_BRANCH_SVD _DOSELF1_BRANCH_SVD(FUNCTION)

#define _DOPAIR1_BRANCH_SVD(f) PASTE(runner_dopair_branch_stars, f)
#define DOPAIR1_BRANCH_SVD _DOPAIR1_BRANCH_SVD(FUNCTION)

#define _DOSUB_PAIR1_SVD(f) PASTE(runner_dosub_pair_stars, f)
#define DOSUB_PAIR1_SVD _DOSUB_PAIR1_SVD(FUNCTION)

#define _DOSUB_SELF1_SVD(f) PASTE(runner_dosub_self_stars, f)
#define DOSUB_SELF1_SVD _DOSUB_SELF1_SVD(FUNCTION)

#define _TIMER_DOSELF_SVD(f) PASTE(timer_doself_stars, f)
#define TIMER_DOSELF_SVD _TIMER_DOSELF_SVD(FUNCTION)

#define _TIMER_DOPAIR_SVD(f) PASTE(timer_dopair_stars, f)
#define TIMER_DOPAIR_SVD _TIMER_DOPAIR_SVD(FUNCTION)

#define _TIMER_DOSUB_SELF_SVD(f) PASTE(timer_dosub_self_stars, f)
#define TIMER_DOSUB_SELF_SVD _TIMER_DOSUB_SELF_SVD(FUNCTION)

#define _TIMER_DOSUB_PAIR_SVD(f) PASTE(timer_dosub_pair_stars, f)
#define TIMER_DOSUB_PAIR_SVD _TIMER_DOSUB_PAIR_SVD(FUNCTION)

#define _IACT_SVD(f) PASTE(runner_iact_nonsym_stars, f)
#define IACT_SVD _IACT_SVD(FUNCTION)

void DOSELF1_BRANCH_SVD(struct runner *r, struct cell *c);
void DOPAIR1_BRANCH_SVD(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1_SVD(struct runner *r, struct cell *ci, int gettimer);
void DOSUB_PAIR1_SVD(struct runner *r, struct cell *ci, struct cell *cj,
                       int gettimer);
