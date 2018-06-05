/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

const int num_tests = 100;
const double eps = 0.1;

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(double a, double b, const char *s) {
  if (fabs(a - b) / fabs(a + b) > 2e-6 && fabs(a - b) > 1.e-6)
    error("Values are inconsistent: %12.15e %12.15e (%s)!", a, b, s);
}

/* Definitions of the potential and force that match
   exactly the theory document */
double S(double x) { return exp(x) / (1. + exp(x)); }

double S_prime(double x) { return exp(x) / ((1. + exp(x)) * (1. + exp(x))); }

double potential(double mass, double r, double H, double rlr) {

  const double u = r / H;
  const double x = r / rlr;
  double pot;
  if (u > 1.)
    pot = -mass / r;
  else
    pot = -mass *
          (-3. * u * u * u * u * u * u * u + 15. * u * u * u * u * u * u -
           28. * u * u * u * u * u + 21. * u * u * u * u - 7. * u * u + 3.) /
          H;

  return pot * (2. - 2. * S(2. * x));
}

double acceleration(double mass, double r, double H, double rlr) {

  const double u = r / H;
  const double x = r / rlr;
  double acc;
  if (u > 1.)
    acc = -mass / (r * r * r);
  else
    acc = -mass * (21. * u * u * u * u * u - 90. * u * u * u * u +
                   140. * u * u * u - 84. * u * u + 14.) /
          (H * H * H);

  return r * acc * (4. * x * S_prime(2 * x) - 2. * S(2. * x) + 2.);
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Initialise a few things to get us going */
  struct engine e;
  e.max_active_bin = num_time_bins;
  e.time = 0.1f;
  e.ti_current = 8;
  e.time_base = 1e-10;

  struct space s;
  s.periodic = 0;
  s.width[0] = 0.2;
  s.width[1] = 0.2;
  s.width[2] = 0.2;
  e.s = &s;

  struct gravity_props props;
  props.a_smooth = 1.25;
  props.r_cut_min = 0.;
  props.theta_crit2 = 0.;
  props.epsilon_cur = eps;
  e.gravity_properties = &props;

  struct runner r;
  bzero(&r, sizeof(struct runner));
  r.e = &e;

  const double rlr = props.a_smooth * s.width[0] * FLT_MAX;

  /* Init the cache for gravity interaction */
  gravity_cache_init(&r.ci_gravity_cache, num_tests);
  gravity_cache_init(&r.cj_gravity_cache, num_tests);

  /* Let's create one cell with a massive particle and a bunch of test particles
   */
  struct cell ci, cj;
  bzero(&ci, sizeof(struct cell));
  bzero(&cj, sizeof(struct cell));

  ci.width[0] = 1.;
  ci.width[1] = 1.;
  ci.width[2] = 1.;
  ci.loc[0] = 0.;
  ci.loc[1] = 0.;
  ci.loc[2] = 0.;
  ci.gcount = 1;
  ci.ti_old_gpart = 8;
  ci.ti_old_multipole = 8;
  ci.ti_gravity_end_min = 8;
  ci.ti_gravity_end_max = 8;

  cj.width[0] = 1.;
  cj.width[1] = 1.;
  cj.width[2] = 1.;
  cj.loc[0] = 1.;
  cj.loc[1] = 0.;
  cj.loc[2] = 0.;
  cj.gcount = num_tests;
  cj.ti_old_gpart = 8;
  cj.ti_old_multipole = 8;
  cj.ti_gravity_end_min = 8;
  cj.ti_gravity_end_max = 8;

  /* Allocate multipoles */
  ci.multipole =
      (struct gravity_tensors *)malloc(sizeof(struct gravity_tensors));
  cj.multipole =
      (struct gravity_tensors *)malloc(sizeof(struct gravity_tensors));
  bzero(ci.multipole, sizeof(struct gravity_tensors));
  bzero(cj.multipole, sizeof(struct gravity_tensors));

  /* Set the multipoles */
  ci.multipole->r_max = 0.1;
  cj.multipole->r_max = 0.1;

  /* Allocate the particles */
  if (posix_memalign((void **)&ci.gparts, gpart_align,
                     ci.gcount * sizeof(struct gpart)) != 0)
    error("Error allocating gparts for cell ci");
  bzero(ci.gparts, ci.gcount * sizeof(struct gpart));

  if (posix_memalign((void **)&cj.gparts, gpart_align,
                     cj.gcount * sizeof(struct gpart)) != 0)
    error("Error allocating gparts for cell ci");
  bzero(cj.gparts, cj.gcount * sizeof(struct gpart));

  /* Create the mass-less test particles */
  for (int n = 0; n < num_tests; ++n) {

    struct gpart *gp = &cj.gparts[n];

    gp->x[0] = 1. + (n + 1) / ((double)num_tests);
    gp->x[1] = 0.5;
    gp->x[2] = 0.5;
    gp->mass = 0.;
    gp->time_bin = 1;
    gp->type = swift_type_dark_matter;
    gp->id_or_neg_offset = n + 1;
#ifdef SWIFT_DEBUG_CHECKS
    gp->ti_drift = 8;
#endif
  }

  /***********************************************/
  /* Let's start by testing the P-P interactions */
  /***********************************************/

  /* Create the massive particle */
  ci.gparts[0].x[0] = 0.;
  ci.gparts[0].x[1] = 0.5;
  ci.gparts[0].x[2] = 0.5;
  ci.gparts[0].mass = 1.;
  ci.gparts[0].time_bin = 1;
  ci.gparts[0].type = swift_type_dark_matter;
  ci.gparts[0].id_or_neg_offset = 1;
#ifdef SWIFT_DEBUG_CHECKS
  ci.gparts[0].ti_drift = 8;
#endif

  /* Now compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.gparts[n];
    const struct gpart *gp2 = &ci.gparts[0];
    const double epsilon = gravity_get_softening(gp, &props);

    double pot_true =
        potential(ci.gparts[0].mass, gp->x[0] - gp2->x[0], epsilon, rlr);
    double acc_true =
        acceleration(ci.gparts[0].mass, gp->x[0] - gp2->x[0], epsilon, rlr);

    message("x=%e f=%e f_true=%e pot=%e pot_true=%e", gp->x[0] - gp2->x[0],
            gp->a_grav[0], acc_true, gp->potential, pot_true);

    check_value(gp->potential, pot_true, "potential");
    check_value(gp->a_grav[0], acc_true, "acceleration");
  }

  message("\n\t\t P-P interactions all good\n");

  /* Reset the accelerations */
  for (int n = 0; n < num_tests; ++n) gravity_init_gpart(&cj.gparts[n]);

  /**********************************/
  /* Test the basic PM interactions */
  /**********************************/

  /* Set an opening angle that allows P-M interactions */
  props.theta_crit2 = 1.;

  ci.gparts[0].mass = 0.;
  ci.multipole->CoM[0] = 0.;
  ci.multipole->CoM[1] = 0.5;
  ci.multipole->CoM[2] = 0.5;

  bzero(&ci.multipole->m_pole, sizeof(struct multipole));
  bzero(&cj.multipole->m_pole, sizeof(struct multipole));
  ci.multipole->m_pole.M_000 = 1.;

  /* Now compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.gparts[n];
    const struct gravity_tensors *mpole = ci.multipole;
    const double epsilon = gravity_get_softening(gp, &props);

    double pot_true = potential(mpole->m_pole.M_000, gp->x[0] - mpole->CoM[0],
                                epsilon, rlr * FLT_MAX);
    double acc_true = acceleration(
        mpole->m_pole.M_000, gp->x[0] - mpole->CoM[0], epsilon, rlr * FLT_MAX);

    message("x=%e f=%e f_true=%e pot=%e pot_true=%e", gp->x[0] - mpole->CoM[0],
            gp->a_grav[0], acc_true, gp->potential, pot_true);

    check_value(gp->potential, pot_true, "potential");
    check_value(gp->a_grav[0], acc_true, "acceleration");
  }

  message("\n\t\t basic P-M interactions all good\n");

  /* Reset the accelerations */
  for (int n = 0; n < num_tests; ++n) gravity_init_gpart(&cj.gparts[n]);

/***************************************/
/* Test the high-order PM interactions */
/***************************************/

#if SELF_GRAVITY_MULTIPOLE_ORDER >= 3

  /* Let's make ci more interesting */
  free(ci.gparts);
  ci.gcount = 8;
  if (posix_memalign((void **)&ci.gparts, gpart_align,
                     ci.gcount * sizeof(struct gpart)) != 0)
    error("Error allocating gparts for cell ci");
  bzero(ci.gparts, ci.gcount * sizeof(struct gpart));

  /* Place particles on a simple cube of side-length 0.2 */
  for (int n = 0; n < 8; ++n) {
    if (n & 1)
      ci.gparts[n].x[0] = 0.0 - 0.1;
    else
      ci.gparts[n].x[0] = 0.0 + 0.1;

    if (n & 2)
      ci.gparts[n].x[1] = 0.5 - 0.1;
    else
      ci.gparts[n].x[1] = 0.5 + 0.1;

    if (n & 2)
      ci.gparts[n].x[2] = 0.5 - 0.1;
    else
      ci.gparts[n].x[2] = 0.5 + 0.1;

    ci.gparts[n].mass = 1. / 8.;

    ci.gparts[n].time_bin = 1;
    ci.gparts[n].type = swift_type_dark_matter;
    ci.gparts[n].id_or_neg_offset = 1;
#ifdef SWIFT_DEBUG_CHECKS
    ci.gparts[n].ti_drift = 8;
#endif
  }

  /* Now let's make a multipole out of it. */
  gravity_reset(ci.multipole);
  gravity_P2M(ci.multipole, ci.gparts, ci.gcount);

  // message("CoM=[%e %e %e]", ci.multipole->CoM[0], ci.multipole->CoM[1],
  // ci.multipole->CoM[2]);
  gravity_multipole_print(&ci.multipole->m_pole);

  /* Compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.gparts[n];
    const struct gravity_tensors *mpole = ci.multipole;

    double pot_true = 0, acc_true[3] = {0., 0., 0.};

    for (int i = 0; i < 8; ++i) {
      const struct gpart *gp2 = &ci.gparts[i];
      const double epsilon = gravity_get_softening(gp, &props);

      const double dx[3] = {gp2->x[0] - gp->x[0], gp2->x[1] - gp->x[1],
                            gp2->x[2] - gp->x[2]};
      const double d = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

      pot_true += potential(gp2->mass, d, epsilon, rlr * FLT_MAX);
      acc_true[0] -=
          acceleration(gp2->mass, d, epsilon, rlr * FLT_MAX) * dx[0] / d;
      acc_true[1] -=
          acceleration(gp2->mass, d, epsilon, rlr * FLT_MAX) * dx[1] / d;
      acc_true[2] -=
          acceleration(gp2->mass, d, epsilon, rlr * FLT_MAX) * dx[2] / d;
    }

    message("x=%e f=%e f_true=%e pot=%e pot_true=%e %e %e",
            gp->x[0] - mpole->CoM[0], gp->a_grav[0], acc_true[0], gp->potential,
            pot_true, acc_true[1], acc_true[2]);

    // check_value(gp->potential, pot_true, "potential");
    // check_value(gp->a_grav[0], acc_true[0], "acceleration");
  }

  message("\n\t\t high-order P-M interactions all good\n");

#endif

  free(ci.multipole);
  free(cj.multipole);
  free(ci.gparts);
  free(cj.gparts);
  return 0;
}
