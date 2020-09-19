/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "black_holes.h"
#include "cell.h"
#include "engine.h"
#include "feedback.h"
#include "kick.h"
#include "multipole.h"
#include "timers.h"
#include "timestep.h"
#include "timestep_limiter.h"
#include "timestep_sync.h"
#include "tracers.h"

/**
 * @brief Initialize the multipoles before the gravity calculation.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_init_grav(struct runner *r, struct cell *c, const int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (!(e->policy & engine_policy_self_gravity))
    error("Grav-init task called outside of self-gravity calculation");
#endif

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Reset the gravity acceleration tensors */
  gravity_field_tensors_init(&c->grav.multipole->pot, e->ti_current);

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        runner_do_init_grav(r, c->progeny[k], /*timer=*/0);
    }
  }

  if (timer) TIMER_TOC(timer_init_grav);
}

/**
 * @brief Perform the first half-kick on all the active particles in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_kick1(struct runner *r, struct cell *c, const int timer) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const int periodic = e->s->periodic;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  struct sink *restrict sinks = c->sinks.parts;
  struct bpart *restrict bparts = c->black_holes.parts;
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int sink_count = c->sinks.count;
  const int bcount = c->black_holes.count;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(c, e) && !cell_is_starting_gravity(c, e) &&
      !cell_is_starting_stars(c, e) && !cell_is_starting_sinks(c, e) &&
      !cell_is_starting_black_holes(c, e))
    return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_kick1(r, c->progeny[k], /*timer=*/0);
  } else {

    integertime_t ti_begin_mesh = -1;
    integertime_t ti_end_mesh = -1;
    double dt_kick_mesh_grav = 0.;

    /* Are we at a step where we do mesh-gravity time-integration? */
    if (periodic && e->mesh->ti_beg_mesh_next == e->ti_current) {

      const integertime_t ti_step =
          e->mesh->ti_end_mesh_next - e->mesh->ti_beg_mesh_next;
      ti_begin_mesh = e->mesh->ti_beg_mesh_next;
      ti_end_mesh = e->mesh->ti_beg_mesh_next + ti_step / 2;

      /* Time interval for this mesh gravity half-kick */
      dt_kick_mesh_grav = kick_get_grav_kick_dt(
          ti_begin_mesh, ti_end_mesh, time_base, with_cosmology, cosmo);
    }

    /* Loop over the parts in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be kicked */
      if (part_is_starting(p, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        if (p->limiter_data.wakeup != time_bin_not_awake)
          error("Woken-up particle that has not been processed in kick1");
#endif

        const integertime_t ti_step = get_integer_timestep(p->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, p->time_bin);
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_begin != ti_current)
          error(
              "Particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d wakeup=%d ti_current=%lld",
              ti_end, ti_begin, ti_step, p->time_bin, p->limiter_data.wakeup,
              ti_current);
#endif

        /* Time intervals for this half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);
        const double dt_kick_hydro = kick_get_hydro_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);
        const double dt_kick_therm = kick_get_therm_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);
        const double dt_kick_corr = kick_get_corr_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Do the kick */
        kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm,
                  dt_kick_corr, cosmo, hydro_props, entropy_floor, ti_begin,
                  ti_end);

        /* Update the accelerations to be used in the drift for hydro */
        if (p->gpart != NULL) {

          xp->a_grav[0] = p->gpart->a_grav[0];
          xp->a_grav[1] = p->gpart->a_grav[1];
          xp->a_grav[2] = p->gpart->a_grav[2];
        }
      }
    }

    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

#ifdef SWIFT_DEBUG_CHECKS
      if (ti_end_mesh != -1 && !gpart_is_starting(gp, e)) {
        error(
            "Particle on a time-step longer than the mesh synchronization "
            "step!");
      }
#endif

      /* If the g-particle has no counterpart and needs to be kicked */
      if ((gp->type == swift_type_dark_matter ||
           gp->type == swift_type_dark_matter_background) &&
          gpart_is_starting(gp, e)) {

        const integertime_t ti_step = get_integer_timestep(gp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, gp->time_bin);
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end_check =
            get_integer_time_end(ti_current + 1, gp->time_bin);

        if (ti_begin != ti_current)
          error(
              "G-particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end_check, ti_begin, ti_step, gp->time_bin, ti_current);
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Do the kick */
        kick_gpart(gp, dt_kick_grav, ti_begin, ti_end, dt_kick_mesh_grav,
                   ti_begin_mesh, ti_end_mesh);
      }
    }

    /* Loop over the stars particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the s-part. */
      struct spart *restrict sp = &sparts[k];

      /* If particle needs to be kicked */
      if (spart_is_starting(sp, e)) {

        const integertime_t ti_step = get_integer_timestep(sp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, sp->time_bin);
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end_check =
            get_integer_time_end(ti_current + 1, sp->time_bin);

        if (ti_begin != ti_current)
          error(
              "S-particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end_check, ti_begin, ti_step, sp->time_bin, ti_current);
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Do the kick */
        kick_spart(sp, dt_kick_grav, ti_begin, ti_end);
      }
    }

    /* Loop over the sink particles in this cell. */
    for (int k = 0; k < sink_count; k++) {

      /* Get a handle on the s-part. */
      struct sink *restrict sink = &sinks[k];

      /* If particle needs to be kicked */
      if (sink_is_starting(sink, e)) {

        const integertime_t ti_step = get_integer_timestep(sink->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, sink->time_bin);
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end_check =
            get_integer_time_end(ti_current + 1, sink->time_bin);

        if (ti_begin != ti_current)
          error(
              "Sink-particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end_check, ti_begin, ti_step, sink->time_bin, ti_current);
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Do the kick */
        kick_sink(sink, dt_kick_grav, ti_begin, ti_end);
      }
    }

    /* Loop over the black hole particles in this cell. */
    for (int k = 0; k < bcount; k++) {

      /* Get a handle on the s-part. */
      struct bpart *restrict bp = &bparts[k];

      /* If particle needs to be kicked */
      if (bpart_is_starting(bp, e)) {

        const integertime_t ti_step = get_integer_timestep(bp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, bp->time_bin);
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end_check =
            get_integer_time_end(ti_current + 1, bp->time_bin);

        if (ti_begin != ti_current)
          error(
              "B-particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end_check, ti_begin, ti_step, bp->time_bin, ti_current);
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Do the kick */
        kick_bpart(bp, dt_kick_grav, ti_begin, ti_end);
      }
    }
  }

  if (timer) TIMER_TOC(timer_kick1);
}

/**
 * @brief Perform the second half-kick on all the active particles in a cell.
 *
 * Also prepares particles to be drifted.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_kick2(struct runner *r, struct cell *c, const int timer) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int periodic = e->s->periodic;
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int sink_count = c->sinks.count;
  const int bcount = c->black_holes.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  struct sink *restrict sinks = c->sinks.parts;
  struct bpart *restrict bparts = c->black_holes.parts;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e) &&
      !cell_is_active_stars(c, e) && !cell_is_active_sinks(c, e) &&
      !cell_is_active_black_holes(c, e))
    return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_kick2(r, c->progeny[k], /*timer=*/0);
  } else {

    integertime_t ti_begin_mesh = -1;
    integertime_t ti_end_mesh = -1;
    double dt_kick_mesh_grav = 0.;

    /* Are we at a step where we do mesh-gravity time-integration? */
    if (periodic && e->mesh->ti_end_mesh_last == e->ti_current) {

      const integertime_t ti_step =
          e->mesh->ti_end_mesh_next - e->mesh->ti_beg_mesh_next;
      ti_begin_mesh = e->mesh->ti_beg_mesh_last + ti_step / 2;
      ti_end_mesh = e->mesh->ti_end_mesh_last;

      /* Time interval for this mesh gravity half-kick */
      dt_kick_mesh_grav = kick_get_grav_kick_dt(
          ti_begin_mesh, ti_end_mesh, time_base, with_cosmology, cosmo);
    }

    /* Loop over the particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be kicked */
      if (part_is_active(p, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        if (p->limiter_data.wakeup != time_bin_not_awake)
          error("Woken-up particle that has not been processed in kick1");
#endif
        /* Time-step length on the integer timeline */
        const integertime_t ti_step = get_integer_timestep(p->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, p->time_bin) + ti_step / 2;
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_end != ti_current)
          error(
              "Particle in wrong time-bin, ti_begin=%lld, ti_step=%lld "
              "time_bin=%d wakeup=%d ti_current=%lld",
              ti_begin, ti_step, p->time_bin, p->limiter_data.wakeup,
              ti_current);
#endif

        /* Time intervals for this half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);
        const double dt_kick_hydro = kick_get_hydro_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);
        const double dt_kick_therm = kick_get_therm_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);
        const double dt_kick_corr = kick_get_corr_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Finish the time-step with a second half-kick */
        kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm,
                  dt_kick_corr, cosmo, hydro_props, entropy_floor, ti_begin,
                  ti_end);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (p->ti_drift != p->ti_kick) error("Error integrating part in time.");
#endif

        /* Prepare the values to be drifted */
        hydro_reset_predicted_values(p, xp, cosmo);
      }
    }

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

#ifdef SWIFT_DEBUG_CHECKS
      if (ti_end_mesh != -1 && !gpart_is_active(gp, e)) {
        error(
            "Particle on a time-step longer than the mesh synchronization "
            "step!");
      }
#endif

      /* If the g-particle has no counterpart and needs to be kicked */
      if ((gp->type == swift_type_dark_matter ||
           gp->type == swift_type_dark_matter_background) &&
          gpart_is_active(gp, e)) {

        const integertime_t ti_step = get_integer_timestep(gp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, gp->time_bin) + ti_step / 2;
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_end != ti_current) error("Particle in wrong time-bin");
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Finish the time-step with a second half-kick */
        kick_gpart(gp, dt_kick_grav, ti_begin, ti_end, dt_kick_mesh_grav,
                   ti_begin_mesh, ti_end_mesh);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (gp->ti_drift != gp->ti_kick)
          error("Error integrating g-part in time.");
#endif

        /* Prepare the values to be drifted */
        gravity_reset_predicted_values(gp);
      }
    }

    /* Loop over the particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the part. */
      struct spart *restrict sp = &sparts[k];

      /* If particle needs to be kicked */
      if (spart_is_active(sp, e)) {

        const integertime_t ti_step = get_integer_timestep(sp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, sp->time_bin) + ti_step / 2;
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_end != ti_current) error("Particle in wrong time-bin");
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Finish the time-step with a second half-kick */
        kick_spart(sp, dt_kick_grav, ti_begin, ti_end);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (sp->ti_drift != sp->ti_kick)
          error("Error integrating s-part in time.");
#endif

        /* Prepare the values to be drifted */
        stars_reset_predicted_values(sp);
      }
    }

    /* Loop over the sink particles in this cell. */
    for (int k = 0; k < sink_count; k++) {

      /* Get a handle on the part. */
      struct sink *restrict sink = &sinks[k];

      /* If particle needs to be kicked */
      if (sink_is_active(sink, e)) {

        const integertime_t ti_step = get_integer_timestep(sink->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, sink->time_bin) + ti_step / 2;
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_end != ti_current) error("Particle in wrong time-bin");
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Finish the time-step with a second half-kick */
        kick_sink(sink, dt_kick_grav, ti_begin, ti_end);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (sink->ti_drift != sink->ti_kick)
          error("Error integrating sink-part in time.");
#endif

        /* Prepare the values to be drifted */
        sink_reset_predicted_values(sink);
      }
    }

    /* Loop over the particles in this cell. */
    for (int k = 0; k < bcount; k++) {

      /* Get a handle on the part. */
      struct bpart *restrict bp = &bparts[k];

      /* If particle needs to be kicked */
      if (bpart_is_active(bp, e)) {

        const integertime_t ti_step = get_integer_timestep(bp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, bp->time_bin) + ti_step / 2;
        const integertime_t ti_end = ti_begin + ti_step / 2;

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_end != ti_current) error("Particle in wrong time-bin");
#endif

        /* Time interval for this gravity half-kick */
        const double dt_kick_grav = kick_get_grav_kick_dt(
            ti_begin, ti_end, time_base, with_cosmology, cosmo);

        /* Finish the time-step with a second half-kick */
        kick_bpart(bp, dt_kick_grav, ti_begin, ti_end);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (bp->ti_drift != bp->ti_kick)
          error("Error integrating b-part in time.");
#endif

        /* Prepare the values to be drifted */
        black_holes_reset_predicted_values(bp);
      }
    }
  }
  if (timer) TIMER_TOC(timer_kick2);
}

/**
 * @brief Computes the next time-step of all active particles in this cell
 * and update the cell's statistics.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_timestep(struct runner *r, struct cell *c, const int timer) {

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int sink_count = c->sinks.count;
  const int bcount = c->black_holes.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  struct sink *restrict sinks = c->sinks.parts;
  struct bpart *restrict bparts = c->black_holes.parts;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e) &&
      !cell_is_active_stars(c, e) && !cell_is_active_sinks(c, e) &&
      !cell_is_active_black_holes(c, e)) {
    c->hydro.updated = 0;
    c->grav.updated = 0;
    c->stars.updated = 0;
    c->sinks.updated = 0;
    c->black_holes.updated = 0;
    return;
  }

  int updated = 0, g_updated = 0, s_updated = 0, sink_updated = 0,
      b_updated = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_end_max = 0,
                ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_end_max = 0,
                ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_end_max = 0, ti_black_holes_beg_max = 0;

  /* No children? */
  if (!c->split) {

    /* Loop over the particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs updating */
      if (part_is_active(p, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Current end of time-step */
        const integertime_t ti_end =
            get_integer_time_end(ti_current, p->time_bin);

        if (ti_end != ti_current)
          error("Computing time-step of rogue particle.");
#endif

        /* Get new time-step */
        const integertime_t ti_new_step = get_part_timestep(p, xp, e);

        /* Update particle */
        p->time_bin = get_time_bin(ti_new_step);
        if (p->gpart != NULL) p->gpart->time_bin = p->time_bin;

        /* Update the tracers properties */
        tracers_after_timestep(p, xp, e->internal_units, e->physical_constants,
                               with_cosmology, e->cosmology,
                               e->hydro_properties, e->cooling_func, e->time);

        /* Number of updated particles */
        updated++;
        if (p->gpart != NULL) g_updated++;

        /* What is the next sync-point ? */
        ti_hydro_end_min = min(ti_current + ti_new_step, ti_hydro_end_min);
        ti_hydro_end_max = max(ti_current + ti_new_step, ti_hydro_end_max);

        /* What is the next starting point for this cell ? */
        ti_hydro_beg_max = max(ti_current, ti_hydro_beg_max);

        if (p->gpart != NULL) {

          /* What is the next sync-point ? */
          ti_gravity_end_min =
              min(ti_current + ti_new_step, ti_gravity_end_min);
          ti_gravity_end_max =
              max(ti_current + ti_new_step, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);
        }
      }

      else { /* part is inactive */

        if (!part_is_inhibited(p, e)) {

          const integertime_t ti_end =
              get_integer_time_end(ti_current, p->time_bin);

          const integertime_t ti_beg =
              get_integer_time_begin(ti_current + 1, p->time_bin);

          /* What is the next sync-point ? */
          ti_hydro_end_min = min(ti_end, ti_hydro_end_min);
          ti_hydro_end_max = max(ti_end, ti_hydro_end_max);

          /* What is the next starting point for this cell ? */
          ti_hydro_beg_max = max(ti_beg, ti_hydro_beg_max);

          if (p->gpart != NULL) {

            /* What is the next sync-point ? */
            ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
            ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

            /* What is the next starting point for this cell ? */
            ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
          }
        }
      }
    }

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* If the g-particle has no counterpart */
      if (gp->type == swift_type_dark_matter ||
          gp->type == swift_type_dark_matter_background) {

        /* need to be updated ? */
        if (gpart_is_active(gp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
          /* Current end of time-step */
          const integertime_t ti_end =
              get_integer_time_end(ti_current, gp->time_bin);

          if (ti_end != ti_current)
            error("Computing time-step of rogue particle.");
#endif

          /* Get new time-step */
          const integertime_t ti_new_step = get_gpart_timestep(gp, e);

          /* Update particle */
          gp->time_bin = get_time_bin(ti_new_step);

          /* Number of updated g-particles */
          g_updated++;

          /* What is the next sync-point ? */
          ti_gravity_end_min =
              min(ti_current + ti_new_step, ti_gravity_end_min);
          ti_gravity_end_max =
              max(ti_current + ti_new_step, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

        } else { /* gpart is inactive */

          if (!gpart_is_inhibited(gp, e)) {

            const integertime_t ti_end =
                get_integer_time_end(ti_current, gp->time_bin);

            /* What is the next sync-point ? */
            ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
            ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

            const integertime_t ti_beg =
                get_integer_time_begin(ti_current + 1, gp->time_bin);

            /* What is the next starting point for this cell ? */
            ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
          }
        }
      }
    }

    /* Loop over the star particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the part. */
      struct spart *restrict sp = &sparts[k];

      /* need to be updated ? */
      if (spart_is_active(sp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Current end of time-step */
        const integertime_t ti_end =
            get_integer_time_end(ti_current, sp->time_bin);

        if (ti_end != ti_current)
          error("Computing time-step of rogue particle.");
#endif
        /* Get new time-step */
        const integertime_t ti_new_step = get_spart_timestep(sp, e);

        /* Update particle */
        sp->time_bin = get_time_bin(ti_new_step);
        sp->gpart->time_bin = get_time_bin(ti_new_step);

        /* Update feedback related counters */
        if (with_feedback)
          feedback_will_do_feedback(sp, e->feedback_props, with_cosmology,
                                    e->cosmology, e->time);

        /* Number of updated s-particles */
        s_updated++;
        g_updated++;

        ti_stars_end_min = min(ti_current + ti_new_step, ti_stars_end_min);
        ti_stars_end_max = max(ti_current + ti_new_step, ti_stars_end_max);
        ti_gravity_end_min = min(ti_current + ti_new_step, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_current + ti_new_step, ti_gravity_end_max);

        /* What is the next starting point for this cell ? */
        ti_stars_beg_max = max(ti_current, ti_stars_beg_max);
        ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

        /* star particle is inactive but not inhibited */
      } else {

        if (!spart_is_inhibited(sp, e)) {

          const integertime_t ti_end =
              get_integer_time_end(ti_current, sp->time_bin);

          const integertime_t ti_beg =
              get_integer_time_begin(ti_current + 1, sp->time_bin);

          ti_stars_end_min = min(ti_end, ti_stars_end_min);
          ti_stars_end_max = max(ti_end, ti_stars_end_max);
          ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
          ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_stars_beg_max = max(ti_beg, ti_stars_beg_max);
          ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
        }
      }
    }

    /* Loop over the sink particles in this cell. */
    for (int k = 0; k < sink_count; k++) {

      /* Get a handle on the part. */
      struct sink *restrict sink = &sinks[k];

      /* need to be updated ? */
      if (sink_is_active(sink, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Current end of time-step */
        const integertime_t ti_end =
            get_integer_time_end(ti_current, sink->time_bin);

        if (ti_end != ti_current)
          error("Computing time-step of rogue particle.");
#endif
        /* Get new time-step */
        const integertime_t ti_new_step = get_sink_timestep(sink, e);

        /* Update particle */
        sink->time_bin = get_time_bin(ti_new_step);
        sink->gpart->time_bin = get_time_bin(ti_new_step);

        /* Number of updated sink-particles */
        sink_updated++;
        g_updated++;

        ti_sinks_end_min = min(ti_current + ti_new_step, ti_sinks_end_min);
        ti_sinks_end_max = max(ti_current + ti_new_step, ti_sinks_end_max);
        ti_gravity_end_min = min(ti_current + ti_new_step, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_current + ti_new_step, ti_gravity_end_max);

        /* What is the next starting point for this cell ? */
        ti_sinks_beg_max = max(ti_current, ti_sinks_beg_max);
        ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

        /* sink particle is inactive but not inhibited */
      } else {

        if (!sink_is_inhibited(sink, e)) {

          const integertime_t ti_end =
              get_integer_time_end(ti_current, sink->time_bin);

          const integertime_t ti_beg =
              get_integer_time_begin(ti_current + 1, sink->time_bin);

          ti_sinks_end_min = min(ti_end, ti_sinks_end_min);
          ti_sinks_end_max = max(ti_end, ti_sinks_end_max);
          ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
          ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_sinks_beg_max = max(ti_beg, ti_sinks_beg_max);
          ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
        }
      }
    }

    /* Loop over the black hole particles in this cell. */
    for (int k = 0; k < bcount; k++) {

      /* Get a handle on the part. */
      struct bpart *restrict bp = &bparts[k];

      /* need to be updated ? */
      if (bpart_is_active(bp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Current end of time-step */
        const integertime_t ti_end =
            get_integer_time_end(ti_current, bp->time_bin);

        if (ti_end != ti_current)
          error("Computing time-step of rogue particle.");
#endif
        /* Get new time-step */
        const integertime_t ti_new_step = get_bpart_timestep(bp, e);

        /* Update particle */
        bp->time_bin = get_time_bin(ti_new_step);
        bp->gpart->time_bin = get_time_bin(ti_new_step);

        /* Number of updated s-particles */
        b_updated++;
        g_updated++;

        ti_black_holes_end_min =
            min(ti_current + ti_new_step, ti_black_holes_end_min);
        ti_black_holes_end_max =
            max(ti_current + ti_new_step, ti_black_holes_end_max);
        ti_gravity_end_min = min(ti_current + ti_new_step, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_current + ti_new_step, ti_gravity_end_max);

        /* What is the next starting point for this cell ? */
        ti_black_holes_beg_max = max(ti_current, ti_black_holes_beg_max);
        ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

        /* star particle is inactive but not inhibited */
      } else {

        if (!bpart_is_inhibited(bp, e)) {

          const integertime_t ti_end =
              get_integer_time_end(ti_current, bp->time_bin);

          const integertime_t ti_beg =
              get_integer_time_begin(ti_current + 1, bp->time_bin);

          ti_black_holes_end_min = min(ti_end, ti_black_holes_end_min);
          ti_black_holes_end_max = max(ti_end, ti_black_holes_end_max);
          ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
          ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_black_holes_beg_max = max(ti_beg, ti_black_holes_beg_max);
          ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
        }
      }
    }

  } else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        /* Recurse */
        runner_do_timestep(r, cp, /*timer=*/0);

        /* And aggregate */
        updated += cp->hydro.updated;
        g_updated += cp->grav.updated;
        sink_updated += cp->sinks.updated;
        s_updated += cp->stars.updated;
        b_updated += cp->black_holes.updated;

        ti_hydro_end_min = min(cp->hydro.ti_end_min, ti_hydro_end_min);
        ti_hydro_end_max = max(cp->hydro.ti_end_max, ti_hydro_end_max);
        ti_hydro_beg_max = max(cp->hydro.ti_beg_max, ti_hydro_beg_max);

        ti_gravity_end_min = min(cp->grav.ti_end_min, ti_gravity_end_min);
        ti_gravity_end_max = max(cp->grav.ti_end_max, ti_gravity_end_max);
        ti_gravity_beg_max = max(cp->grav.ti_beg_max, ti_gravity_beg_max);

        ti_stars_end_min = min(cp->stars.ti_end_min, ti_stars_end_min);
        ti_stars_end_max = max(cp->stars.ti_end_max, ti_stars_end_max);
        ti_stars_beg_max = max(cp->stars.ti_beg_max, ti_stars_beg_max);

        ti_sinks_end_min = min(cp->sinks.ti_end_min, ti_sinks_end_min);
        ti_sinks_end_max = max(cp->sinks.ti_end_max, ti_sinks_end_max);
        ti_sinks_beg_max = max(cp->sinks.ti_beg_max, ti_sinks_beg_max);

        ti_black_holes_end_min =
            min(cp->black_holes.ti_end_min, ti_black_holes_end_min);
        ti_black_holes_end_max =
            max(cp->grav.ti_end_max, ti_black_holes_end_max);
        ti_black_holes_beg_max =
            max(cp->grav.ti_beg_max, ti_black_holes_beg_max);
      }
    }
  }

  /* Store the values. */
  c->hydro.updated = updated;
  c->grav.updated = g_updated;
  c->stars.updated = s_updated;
  c->sinks.updated = sink_updated;
  c->black_holes.updated = b_updated;

  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_end_max = ti_hydro_end_max;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_end_max = ti_gravity_end_max;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_end_max = ti_stars_end_max;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_end_max = ti_sinks_end_max;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_end_max = ti_black_holes_end_max;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.ti_end_min == e->ti_current &&
      c->hydro.ti_end_min < max_nr_timesteps)
    error("End of next hydro step is current time!");
  if (c->grav.ti_end_min == e->ti_current &&
      c->grav.ti_end_min < max_nr_timesteps)
    error("End of next gravity step is current time!");
  if (c->stars.ti_end_min == e->ti_current &&
      c->stars.ti_end_min < max_nr_timesteps)
    error("End of next stars step is current time!");
  if (c->sinks.ti_end_min == e->ti_current &&
      c->sinks.ti_end_min < max_nr_timesteps)
    error("End of next sinks step is current time!");
  if (c->black_holes.ti_end_min == e->ti_current &&
      c->black_holes.ti_end_min < max_nr_timesteps)
    error("End of next black holes step is current time!");
#endif

  if (timer) TIMER_TOC(timer_timestep);
}

/**
 * @brief Apply the time-step limiter to all awaken particles in a cell
 * hierarchy.
 *
 * @param r The task #runner.
 * @param c The #cell.
 * @param force Limit the particles irrespective of the #cell flags.
 * @param timer Are we timing this ?
 */
void runner_do_limiter(struct runner *r, struct cell *c, int force,
                       const int timer) {

  const struct engine *e = r->e;
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only limit local cells. */
  if (c->nodeID != engine_rank) error("Limiting dt of a foreign cell is nope.");
#endif

  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;

  /* Limit irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_hydro_limiter));

  /* Early abort? */
  if (c->hydro.count == 0) {

    /* Clear the limiter flags. */
    cell_clear_flag(
        c, cell_flag_do_hydro_limiter | cell_flag_do_hydro_sub_limiter);
    return;
  }

  /* Loop over the progeny ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_hydro_sub_limiter))) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        /* Recurse */
        runner_do_limiter(r, cp, force, /*timer=*/0);

        /* And aggregate */
        ti_hydro_end_min = min(cp->hydro.ti_end_min, ti_hydro_end_min);
        ti_hydro_end_max = max(cp->hydro.ti_end_max, ti_hydro_end_max);
        ti_hydro_beg_max = max(cp->hydro.ti_beg_max, ti_hydro_beg_max);
        ti_gravity_end_min = min(cp->grav.ti_end_min, ti_gravity_end_min);
        ti_gravity_end_max = max(cp->grav.ti_end_max, ti_gravity_end_max);
        ti_gravity_beg_max = max(cp->grav.ti_beg_max, ti_gravity_beg_max);
      }
    }

    /* Store the updated values */
    c->hydro.ti_end_min = min(c->hydro.ti_end_min, ti_hydro_end_min);
    c->hydro.ti_end_max = max(c->hydro.ti_end_max, ti_hydro_end_max);
    c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, ti_hydro_beg_max);
    c->grav.ti_end_min = min(c->grav.ti_end_min, ti_gravity_end_min);
    c->grav.ti_end_max = max(c->grav.ti_end_max, ti_gravity_end_max);
    c->grav.ti_beg_max = max(c->grav.ti_beg_max, ti_gravity_beg_max);

  } else if (!c->split && force) {

    ti_hydro_end_min = c->hydro.ti_end_min;
    ti_hydro_end_max = c->hydro.ti_end_max;
    ti_hydro_beg_max = c->hydro.ti_beg_max;
    ti_gravity_end_min = c->grav.ti_end_min;
    ti_gravity_end_max = c->grav.ti_end_max;
    ti_gravity_beg_max = c->grav.ti_beg_max;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Avoid inhibited particles */
      if (part_is_inhibited(p, e)) continue;

      /* Bip, bip, bip... wake-up time */
      if (p->limiter_data.wakeup != time_bin_not_awake) {

        // message("Limiting particle %lld in cell %d", p->id, c->cellID);

        /* Apply the limiter and get the new end of time-step */
        const integertime_t ti_end_new = timestep_limit_part(p, xp, e);
        const timebin_t new_bin = p->time_bin;
        const integertime_t ti_beg_new =
            ti_end_new - get_integer_timestep(new_bin);

        /* Mark this particle has not needing synchronization */
        p->limiter_data.to_be_synchronized = 0;

        /* What is the next sync-point ? */
        ti_hydro_end_min = min(ti_end_new, ti_hydro_end_min);
        ti_hydro_end_max = max(ti_end_new, ti_hydro_end_max);

        /* What is the next starting point for this cell ? */
        ti_hydro_beg_max = max(ti_beg_new, ti_hydro_beg_max);

        /* Also limit the gpart counter-part */
        if (p->gpart != NULL) {

          /* Register the time-bin */
          p->gpart->time_bin = p->time_bin;

          /* What is the next sync-point ? */
          ti_gravity_end_min = min(ti_end_new, ti_gravity_end_min);
          ti_gravity_end_max = max(ti_end_new, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_beg_new, ti_gravity_beg_max);
        }
      }
    }

    /* Store the updated values */
    c->hydro.ti_end_min = min(c->hydro.ti_end_min, ti_hydro_end_min);
    c->hydro.ti_end_max = max(c->hydro.ti_end_max, ti_hydro_end_max);
    c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, ti_hydro_beg_max);
    c->grav.ti_end_min = min(c->grav.ti_end_min, ti_gravity_end_min);
    c->grav.ti_end_max = max(c->grav.ti_end_max, ti_gravity_end_max);
    c->grav.ti_beg_max = max(c->grav.ti_beg_max, ti_gravity_beg_max);
  }

  /* Clear the limiter flags. */
  cell_clear_flag(c,
                  cell_flag_do_hydro_limiter | cell_flag_do_hydro_sub_limiter);

  if (timer) TIMER_TOC(timer_do_limiter);
}

/**
 * @brief Apply the time-step synchronization proceduere to all flagged
 * particles in a cell hierarchy.
 *
 * @param r The task #runner.
 * @param c The #cell.
 * @param force Limit the particles irrespective of the #cell flags.
 * @param timer Are we timing this ?
 */
void runner_do_sync(struct runner *r, struct cell *c, int force,
                    const int timer) {

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only sync local cells. */
  if (c->nodeID != engine_rank) error("Syncing of a foreign cell is nope.");
#endif

  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;

  /* Limit irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_hydro_sync));

  /* Early abort? */
  if (c->hydro.count == 0) {

    /* Clear the sync flags. */
    cell_clear_flag(c, cell_flag_do_hydro_sync | cell_flag_do_hydro_sub_sync);
    return;
  }

  /* Loop over the progeny ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_hydro_sub_sync))) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        /* Recurse */
        runner_do_sync(r, cp, force, /*timer=*/0);

        /* And aggregate */
        ti_hydro_end_min = min(cp->hydro.ti_end_min, ti_hydro_end_min);
        ti_hydro_end_max = max(cp->hydro.ti_end_max, ti_hydro_end_max);
        ti_hydro_beg_max = max(cp->hydro.ti_beg_max, ti_hydro_beg_max);
        ti_gravity_end_min = min(cp->grav.ti_end_min, ti_gravity_end_min);
        ti_gravity_end_max = max(cp->grav.ti_end_max, ti_gravity_end_max);
        ti_gravity_beg_max = max(cp->grav.ti_beg_max, ti_gravity_beg_max);
      }
    }

    /* Store the updated values */
    c->hydro.ti_end_min = min(c->hydro.ti_end_min, ti_hydro_end_min);
    c->hydro.ti_end_max = max(c->hydro.ti_end_max, ti_hydro_end_max);
    c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, ti_hydro_beg_max);
    c->grav.ti_end_min = min(c->grav.ti_end_min, ti_gravity_end_min);
    c->grav.ti_end_max = max(c->grav.ti_end_max, ti_gravity_end_max);
    c->grav.ti_beg_max = max(c->grav.ti_beg_max, ti_gravity_beg_max);

  } else if (!c->split && force) {

    ti_hydro_end_min = c->hydro.ti_end_min;
    ti_hydro_end_max = c->hydro.ti_end_max;
    ti_hydro_beg_max = c->hydro.ti_beg_max;
    ti_gravity_end_min = c->grav.ti_end_min;
    ti_gravity_end_max = c->grav.ti_end_max;
    ti_gravity_beg_max = c->grav.ti_beg_max;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Avoid inhibited particles */
      if (part_is_inhibited(p, e)) continue;

      /* If the particle is active no need to sync it */
      if (part_is_active(p, e) && p->limiter_data.to_be_synchronized) {
        p->limiter_data.to_be_synchronized = 0;
      }

      if (p->limiter_data.to_be_synchronized) {

        /* Finish this particle's time-step */
        timestep_process_sync_part(p, xp, e, cosmo);

        /* Get new time-step */
        integertime_t ti_new_step = get_part_timestep(p, xp, e);
        timebin_t new_time_bin = get_time_bin(ti_new_step);

        /* Limit the time-bin to what is allowed in this step */
        new_time_bin = min(new_time_bin, e->max_active_bin);
        ti_new_step = get_integer_timestep(new_time_bin);

        /* Update particle */
        p->time_bin = new_time_bin;
        if (p->gpart != NULL) p->gpart->time_bin = new_time_bin;

        /* Update the tracers properties */
        tracers_after_timestep(p, xp, e->internal_units, e->physical_constants,
                               with_cosmology, e->cosmology,
                               e->hydro_properties, e->cooling_func, e->time);

        /* What is the next sync-point ? */
        ti_hydro_end_min = min(ti_current + ti_new_step, ti_hydro_end_min);
        ti_hydro_end_max = max(ti_current + ti_new_step, ti_hydro_end_max);

        /* What is the next starting point for this cell ? */
        ti_hydro_beg_max = max(ti_current, ti_hydro_beg_max);

        /* Also limit the gpart counter-part */
        if (p->gpart != NULL) {

          /* Register the time-bin */
          p->gpart->time_bin = p->time_bin;

          /* What is the next sync-point ? */
          ti_gravity_end_min =
              min(ti_current + ti_new_step, ti_gravity_end_min);
          ti_gravity_end_max =
              max(ti_current + ti_new_step, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);
        }
      }
    }

    /* Store the updated values */
    c->hydro.ti_end_min = min(c->hydro.ti_end_min, ti_hydro_end_min);
    c->hydro.ti_end_max = max(c->hydro.ti_end_max, ti_hydro_end_max);
    c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, ti_hydro_beg_max);
    c->grav.ti_end_min = min(c->grav.ti_end_min, ti_gravity_end_min);
    c->grav.ti_end_max = max(c->grav.ti_end_max, ti_gravity_end_max);
    c->grav.ti_beg_max = max(c->grav.ti_beg_max, ti_gravity_beg_max);
  }

  /* Clear the sync flags. */
  cell_clear_flag(c, cell_flag_do_hydro_sync | cell_flag_do_hydro_sub_sync);

  if (timer) TIMER_TOC(timer_do_sync);
}
