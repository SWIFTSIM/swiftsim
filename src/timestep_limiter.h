/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TIMESTEP_LIMITER_H
#define SWIFT_TIMESTEP_LIMITER_H

/* Config parameters. */
#include "../config.h"

__attribute__((always_inline)) INLINE static integertime_t timestep_limit_part(
    struct part *restrict p, struct xpart *restrict xp,
    const struct engine *e) {

  /* Let's see when this particle started and used to end */
  const integertime_t old_ti_beg =
      get_integer_time_begin(e->ti_current, p->time_bin);
  const integertime_t old_ti_end =
      get_integer_time_end(e->ti_current, p->time_bin);
  const integertime_t old_dti = old_ti_end - old_ti_beg;

  /* The new time-step */
  const integertime_t new_fake_ti_step =
      get_integer_timestep(e->min_active_bin);
  const integertime_t new_ti_beg = old_ti_beg;
  const integertime_t new_ti_end = e->ti_current + new_fake_ti_step;
  const integertime_t new_dti = new_ti_end - new_ti_beg;

#ifdef SWIFT_DEBUG_CHECKS
  /* Some basic safety checks */
  if (old_ti_beg >= e->ti_current)
    error(
        "Incorrect value for old time-step beginning ti_current=%lld, "
        "old_ti_beg=%lld",
        e->ti_current, old_ti_beg);
  if (old_ti_end <= e->ti_current)
    error(
        "Incorrect value for old time-step end ti_current=%lld, "
        "old_ti_end=%lld",
        e->ti_current, old_ti_end);
  if (new_ti_end > old_ti_end)
    error("New end of time-step longer than old one");
  if (new_dti > old_dti) error("New time-step larger than old one");
  if (new_fake_ti_step == 0) error("Wakeup call too early");
#endif

  /* if (p->id == ICHECK) */
  /*   message( */
  /*       "ti_current=%lld, ti_drift=%lld ti_kick=%lld, time_bin=%d wakeup=%d " */
  /*       "new_ti_end=%lld", */
  /*       e->ti_current, p->ti_drift, p->ti_kick, p->time_bin, p->wakeup, */
  /*       new_ti_end); */

  /* Now we need to reverse the kick1...*/
  kick_part(p, xp, old_ti_beg + old_dti / 2, old_ti_beg, e->timeBase);

  /* if (p->id == ICHECK) */
  /*   message( */
  /*       "ti_current=%lld, ti_drift=%lld ti_kick=%lld, time_bin=%d wakeup=%d " */
  /*       "new_ti_end=%lld", */
  /*       e->ti_current, p->ti_drift, p->ti_kick, p->time_bin, p->wakeup, */
  /*       new_ti_end); */

  /* ...and apply the new one */
  kick_part(p, xp, new_ti_beg, new_ti_beg + new_dti / 2, e->timeBase);

  /* if (p->id == ICHECK) */
  /*   message( */
  /*       "ti_current=%lld, ti_drift=%lld ti_kick=%lld, time_bin=%d wakeup=%d", */
  /*       e->ti_current, p->ti_drift, p->ti_kick, p->time_bin, p->wakeup); */

  /* Remember the old time-bin */
  p->wakeup = p->time_bin;

  /* Update the time bin of this particle */
  p->time_bin = e->min_active_bin;

  /* if (p->id == ICHECK) */
  /*   message( */
  /*       "ti_current=%lld, ti_drift=%lld ti_kick=%lld, time_bin=%d wakeup=%d", */
  /*       e->ti_current, p->ti_drift, p->ti_kick, p->time_bin, p->wakeup); */

  return new_fake_ti_step;
}

#endif /* SWIFT_TIMESTEP_LIMITER_H */
