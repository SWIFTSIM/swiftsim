/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CELL_SINKS_H
#define SWIFT_CELL_SINKS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "timeline.h"

struct cell_sinks {

#ifdef SINK_NONE
  union {
#endif

    /*! Pointer to the #sink data. */
    struct sink *parts;

    /*! The drift task for sinks */
    struct task *drift;

    /*! Last (integer) time the cell's sink were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Maximum end of (integer) time step in this cell for sink tasks. */
    integertime_t ti_end_max;

    /*! Nr of #sink this cell can hold after addition of new one. */
    int count_total;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

#ifdef SINK_NONE
  };
#endif

  /*! Minimum end of (integer) time step in this cell for sink tasks. */
  integertime_t ti_end_min;

  /*! Maximum beginning of (integer) time step in this cell for sink
   * tasks. */
  integertime_t ti_beg_max;

  /*! Spin lock for various uses (#sink case). */
  swift_lock_type lock;

  /*! Number of #sink updated in this cell. */
  int updated;

  /*! Is the #sink data of this cell being used in a sub-cell? */
  int hold;

  /*! Nr of #sink in this cell. */
  int count;
};

#endif /* SWIFT_CELL_SINKS_H */
