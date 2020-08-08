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

    /*! Nr of #sink in this cell. */
    int count;

    /*! Nr of #sink this cell can hold after addition of new one. */
    int count_total;

    /*! Is the #sink data of this cell being used in a sub-cell? */
    int hold;

    /*! Spin lock for various uses (#sink case). */
    swift_lock_type lock;

    /*! Last (integer) time the cell's sink were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Minimum end of (integer) time step in this cell for sink tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for sink tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for sink
     * tasks. */
    integertime_t ti_beg_max;

#ifdef SINK_NONE
  };
#endif
};

#endif /* SWIFT_CELL_SINKS_H */
