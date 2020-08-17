#ifndef SWIFT_CELL_BLACK_HOLES_H
#define SWIFT_CELL_BLACK_HOLES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "task.h"
#include "timeline.h"

struct cell_black_holes {

#ifdef BLACK_HOLES_NONE
  union {
#endif

    /*! Pointer to the #bpart data. */
    struct bpart *parts;

    /*! The drift task for bparts */
    struct task *drift;

    /*! Implicit tasks marking the entry of the BH physics block of tasks
     */
    struct task *black_holes_in;

    /*! Implicit tasks marking the exit of the BH physics block of tasks */
    struct task *black_holes_out;

    /*! The density ghost task itself */
    struct task *density_ghost;

    /*! The other ghost tasks themselves */
    struct task *swallow_ghost_0;
    struct task *swallow_ghost_1;
    struct task *swallow_ghost_2;

    /*! Linked list of the tasks computing this cell's BH density. */
    struct link *density;

    /*! Linked list of the tasks computing this cell's BH swallowing and
     * merging. */
    struct link *swallow;

    /*! Linked list of the tasks processing the particles to swallow */
    struct link *do_gas_swallow;

    /*! Linked list of the tasks processing the particles to swallow */
    struct link *do_bh_swallow;

    /*! Linked list of the tasks computing this cell's BH feedback. */
    struct link *feedback;

    /*! Last (integer) time the cell's bpart were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Maximum end of (integer) time step in this cell for black hole tasks. */
    integertime_t ti_end_max;

    /*! Max smoothing length in this cell. */
    float h_max;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

#ifdef BLACK_HOLES_NONE
  };
#endif
  
  /*! Maximum end of (integer) time step in this cell for black tasks. */
  integertime_t ti_end_min;
    
  /*! Maximum beginning of (integer) time step in this cell for black hole
   * tasks. */
  integertime_t ti_beg_max;

  /*! Spin lock for various uses (#bpart case). */
  swift_lock_type lock;
  
  /*! Nr of #bpart this cell can hold after addition of new #bpart. */
  int count_total;
  
  /*! Number of #bpart updated in this cell. */
  int updated;
  
  /*! Is the #bpart data of this cell being used in a sub-cell? */
  int hold;
  
  /*! Nr of #bpart in this cell. */
  int count;
};

#endif /* SWIFT_CELL_BLACK_HOLES_H */
