/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_STRUCT_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_STRUCT_H

#include "event_logger_core.h"

/* feedback history struct for SNIa */
struct feedback_history_SNIa {

  /* Load the core of logging functions */
  struct event_history_logger core;

  /*! Total new SNIa injected energy */
  double SNIa_energy;

  /*! Number of heating events */
  int events;
};

/* feedback history struct for SNII */
struct feedback_history_SNII {

  /* Load the core of logging functions */
  struct event_history_logger core;

  /*! Total new SNIa injected energy */
  double SNII_energy;

  /*! Total new SNIas in the simulation */
  double N_SNII;

  /*! Number of heating events */
  int events;
};

/* feedback history struct for r-processes */
struct feedback_history_r_processes {

  /* Load the core of logging functions */
  struct event_history_logger core;

  /*! Total new r-processes mass in the simulation */
  double enrichment_mass;

  /*! Number of r-processes events */
  int events;
};

#ifdef SWIFT_DEBUG_CHECKS
/* feedback history debug for SNIa */
struct feedback_history_debug {

  /* The file pointer of the struct */
  FILE *fp;

  /* The lock of the debugger struct */
  swift_lock_type lock;
};
#endif /* SWIFT_DEBUG_CHECKS */

/* Global variable for the SNIa events */
extern struct feedback_history_SNIa log_SNIa;

/* Global variable for the SNII events */
extern struct feedback_history_SNII log_SNII;

/* Global variable for the r-processes */
extern struct feedback_history_r_processes log_r_processes;

#ifdef SWIFT_DEBUG_CHECKS
/* Global variable for SNIa events in debugging mode */
extern struct feedback_history_debug log_SNIa_debug;
/* Global variable for SNII events in debugging mode */
extern struct feedback_history_debug log_SNII_debug;
#endif

#ifdef SWIFT_DEBUG_CHECKS
struct feedback_history_debug log_SNIa_debug;
struct feedback_history_debug log_SNII_debug;
#endif

/* Declare the feedback structures */
struct feedback_history_SNIa log_SNIa;
struct feedback_history_SNII log_SNII;
struct feedback_history_r_processes log_r_processes;

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_STRUCT_H */
