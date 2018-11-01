/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_LOGGER_H
#define SWIFT_LOGGER_H

#ifdef WITH_LOGGER

/* Includes. */
#include "common_io.h"
#include "inline.h"
#include "timeline.h"
#include "units.h"

/* Forward declaration */
struct dump;
struct gpart;
struct part;
/* TODO remove dependency */
struct engine;

/**
 * Logger entries contain messages representing the particle data at a given
 * point in time during the simulation.
 *
 * The logger messages always start with an 8-byte header structured as
 * follows:
 *
 *   data: [ mask |                     offset                     ]
 *   byte: [  01  |  02  |  03  |  04  |  05  |  06  |  07  |  08  ]
 *
 * I.e. a first "mask" byte followed by 7 "offset" bytes. The mask contains
 * information on what kind of data is packed after the header. The mask
 * bits correspond to the following data:
 *
 *   bit | name   | size | comment
 *   -------------------------------------------------------------------------
 *   0   | x      | 24   | The particle position, in absolute coordinates,
 *       |        |      | stored as three doubles.
 *   1   | v      | 12   | Particle velocity, stored as three floats.
 *   2   | a      | 12   | Particle acceleration, stored as three floats.
 *   3   | u      | 4    | Particle internal energy (or entropy, if Gadget-SPH
 *       |        |      | is used), stored as a single float.
 *   4   | h      | 4    | Particle smoothing length (or epsilon, if a gpart),
 *       |        |      | stored as a single float.
 *   5   | rho    | 4    | Particle density, stored as a single float.
 *   6   | consts | 12   | Particle constants, i.e. mass and ID.
 *   7   | time   | 8    | Timestamp, not associated with a particle, just
 *       |        |      | marks the transitions from one timestep to another.
 *
 * There is no distinction between gravity and SPH particles.
 *
 * The offset refers to the relative location of the previous message for the
 * same particle or for the previous timestamp (if mask bit 7 is set). I.e.
 * the previous log entry will be at the address of the current mask byte minus
 * the unsigned value stored in the offset. An offset equal to the chunk offset
 * indicated that this is the first message for the given particle/timestamp.
 */

/* Some constants. */
enum logger_masks {
  logger_mask_x = (1 << 0),
  logger_mask_v = (1 << 1),
  logger_mask_a = (1 << 2),
  logger_mask_u = (1 << 3),
  logger_mask_h = (1 << 4),
  logger_mask_rho = (1 << 5),
  logger_mask_consts = (1 << 6),
  logger_mask_timestamp = (1 << 7),
};

/* Size of the strings. */
#define logger_string_length 200

/* parameters of the logger */
struct logger_parameters {
  /* size of a label in bytes */
  size_t label_size;

  /* size of an offset in bytes */
  size_t offset_size;

  /* size of a mask in bytes */
  size_t mask_size;

  /* size of a number in bytes */
  size_t number_size;

  /* size of a data type in bytes */
  size_t data_type_size;

  /* number of different mask */
  size_t number_mask;

  /* value of each masks */
  size_t *masks;

  /* data size of each mask */
  size_t *masks_data_size;

  /* label of each mask */
  char *masks_name;

  /* Size of a chunk if every mask are activated */
  size_t total_size;
};

/* structure containing global data */
struct logger {
  /* Number of particle steps between dumping a chunk of data */
  short int delta_step;

  /* Logger basename */
  char base_name[logger_string_length];

  /* File name of the dump file */
  struct dump *dump;

  /* timestamp offset for logger*/
  size_t timestamp_offset;

  /* scaling factor when buffer is too small */
  float buffer_scale;

  /* logger parameters */
  struct logger_parameters *params;

} SWIFT_STRUCT_ALIGN;

/* required structure for each particle type */
struct logger_part_data {
  /* Number of particle updates since last output */
  int steps_since_last_output;

  /* offset of last particle log entry */
  size_t last_offset;
};

enum logger_datatype {
  logger_data_int,
  logger_data_float,
  logger_data_double,
  logger_data_char,
  logger_data_longlong,
  logger_data_bool,
  logger_data_count /* should be last */
};

extern const unsigned int logger_datatype_size[];

/* Function prototypes. */
int logger_compute_chunk_size(unsigned int mask);
void logger_log_all(struct logger *log, const struct engine *e);
void logger_log_part(struct logger *log, const struct part *p,
                     unsigned int mask, size_t *offset);
void logger_log_gpart(struct logger *log, const struct gpart *p,
                      unsigned int mask, size_t *offset);
void logger_init(struct logger *log, struct swift_params *params);
void logger_clean(struct logger *log);
void logger_log_timestamp(struct logger *log, integertime_t t, size_t *offset);
void logger_ensure_size(struct logger *log, size_t total_nr_parts,
                        size_t total_nr_gparts, size_t total_nr_sparts);
void logger_write_file_header(struct logger *log, const struct engine *e);

int logger_read_part(struct part *p, size_t *offset, const char *buff);
int logger_read_gpart(struct gpart *p, size_t *offset, const char *buff);
int logger_read_timestamp(unsigned long long int *t, size_t *offset,
                          const char *buff);

void logger_parameters_init(struct logger_parameters *log_params);
void logger_parameters_clean(struct logger_parameters *log_params);

/**
 * @brief Initialize the logger data for a particle.
 *
 * @param logger The #logger_part_data.
 */
INLINE static void logger_part_data_init(struct logger_part_data *logger) {
  logger->last_offset = 0;
  logger->steps_since_last_output = SHRT_MAX;
}

/**
 * @brief Should this particle write its data now ?
 *
 * @param xp The #xpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part should write, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int logger_should_write(
    const struct logger_part_data *logger_data, const struct logger *log) {

  return (logger_data->steps_since_last_output > log->delta_step);
}

#endif /* WITH_LOGGER */

#endif /* SWIFT_LOGGER_H */
