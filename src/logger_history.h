/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_LOGGER_HISTORY_H
#define SWIFT_LOGGER_HISTORY_H

#include "../config.h"

/* Standard includes */
#include <stdint.h>

/* Local include */
#include "error.h"
#include "part_type.h"

#if defined(WITH_LOGGER)

/* Forward declaration */
struct xpart;
struct part;
struct gpart;
struct spart;
struct bpart;
struct engine;
struct swift_params;

/**
 * @brief Contains the information concerning
 * a particle for the index files.
 */
struct logger_index_data {
  /* Id of the particle. */
  int64_t id;

  /* Offset of the particle in the file. */
  uint64_t offset;
};

/**
 * @brief Structure dealing with the changes in the number
 * of particles (e.g. creation, deletion, transformation).
 */
struct logger_history {

  /* Number of elements currently stored */
  uint64_t size[swift_type_count];

  /* Size of the current buffer */
  size_t capacity[swift_type_count];

  /* Buffer containing the particles */
  struct logger_index_data *data[swift_type_count];
};

void logger_history_first_init(struct logger_history *hist);
void logger_history_init(struct logger_history *hist);
void logger_history_clean(struct logger_history *hist);
void logger_history_log_part(struct logger_history *hist, const struct part *p,
                             const struct xpart *xp);
void logger_history_log_spart(struct logger_history *hist,
                              const struct spart *sp);
void logger_history_log_gpart(struct logger_history *hist,
                              const struct gpart *gp);
void logger_history_log_bpart(struct logger_history *hist,
                              const struct bpart *bp);
void logger_history_write(struct logger_history *hist, struct engine *e,
                          FILE *f);
size_t logger_history_get_size(const struct logger_history *hist);

void logger_history_dump(const struct logger_history *hist, FILE *stream);
void logger_history_restore(struct logger_history *hist, FILE *stream);

#endif  // WITH_LOGGER
#endif  // SWIFT_LOGGER_HISTORY_H
