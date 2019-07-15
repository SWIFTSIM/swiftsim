/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_LOGGER_IO_H
#define SWIFT_LOGGER_IO_H

/* Config parameters. */
#include "../config.h"

#ifdef WITH_LOGGER

/* Includes. */
#include "engine.h"
#include "io_properties.h"
#include "part.h"
#include "units.h"


void logger_write_index_file(struct logger *log, struct engine* e);
void logger_write_description(struct logger *log, struct engine* e);

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param N number of particles.
 * @param list (out) The parameters to write.
 *
 * In this version, we only want the ids and the offset.
 */
__attribute__((always_inline)) INLINE static int hydro_write_index(
    const struct part* parts, const struct xpart* xparts,
    struct io_props *list) {

  /* List what we want to write */
  list[0] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[1] = io_make_output_field("Offset", SIZE_T, 1, UNIT_CONV_NO_UNITS,
				 xparts, logger_data.last_offset);

  return 2;
}


/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param gparts The gparticle array.
 * @param N number of particles.
 * @param list (out) The parameters to write.
 *
 * In this version, we only want the ids and the offset.
 */
__attribute__((always_inline)) INLINE static int darkmatter_write_index(
    const struct gpart* gparts, struct io_props *list) {

  /* List what we want to write */
  list[0] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                UNIT_CONV_NO_UNITS, gparts, id_or_neg_offset);
  list[1] = io_make_output_field("Offset", SIZE_T, 1, UNIT_CONV_NO_UNITS,
				 gparts, logger_data.last_offset);

  return 2;
}

#endif

#endif /* SWIFT_LOGGER_IO_H */
