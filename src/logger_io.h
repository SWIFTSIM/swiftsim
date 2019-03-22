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

/* Overwrite a few HDF5 functions in order to avoid code duplication */
#define hid_t FILE *
#define io_write_attribute_s(file, params, value) {		\
    fprintf(file, "    %s: %s\n", params, value);		\
  }
#define io_write_attribute_d(file, params, value) {		\
    fprintf(file, "    %s: %g\n", params, value);		\
  }
#define io_write_attribute_f(file, params, value) {		\
    fprintf(file, "    %s: %g\n", params, value);		\
  }
#define io_write_attribute_i(file, params, value) {		\
    fprintf(file, "    %s: %i\n", params, value);		\
  }
#define io_write_attribute_l(file, params, value) {		\
    fprintf(file, "    %s: %li\n", params, value);		\
  }

#define io_write_attribute(file, params, type, value, dim) {	\
    fprintf(file, "    %s: [", params);				\
    for(int i = 0; i < dim-1; i++) {				\
      fprintf(file, type ", ", value[i]);			\
    }								\
    fprintf(file, type "]\n", value[dim-1]);			\
}

void logger_write_index_file(struct logger *log, struct engine* e);
void logger_write_description(struct logger *log, struct engine* e);

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param N number of particles.
 * @param f File to use.
 *
 * In this version, we only want the ids and the offset.
 */
__attribute__((always_inline)) INLINE static void hydro_write_index(
    const struct part* parts, const struct xpart* xparts,
    struct io_props *list) {

  /* List what we want to write */
  list[0] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[1] = io_make_output_field("Offset", SIZE_T, 1, UNIT_CONV_NO_UNITS,
				 xparts, logger_data.last_offset);
}
#endif

#endif /* SWIFT_LOGGER_IO_H */
