/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

#ifndef SWIFT_SPACE_UNIQUE_ID_H
#define SWIFT_SPACE_UNIQUE_ID_H

/* Config parameters. */
#include "../config.h"

/* Predefine the space structure */
struct space;

/**
 * @brief Slice of unique IDs for particle creation.
 */
struct slice {
  /*! Current free unique id */
  long long current;

  /*! Maximal unique id in this slice  (not included) */
  long long max;
};

void space_update_unique_id(struct space *s);
long long space_get_new_unique_id(struct space *s);

#endif  // SWIFT_SPACE_UNIQUE_ID_H
