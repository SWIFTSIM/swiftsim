/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Josh Borrow (josh@joshborrow.com)
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

/**
 * @file particle_splitting.c
 * @brief Particle splitting functions, used for tracking and
 *        outputting split data.
 */

/* Local headers. */
#include "particle_splitting.h"

/**
 * @brief Initialise a particle_splitting_data struct
 *        at the start of a run, given an initial
 *        progenitor ID.
 *
 * @param splitting_data the uninitialised particle struct.
 * @param id the ID of the particle (used for progenitor_id)
 */
__attribute__((always_inline)) INLINE static void
particle_splitting_mark_part_as_not_split(
    struct particle_splitting_data *restrict splitting_data, int id) {

  splitting_data->progenitor_id = id;
  splitting_data->split_tree = 0;
  splitting_data->split_count = 0;

  return;
}