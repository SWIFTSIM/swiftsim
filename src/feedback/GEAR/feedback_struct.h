/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_STRUCT_GEAR_H
#define SWIFT_FEEDBACK_STRUCT_GEAR_H

#include "chemistry_struct.h"

/**
 * @brief Feedback fields carried by each hydro particles
 */

struct feedback_part_data {
  /*! mass received from supernovae */
  float delta_mass;

  /*! specific energy received from supernovae */
  float delta_u;
};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  /*! Inverse of normalisation factor used for the enrichment. */
  float enrichment_weight;

  union {

    /**
     * @brief Values collected from the gas neighbours.
     */
    struct {


    } to_collect;

    /**
     * @brief Values to be distributed to the gas neighbours.
     */
    struct {

      /*! Mass released. */
      float mass;

      /*! Energy change due to ejectas. */
      float energy;

    } to_distribute;
  };

    /* Does this star already exploded */
  integertime_t explosion_time;

};

#endif /* SWIFT_FEEDBACK_STRUCT_GEAR_H */
