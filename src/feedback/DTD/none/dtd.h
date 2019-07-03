/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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

#include "physical_constants.h"
#include "feedback_properties.h"
#include "parser.h"
#include "units.h"
#include "snia_dtd_struct.h"

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * This model assumes that the SNIa DTD is zero
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param fp The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                      const double t1,
                                      const struct feedback_props* fp) {

  return 0.;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {}

