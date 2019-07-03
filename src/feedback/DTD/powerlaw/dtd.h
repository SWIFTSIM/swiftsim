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
 * This model assumes that the SNIa DTD is constant
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                      const double t1,
                                      const struct feedback_props* fp) {

/* The calculation is written as the integral between t0 and t1 of
   * a constant DTD given by \nu / \tau */
  const double norm = fp->dtd_data.norm;
  const double power = fp->dtd_data.power;
  const double num_SNIa_per_Msun = norm * (pow(t1,1. - power) - pow(t0, 1. - power));

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
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
    const struct unit_system* us, struct swift_params* params) {

  fp->dtd_data.SNIa_efficiency = parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  fp->dtd_data.SNIa_timescale_Gyr = parser_get_param_float(params, "SNIaDTD:SNIa_timescale_Gyr");

  fp->dtd_data.power = parser_get_param_double(params, "SNIaDTD:power_law_slope");

  fp->dtd_data.normalization_timescale_Gyr = parser_get_param_double(params, "SNIaDTD:Normalization_timescale_Gyr");

  fp->dtd_data.delay_time_Gyr = parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  const double below_frac = pow(fp->dtd_data.normalization_timescale_Gyr,1-fp->dtd_data.power) - pow(fp->dtd_data.delay_time_Gyr,1-fp->dtd_data.power);
 
  fp->dtd_data.norm = fp->dtd_data.SNIa_efficiency / below_frac;
  
}

