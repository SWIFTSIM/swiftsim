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

#include <math.h>
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
  const double norm = fp->dtd_data.norm_const;
  const double nu_gauss = fp->dtd_data.SNIa_efficiency_gauss;
  const double inv_std = fp->dtd_data.std_characteristic_time_Gyr_inv;
  const double tdif0 = t0 - fp->dtd_data.characteristic_time_Gyr;
  const double tdif1 = t1 - fp->dtd_data.characteristic_time_Gyr;
  const double num_SNIa_per_Msun = norm * (t1 - t0) + .5 * nu_gauss * ( erf(tdif1 * inv_std) - erf(tdif0 * inv_std));

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

  /* Set the SNIa efficiency for the constant part of the DTD */
  fp->dtd_data.SNIa_efficiency_const = parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_const_p_Msun");

  /* Set the SNIa efficiency for the Gaussian part of the DTD */
  fp->dtd_data.SNIa_efficiency_gauss = parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_gauss_p_Msun");

  /* Set the normalization time for the constant part of the DTD */
  fp->dtd_data.normalization_timescale_Gyr = parser_get_param_double(params, "SNIaDTD:normalization_timescale_Gyr");

  /* Set the delay time of the DTD */
  fp->dtd_data.delay_time_Gyr = parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  /* Set the characteristic time of the Gaussian part of the DTD */
  fp->dtd_data.characteristic_time_Gyr = parser_get_param_double(params, "SNIaDTD:characteristic_time_Gyr"); 

  /* Set the standard deviation of the Gaussian part of the DTD */
  fp->dtd_data.std_characteristic_time_Gyr = parser_get_param_double(params, "SNIaDTD:STD_characteristic_time_Gyr");

  /* Calculate the inverse of the standard deviation of the Gaussian part divided by sqrt 2 */
  fp->dtd_data.std_characteristic_time_Gyr_inv = 1. / fp->dtd_data.std_characteristic_time_Gyr / sqrt(2);
  
  /* Calculate the norm of the constant part of the DTD */
  fp->dtd_data.norm_const = fp->dtd_data.SNIa_efficiency_const/fp->dtd_data.normalization_timescale_Gyr;
}

