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

/*
double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                      const double t1,
                                      const struct feedback_props* props);

void dtd_init(struct feedback_props* fp, const struct phys_const* phys_const,
     const struct unit_system* us, struct swift_params* params);
*/
double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                      const double t1,
                                      const struct feedback_props* props) {

/* The calculation is written as the integral between t0 and t1 of
   * a constant DTD given by \nu / \tau */
  const double tau_inv = props->dtd_data.SNIa_timescale_Gyr_inv;
  const double nu = props->dtd_data.SNIa_efficiency;
  const double num_SNIa_per_Msun = nu * (t1 - t0) * tau_inv;

  return num_SNIa_per_Msun * sp->mass_init * props->mass_to_solar_mass;
}

void dtd_init(struct feedback_props* fp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params) {

  fp->dtd_data.SNIa_efficiency = parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  fp->dtd_data.SNIa_timescale_Gyr = parser_get_param_float(params, "SNIaDTD:SNIa_timescale_Gyr");

  fp->dtd_data.SNIa_timescale_Gyr_inv = 1.f / fp->dtd.SNIa_timescale_Gyr;

}

