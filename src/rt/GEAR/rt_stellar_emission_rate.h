/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STELLAR_EMISSION_RATE_GEAR_H
#define SWIFT_RT_STELLAR_EMISSION_RATE_GEAR_H

/**
 * @file src/rt/GEAR/rt_stellar_emission_rate.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * stellar radiation emission rates related functions.
 */

/**
 * @brief Main function for getting the stellar emission rate and updating it
 * in the spart.
 *
 * @param sp Star particle to work on.
 * @param age_beg Age of the stars at the beginning of the step
 * @param age_end Age of the stars at the end of the step
 * @param props RT props struct
 */
__attribute__((always_inline)) INLINE static void rt_set_stellar_emission_rate(
    struct spart *restrict sp, double age_beg, double age_end, struct rt_props* props) {

  const double dt = age_end - age_beg;
  const float dtf = (float) dt;
  if (props->use_const_emission_rates) {
    for (int g = 0; g < RT_NGROUPS; g++) {
      sp->rt_data.emission_this_step[g] = props->stellar_const_emission_rates[g] * dtf;
    }
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  sp->rt_data.debug_emission_rate_set += 1;
#endif
}

#endif /* SWIFT_RT_STELLAR_EMISSION_RATE_GEAR_H */
