/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_GEAR_RT_PARAMETERS_H
#define SWIFT_GEAR_RT_PARAMETERS_H

/**
 * @file src/rt/GEAR/rt_parameters.h
 * @brief Global RT parameters.
 */

extern struct rt_parameters rt_params;

/**
 * Some global RT related parameters.
 */
struct rt_parameters {
  float reduced_speed_of_light;
  float reduced_speed_of_light_inverse;
};


#endif /* SWIFT_GEAR_RT_PARAMETERS_H */
