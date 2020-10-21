/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_M1_H
#define SWIFT_RT_M1_H

/**
 * @file src/rt/M1closure/rt.h
 * @brief Main header file for the M1 Closure radiative transfer scheme.
 */

/**
 * @brief Update the photon number of a particle, i.e. compute
 *        E^{n+1} = E^n + dt * dE_* / dt
 */
__attribute__((always_inline)) INLINE static void
rt_injection_update_photon_density(struct part* restrict p) {}

/**
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is initialized
 *        and assumes that the photon emission rate is an intrinsic
 *        stellar property, i.e. doesn't depend on the environment.
 *
 * @param sp star particle to work on
 * @param cosmo cosmology struct in use
 * @param with_cosmology whether we're running with cosmology
 * @param ti_current current system integer time
 * @param time current system time
 * @param time_base time base in use
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(struct spart* restrict sp,
                                 const struct cosmology* cosmo,
                                 int with_cosmology,
                                 const integertime_t ti_current, double time,
                                 double time_base) {}

/**
 * @brief First initialisation of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p) {}

/**
 * @brief Initialisation of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {}

/**
 * @brief First initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

#endif /* SWIFT_RT_M1_H */
