/****************************************************************************
 * This file is part of CHIMES.
 * Copyright (c) 2020 Alexander Richings (alexander.j.richings@durham.ac.uk)
 *
 * CHIMES is free software: you can redistribute it and/or modify
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
 ***************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "chimes_interpol.h"

/* Define CHIMES-specific maths operations, depending
 * on whether we are using single or double precision.
 * Also, the chimes_exp10() function may use the
 * exp10() or exp10f() function if the GNU extensions
 * are enabled, otherwise we will need to use our
 * own definition. */
#ifdef CHIMES_USE_DOUBLE_PRECISION
#define CHIMES_FLT_MIN DBL_MIN

#ifdef CHIMES_ENABLE_GNU_SOURCE
__attribute__((always_inline, const)) inline static double chimes_exp10(
    const double x) {
  return exp10(x);
}
#else
__attribute__((always_inline, const)) inline static double chimes_exp10(
    const double x) {
  return exp(x * M_LN10);
}
#endif

__attribute__((always_inline, const)) inline static double chimes_exp(
    const double x) {
  return exp(x);
}

__attribute__((always_inline, const)) inline static double chimes_pow(
    const double x, const double y) {
  return pow(x, y);
}

__attribute__((always_inline, const)) inline static double chimes_log10(
    const double x) {
  return log10(x);
}

__attribute__((always_inline, const)) inline static double chimes_log(
    const double x) {
  return log(x);
}

__attribute__((always_inline, const)) inline static double chimes_sqrt(
    const double x) {
  return sqrt(x);
}

#else  // CHIMES_USE_DOUBLE_PRECISION
#define CHIMES_FLT_MIN FLT_MIN

#ifdef CHIMES_ENABLE_GNU_SOURCE
__attribute__((always_inline, const)) inline static float chimes_exp10(
    const float x) {
  return exp10f(x);
}
#else   // CHIMES_ENABLE_GNU_SOURCE
__attribute__((always_inline, const)) inline static float chimes_exp10(
    const float x) {
  return expf(x * (float)M_LN10);
}
#endif  // CHIMES_ENABLE_GNU_SOURCE

__attribute__((always_inline, const)) inline static float chimes_exp(
    const float x) {
  return expf(x);
}

__attribute__((always_inline, const)) inline static float chimes_pow(
    const float x, const float y) {
  return powf(x, y);
}

__attribute__((always_inline, const)) inline static float chimes_log10(
    const float x) {
  return log10f(x);
}

__attribute__((always_inline, const)) inline static float chimes_log(
    const float x) {
  return logf(x);
}

__attribute__((always_inline, const)) inline static float chimes_sqrt(
    const float x) {
  return sqrtf(x);
}
#endif  // CHIMES_USE_DOUBLE_PRECISION

/* In some cases, we need an exp10() function
 * that is always in double precision, even
 * when we are running the rest of CHIMES
 * in single precision. This is used in the
 * photoheating shielding functions. */
#ifdef CHIMES_ENABLE_GNU_SOURCE
__attribute__((always_inline, const)) inline static double chimes_exp10_dbl(
    const double x) {
  return exp10(x);
}
#else   // CHIMES_ENABLE_GNU_SOURCE
__attribute__((always_inline, const)) inline static double chimes_exp10_dbl(
    const double x) {
  return exp(x * M_LN10);
}
#endif  // CHIMES_ENABLE_GNU_SOURCE

/*!< Electron mass. Units: g. */
#define ELECTRON_MASS 9.10938356e-28

/*!< Proton mass. Units: g. */
#define PROTON_MASS 1.672621898e-24

/*!< Pi. */
#define PI 3.1415927

/*!< Speed of light. Units: cm / s. */
#define LIGHTSPEED 2.99792458e10

/*!< Boltzmann constant. Units: erg / K. */
#define BOLTZMANNCGS 1.38064852e-16

/*!< Boltzmann constant. Units: eV / K. */
#define BOLTZMANN_EVK 8.61733105e-5

/*!< For dust processes involving G0, e.g. photoelectric heating, we attenuate
 * G0 by exp(- G0_gamma * Av). */
#define G0_GAMMA 2.77

/*!< Effective dust absorption and scattering cross section, to calculate Av
 * from N_Htot. Units: cm^2 / H atom. */
#define DUSTEFFSIZE 4.0e-22

/*!< Dust surface area for H2 formation. Units: cm^2 / H atom. */
#define DUST_CROSS_SECTION 1.0e-21

/*!< All reaction and cooling rates are limited to this maximum temperature. */
#define MAX_TEMPERATURE_FOR_RATES 2.0e9

/*!< Exclude species when their corresponding element abundance is below this
 * threshold. */
#define METALS_MINIMUM_THRESHOLD 1.0e-30

/*!< The maximum number of steps in the CVODE solver. */
#define MAXSTEPS 1e5

/**
 * Structure containing information
 * about each non-equilibrium species
 */
struct Species_Structure {
  int include_species; /*!< Flag indicating whether the species is included in
                          the network. */
  ChimesFloat element_abundance; /*!< Abundance of the corresponding element. */
  ChimesFloat
      creation_rate; /*!< Sum of the rates of all reactions that create this
                        species. */
  ChimesFloat destruction_rate; /*!< Sum of the rates of all reactions that
                                   destroy this species. */
};
