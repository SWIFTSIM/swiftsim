/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_INTERPOLATION_H
#define SWIFT_GEAR_INTERPOLATION_H

struct interpolation_1d {
  /* Data to interpolate */
  float *data;

  /* Minimal x */
  float xmin;

  /* Step size between x points */
  float dx;

  /* Number of element in the data */
  int N;
};

/**
 * @brief Initialize the #interpolation_1d.
 *
 * @params interp The #interpolation_1d.
 * @params xmin Minimal value of x.
 * @params xmax Maximal value of x.
 * @params data The data to interpolate (y).
 * @params N The number of element in data.
 * @params integrate Do you want to interpolate the integral?
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_init(
    struct interpolation_1d *interp, float xmin, float xmax,
    const float *data, int N, int integrate) {

  /* Save the variables */
  interp->N = N;
  interp->xmin = xmin;
  interp->dx = (xmax - xmin / N);

  /* Copy the data */
  interp->data = malloc(sizeof(float) * N);
  if (interp->data == NULL)
    error("Failed to allocate memory for the interpolation");
  memcpy(interp->data, data, sizeof(float) * N);

  /* If we do not need to integrate the data, leave */
  if (!integrate)
    return;

  /* Integrate the yields with the trapezoid rule */
  float *tmp = (float*) malloc(sizeof(float) * N);
  if (tmp == NULL)
    error("Failed to allocate temporary array");

  tmp[0] = 0.;
  for(int i = 1; i < N; i++) {
    tmp[i] = tmp[i-1] + 0.5 * (interp->data[i] + interp->data[i-1]);
  }

  /* Copy the data back to the correct array */
  memcpy(interp->data, tmp, sizeof(float) * N);

  /* Cleanup */
  free(tmp);
}

/**
 * @brief Change the units of the #interpolation_1d.
 *
 * @params interp The #interpolation_1d.
 * @params x_units The new units for x in unit of the previous ones.
 * @params y_units The new units for y in unit of the previous ones.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_change_units(
    struct interpolation_1d *interp, float x_units, float y_units) {

  interp->xmin *= x_units;
  interp->dx *= x_units;

  /* do the data */
  for(int j = 0; j < interp->N; j++) {
    interp->data[j] *= y_units;
  }
  
}


/**
 * @brief Interpolate the data.
 *
 * @params interp The #interpolation_1d.
 * @params x The x value where to interpolate.
 *
 * @return The interpolated value y.
 */
__attribute__((always_inline)) static INLINE float interpolate_1d(
    const struct interpolation_1d *interp, float x) {

  /* Find indice */
  const float i = (x - interp->xmin) / interp->dx;
  const int idx = i;
  const float dx = i - idx;

#ifdef SWIFT_DEBUG_CHECKS
  if (i >= interp->N || i < 0) {
    error("Cannot extrapolate (i=%g; N=%i).",
	  i, interp->N);
  }
#endif

  /* interpolate */
  return interp->data[idx] * (1. - dx) +
    interp->data[idx+1] * dx;
  
}

/**
 * @brief Cleanup the #interpolation_1d structure.
 *
 * @params interp The #interpolation_1d.
 */
__attribute__((always_inline)) static INLINE void interpolate_1d_free(
    struct interpolation_1d *interp) {

  /* Free the allocated memory */
  free(interp->data);
  interp->data = NULL;
}

#endif // SWIFT_GEAR_INTERPOLATION_H
