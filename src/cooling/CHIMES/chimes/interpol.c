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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "chimes_proto.h"

/**
 * @brief Get table index.
 *
 * Gets the index of a given values within a 1-d table,
 * along with the displacement of that value between
 * the discrete table values. This routine assumes
 * that the table is evenly spaced.
 *
 * @param table The 1-dimensional table.
 * @param ntable The length of the table.
 * @param x The value that we are looking up in the table.
 * @param i Output index in the table.
 * @param dx Output displacement of the value between discrete table values.
 */
void chimes_get_table_index(ChimesFloat *table, int ntable, ChimesFloat x,
                            int *i, ChimesFloat *dx) {
  ChimesFloat denominator;

  denominator = (table[ntable - 1] - table[0]) / (ntable - 1.0);

  if (x <= table[0]) {
    *i = 0;
    *dx = 0.0;
  } else if (x >= table[ntable - 1]) {
    *i = ntable - 2;
    *dx = 1.0;
  } else {
    *i = (int)floor((x - table[0]) / denominator);
    *dx = (x - table[*i]) / denominator;
  }
}

/**
 * @brief Perform a linear interpolation.
 *
 * Performs a linear interpolation on a 1-d table,
 * based on the index and displacement from the
 * #chimes_get_table_index() routine.
 *
 * @param table The 1-dimensional table.
 * @param i The position in the table.
 * @param dx The displacement between positions i and i + 1.
 */
ChimesFloat chimes_interpol_1d(ChimesFloat *table, int i, ChimesFloat dx) {
  return (1 - dx) * table[i] + dx * table[i + 1];
}

/**
 * @brief Perform a bi-linear interpolation.
 *
 * Performs a bi-linear interpolation on a 2-d table,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 2-dimensional table.
 * @param i The position in the first dimension.
 * @param j The position in the second dimension.
 * @param dx The displacement in the first dimension.
 * @param dy The displacement in the second dimension.
 */
ChimesFloat chimes_interpol_2d(ChimesFloat **table, int i, int j,
                               ChimesFloat dx, ChimesFloat dy) {
  ChimesFloat output, dx_m, dy_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;

  output = dx_m * dy_m * table[i][j] + dx_m * dy * table[i][j + 1] +
           dx * dy_m * table[i + 1][j] + dx * dy * table[i + 1][j + 1];

  return output;
}

/**
 * @brief Perform a tri-linear interpolation.
 *
 * Performs a tri-linear interpolation on a 3-d table,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 3-dimensional table.
 * @param i The position in the first dimension.
 * @param j The position in the second dimension.
 * @param k The position in the third dimension.
 * @param dx The displacement in the first dimension.
 * @param dy The displacement in the second dimension.
 * @param dz The displacement in the third dimension.
 */
ChimesFloat chimes_interpol_3d(ChimesFloat ***table, int i, int j, int k,
                               ChimesFloat dx, ChimesFloat dy, ChimesFloat dz) {
  ChimesFloat output, dx_m, dy_m, dz_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;
  dz_m = 1.0 - dz;

  output = dx_m * dy_m * dz_m * table[i][j][k] +
           dx_m * dy_m * dz * table[i][j][k + 1] +
           dx_m * dy * dz_m * table[i][j + 1][k] +
           dx_m * dy * dz * table[i][j + 1][k + 1] +
           dx * dy_m * dz_m * table[i + 1][j][k] +
           dx * dy_m * dz * table[i + 1][j][k + 1] +
           dx * dy * dz_m * table[i + 1][j + 1][k] +
           dx * dy * dz * table[i + 1][j + 1][k + 1];

  return output;
}

/**
 * @brief Perform a quadri-linear interpolation.
 *
 * Performs a quadri-linear interpolation on a 4-d table,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 4-dimensional table.
 * @param i The position in the first dimension.
 * @param j The position in the second dimension.
 * @param k The position in the third dimension.
 * @param l The position in the fourth dimension.
 * @param dx The displacement in the first dimension.
 * @param dy The displacement in the second dimension.
 * @param dz The displacement in the third dimension.
 * @param dw The displacement in the fourth dimension.
 */
ChimesFloat chimes_interpol_4d(ChimesFloat ****table, int i, int j, int k,
                               int l, ChimesFloat dx, ChimesFloat dy,
                               ChimesFloat dz, ChimesFloat dw) {
  ChimesFloat output, dx_m, dy_m, dz_m, dw_m;

  dx_m = 1.0 - dx;
  dy_m = 1.0 - dy;
  dz_m = 1.0 - dz;
  dw_m = 1.0 - dw;

  output = dx_m * dy_m * dz_m * dw_m * table[i][j][k][l] +
           dx_m * dy_m * dz_m * dw * table[i][j][k][l + 1] +
           dx_m * dy_m * dz * dw_m * table[i][j][k + 1][l] +
           dx_m * dy_m * dz * dw * table[i][j][k + 1][l + 1] +
           dx_m * dy * dz_m * dw_m * table[i][j + 1][k][l] +
           dx_m * dy * dz_m * dw * table[i][j + 1][k][l + 1] +
           dx_m * dy * dz * dw_m * table[i][j + 1][k + 1][l] +
           dx_m * dy * dz * dw * table[i][j + 1][k + 1][l + 1] +
           dx * dy_m * dz_m * dw_m * table[i + 1][j][k][l] +
           dx * dy_m * dz_m * dw * table[i + 1][j][k][l + 1] +
           dx * dy_m * dz * dw_m * table[i + 1][j][k + 1][l] +
           dx * dy_m * dz * dw * table[i + 1][j][k + 1][l + 1] +
           dx * dy * dz_m * dw_m * table[i + 1][j + 1][k][l] +
           dx * dy * dz_m * dw * table[i + 1][j + 1][k][l + 1] +
           dx * dy * dz * dw_m * table[i + 1][j + 1][k + 1][l] +
           dx * dy * dz * dw * table[i + 1][j + 1][k + 1][l + 1];

  return output;
}
