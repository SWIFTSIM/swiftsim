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

/*
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
__attribute__((always_inline)) inline void chimes_get_table_index(
    ChimesFloat *table, int ntable, ChimesFloat x, int *i, ChimesFloat *dx) {
  ChimesFloat denominator;

  denominator = (table[ntable - 1] - table[0]) / (ntable - 1.0f);

  if (x <= table[0]) {
    *i = 0;
    *dx = 0.0f;
  } else if (x >= table[ntable - 1]) {
    *i = ntable - 2;
    *dx = 1.0f;
  } else {
    *i = (int)floor((x - table[0]) / denominator);
    *dx = (x - table[*i]) / denominator;
  }
}

/*
 * @brief Flatten 2-d array indices.
 *
 * Converts 2-d array indices into an index in
 * the flattened 1-d array.
 *
 * @param x The position in the first dimension.
 * @param y The position in the second dimension.
 * @param Ny The size of the second dimension.
 */
__attribute__((always_inline)) inline int chimes_flatten_index_2d(
    const int x, const int y, const int Ny) {
  return x * Ny + y;
}

/*
 * @brief Flatten 3-d array indices.
 *
 * Converts 3-d array indices into an index in
 * the flattened 1-d array.
 *
 * @param x The position in the first dimension.
 * @param y The position in the second dimension.
 * @param z The position in the third dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 */
__attribute__((always_inline)) inline int chimes_flatten_index_3d(
    const int x, const int y, const int z, const int Ny, const int Nz) {
  return x * Ny * Nz + y * Nz + z;
}

/*
 * @brief Flatten 4-d array indices.
 *
 * Converts 4-d array indices into an index in
 * the flattened 1-d array.
 *
 * @param x The position in the first dimension.
 * @param y The position in the second dimension.
 * @param z The position in the third dimension.
 * @param v The position in the fourth dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 * @param Nv The size of the fourth dimension.
 */
__attribute__((always_inline)) inline int chimes_flatten_index_4d(
    const int x, const int y, const int z, const int v, const int Ny,
    const int Nz, const int Nv) {
  return x * Ny * Nz * Nv + y * Nz * Nv + z * Nv + v;
}

/*
 * @brief Flatten 5-d array indices.
 *
 * Converts 5-d array indices into an index in
 * the flattened 1-d array.
 *
 * @param x The position in the first dimension.
 * @param y The position in the second dimension.
 * @param z The position in the third dimension.
 * @param v The position in the fourth dimension.
 * @param w The position in the fifth dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 * @param Nv The size of the fourth dimension.
 * @param Nw The size of the fifth dimension.
 */
__attribute__((always_inline)) inline int chimes_flatten_index_5d(
    const int x, const int y, const int z, const int v, const int w,
    const int Ny, const int Nz, const int Nv, const int Nw) {
  return x * Ny * Nz * Nv * Nw + y * Nz * Nv * Nw + z * Nv * Nw + v * Nw + w;
}

/*
 * @brief Perform a linear interpolation.
 *
 * Performs a linear interpolation on a 1-d table,
 * based on the index and displacement from the
 * #chimes_get_table_index() routine.
 *
 * @param table The 1-dimensional table.
 * @param x The position in the first dimension.
 * @param dx The displacement in the first dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_1d(
    const ChimesFloat *table, const int x, const ChimesFloat dx) {
  return (1.0f - dx) * table[x] + dx * table[x + 1];
}

/*
 * @brief Perform a bi-linear interpolation.
 *
 * Performs a bi-linear interpolation on a 2-d table,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 2-dimensional table.
 * @param x The position in the first dimension.
 * @param y The position in the second dimension.
 * @param dx The displacement in the first dimension.
 * @param dy The displacement in the second dimension.
 * @param Ny The size of the second dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_2d(
    const ChimesFloat *table, const int x, const int y, const ChimesFloat dx,
    const ChimesFloat dy, const int Ny) {
  ChimesFloat output;

  const ChimesFloat dx_m = 1.0f - dx;
  const ChimesFloat dy_m = 1.0f - dy;

  output = dx_m * dy_m * table[chimes_flatten_index_2d(x, y, Ny)];
  output += dx_m * dy * table[chimes_flatten_index_2d(x, y + 1, Ny)];
  output += dx * dy_m * table[chimes_flatten_index_2d(x + 1, y, Ny)];
  output += dx * dy * table[chimes_flatten_index_2d(x + 1, y + 1, Ny)];

  return output;
}

/*
 * @brief Perform a linear interpolation on a 2-d table at fixed x.
 *
 * Performs a linear interpolation on a 2-d table with a fixed x,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 2-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The position in the second dimension.
 * @param dy The displacement in the second dimension.
 * @param Ny The size of the second dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_2d_fix_x(
    const ChimesFloat *table, const int x, const int y, const ChimesFloat dy,
    const int Ny) {
  ChimesFloat output;

  const ChimesFloat dy_m = (1.0f - dy);

  output = dy_m * table[chimes_flatten_index_2d(x, y, Ny)];
  output += dy * table[chimes_flatten_index_2d(x, y + 1, Ny)];

  return output;
}

/*
 * @brief Perform a bi-linear interpolation on a 3-d table at fixed x.
 *
 * Performs a bi-linear interpolation on a 3-d table at fixed x,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 3-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The position in the second dimension.
 * @param z The position in the third dimension.
 * @param dy The displacement in the second dimension.
 * @param dz The displacement in the third dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_3d_fix_x(
    const ChimesFloat *table, const int x, const int y, const int z,
    const ChimesFloat dy, const ChimesFloat dz, const int Ny, const int Nz) {
  ChimesFloat output;

  const ChimesFloat dy_m = (1.0f - dy);
  const ChimesFloat dz_m = (1.0f - dz);

  output = dy_m * dz_m * table[chimes_flatten_index_3d(x, y, z, Ny, Nz)];
  output += dy_m * dz * table[chimes_flatten_index_3d(x, y, z + 1, Ny, Nz)];
  output += dy * dz_m * table[chimes_flatten_index_3d(x, y + 1, z, Ny, Nz)];
  output += dy * dz * table[chimes_flatten_index_3d(x, y + 1, z + 1, Ny, Nz)];

  return output;
}

/*
 * @brief Perform a linear interpolation on a 3-d table at fixed xy.
 *
 * Performs a linear interpolation on a 3-d table at fixed xy,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 3-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The fixed position in the second dimension.
 * @param z The position in the third dimension.
 * @param dz The displacement in the third dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_3d_fix_xy(
    const ChimesFloat *table, const int x, const int y, const int z,
    const ChimesFloat dz, const int Ny, const int Nz) {
  ChimesFloat output;

  const ChimesFloat dz_m = (1.0f - dz);

  output = dz_m * table[chimes_flatten_index_3d(x, y, z, Ny, Nz)];
  output += dz * table[chimes_flatten_index_3d(x, y, z + 1, Ny, Nz)];

  return output;
}

/*
 * @brief Perform a tri-linear interpolation on a 4-d table at fixed x.
 *
 * Performs a tri-linear interpolation on a 4-d table at fixed x,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 4-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The position in the second dimension.
 * @param z The position in the third dimension.
 * @param v The position in the fourth dimension.
 * @param dy The displacement in the second dimension.
 * @param dz The displacement in the third dimension.
 * @param dv The displacement in the fourth dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 * @param Nv The size of the fourth dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_4d_fix_x(
    const ChimesFloat *table, const int x, const int y, const int z,
    const int v, const ChimesFloat dy, const ChimesFloat dz,
    const ChimesFloat dv, const int Ny, const int Nz, const int Nv) {
  ChimesFloat output;

  const ChimesFloat dy_m = (1.0f - dy);
  const ChimesFloat dz_m = (1.0f - dz);
  const ChimesFloat dv_m = (1.0f - dv);

  output = dy_m * dz_m * dv_m *
           table[chimes_flatten_index_4d(x, y, z, v, Ny, Nz, Nv)];
  output += dy_m * dz_m * dv *
            table[chimes_flatten_index_4d(x, y, z, v + 1, Ny, Nz, Nv)];
  output += dy_m * dz * dv_m *
            table[chimes_flatten_index_4d(x, y, z + 1, v, Ny, Nz, Nv)];
  output += dy * dz_m * dv_m *
            table[chimes_flatten_index_4d(x, y + 1, z, v, Ny, Nz, Nv)];
  output += dy * dz * dv_m *
            table[chimes_flatten_index_4d(x, y + 1, z + 1, v, Ny, Nz, Nv)];
  output += dy * dz_m * dv *
            table[chimes_flatten_index_4d(x, y + 1, z, v + 1, Ny, Nz, Nv)];
  output += dy_m * dz * dv *
            table[chimes_flatten_index_4d(x, y, z + 1, v + 1, Ny, Nz, Nv)];
  output += dy * dz * dv *
            table[chimes_flatten_index_4d(x, y + 1, z + 1, v + 1, Ny, Nz, Nv)];

  return output;
}

/*
 * @brief Perform a linear interpolation on a 4-d table at fixed xyz.
 *
 * Performs a linear interpolation on a 4-d table at fixed xyz,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 4-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The fixed position in the second dimension.
 * @param z The fixed position in the third dimension.
 * @param v The position in the fourth dimension.
 * @param dv The displacement in the fourth dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 * @param Nv The size of the fourth dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_4d_fix_xyz(
    const ChimesFloat *table, const int x, const int y, const int z,
    const int v, const ChimesFloat dv, const int Ny, const int Nz,
    const int Nv) {
  ChimesFloat output;

  const ChimesFloat dv_m = (1.0f - dv);

  output = dv_m * table[chimes_flatten_index_4d(x, y, z, v, Ny, Nz, Nv)];
  output += dv * table[chimes_flatten_index_4d(x, y, z, v + 1, Ny, Nz, Nv)];

  return output;
}

/*
 * @brief Perform a quadri-linear interpolation on a 5-d table at fixed x.
 *
 * Performs a quadri-linear interpolation on a 5-d table at fixed x,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 5-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The position in the second dimension.
 * @param z The position in the third dimension.
 * @param v The position in the fourth dimension.
 * @param w The position in the fifth dimension.
 * @param dy The displacement in the second dimension.
 * @param dz The displacement in the third dimension.
 * @param dv The displacement in the fourth dimension.
 * @param dw The displacement in the fifth dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 * @param Nv The size of the fourth dimension.
 * @param Nw The size of the fifth dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_5d_fix_x(
    const ChimesFloat *table, const int x, const int y, const int z,
    const int v, const int w, const ChimesFloat dy, const ChimesFloat dz,
    const ChimesFloat dv, const ChimesFloat dw, const int Ny, const int Nz,
    const int Nv, const int Nw) {
  ChimesFloat output;

  const ChimesFloat dy_m = (1.0f - dy);
  const ChimesFloat dz_m = (1.0f - dz);
  const ChimesFloat dv_m = (1.0f - dv);
  const ChimesFloat dw_m = (1.0f - dw);

  output = dy_m * dz_m * dv_m * dw_m *
           table[chimes_flatten_index_5d(x, y, z, v, w, Ny, Nz, Nv, Nw)];
  output += dy_m * dz_m * dv * dw_m *
            table[chimes_flatten_index_5d(x, y, z, v + 1, w, Ny, Nz, Nv, Nw)];
  output += dy_m * dz * dv_m * dw_m *
            table[chimes_flatten_index_5d(x, y, z + 1, v, w, Ny, Nz, Nv, Nw)];
  output += dy * dz_m * dv_m * dw_m *
            table[chimes_flatten_index_5d(x, y + 1, z, v, w, Ny, Nz, Nv, Nw)];
  output +=
      dy * dz * dv_m * dw_m *
      table[chimes_flatten_index_5d(x, y + 1, z + 1, v, w, Ny, Nz, Nv, Nw)];
  output +=
      dy * dz_m * dv * dw_m *
      table[chimes_flatten_index_5d(x, y + 1, z, v + 1, w, Ny, Nz, Nv, Nw)];
  output +=
      dy_m * dz * dv * dw_m *
      table[chimes_flatten_index_5d(x, y, z + 1, v + 1, w, Ny, Nz, Nv, Nw)];
  output +=
      dy * dz * dv * dw_m *
      table[chimes_flatten_index_5d(x, y + 1, z + 1, v + 1, w, Ny, Nz, Nv, Nw)];
  output += dy_m * dz_m * dv_m * dw *
            table[chimes_flatten_index_5d(x, y, z, v, w + 1, Ny, Nz, Nv, Nw)];
  output +=
      dy_m * dz_m * dv * dw *
      table[chimes_flatten_index_5d(x, y, z, v + 1, w + 1, Ny, Nz, Nv, Nw)];
  output +=
      dy_m * dz * dv_m * dw *
      table[chimes_flatten_index_5d(x, y, z + 1, v, w + 1, Ny, Nz, Nv, Nw)];
  output +=
      dy * dz_m * dv_m * dw *
      table[chimes_flatten_index_5d(x, y + 1, z, v, w + 1, Ny, Nz, Nv, Nw)];
  output +=
      dy * dz * dv_m * dw *
      table[chimes_flatten_index_5d(x, y + 1, z + 1, v, w + 1, Ny, Nz, Nv, Nw)];
  output +=
      dy * dz_m * dv * dw *
      table[chimes_flatten_index_5d(x, y + 1, z, v + 1, w + 1, Ny, Nz, Nv, Nw)];
  output +=
      dy_m * dz * dv * dw *
      table[chimes_flatten_index_5d(x, y, z + 1, v + 1, w + 1, Ny, Nz, Nv, Nw)];
  output += dy * dz * dv * dw *
            table[chimes_flatten_index_5d(x, y + 1, z + 1, v + 1, w + 1, Ny, Nz,
                                          Nv, Nw)];

  return output;
}

/*
 * @brief Perform a bi-linear interpolation on a 5-d table at fixed xyz.
 *
 * Performs a bi-linear interpolation on a 5-d table at fixed xyz,
 * based on the indices and displacements in each
 * dimension from the #chimes_get_table_index() routine.
 *
 * @param table The 5-dimensional table.
 * @param x The fixed position in the first dimension.
 * @param y The fixed position in the second dimension.
 * @param z The fixed position in the third dimension.
 * @param v The position in the fourth dimension.
 * @param w The position in the fifth dimension.
 * @param dv The displacement in the fourth dimension.
 * @param dw The displacement in the fifth dimension.
 * @param Ny The size of the second dimension.
 * @param Nz The size of the third dimension.
 * @param Nv The size of the fourth dimension.
 * @param Nw The size of the fifth dimension.
 */
__attribute__((always_inline)) inline ChimesFloat chimes_interpol_5d_fix_xyz(
    const ChimesFloat *table, const int x, const int y, const int z,
    const int v, const int w, const ChimesFloat dv, const ChimesFloat dw,
    const int Ny, const int Nz, const int Nv, const int Nw) {
  ChimesFloat output;

  const ChimesFloat dv_m = (1.0f - dv);
  const ChimesFloat dw_m = (1.0f - dw);

  output = dv_m * dw_m *
           table[chimes_flatten_index_5d(x, y, z, v, w, Ny, Nz, Nv, Nw)];
  output += dv_m * dw *
            table[chimes_flatten_index_5d(x, y, z, v, w + 1, Ny, Nz, Nv, Nw)];
  output += dv * dw_m *
            table[chimes_flatten_index_5d(x, y, z, v + 1, w, Ny, Nz, Nv, Nw)];
  output +=
      dv * dw *
      table[chimes_flatten_index_5d(x, y, z, v + 1, w + 1, Ny, Nz, Nv, Nw)];

  return output;
}
