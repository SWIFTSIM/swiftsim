/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_SESAME_EQUATION_OF_STATE_H
#define SWIFT_SESAME_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/sesame.h
 *
 * Contains the SESAME EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 *
 *              WORK IN PROGRESS!
 *
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "units.h"
#include "utilities.h"

// SESAME parameters
struct SESAME_params {
  float *table_log_rho;
  float *table_log_u_rho_T;
  float *table_P_rho_T;
  float *table_c_rho_T;
  int num_rho, num_T;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material (cgs units)
INLINE static void set_SESAME_basalt(struct SESAME_params *mat,
                                     enum eos_planetary_material_id mat_id) {
  // SESAME 7530
  mat->mat_id = mat_id;
}
INLINE static void set_SESAME_water(struct SESAME_params *mat,
                                    enum eos_planetary_material_id mat_id) {
  // SESAME 7154
  mat->mat_id = mat_id;
}

// Read the tables from file
INLINE static void load_table_SESAME(struct SESAME_params *mat,
                                     char *table_file) {

  // Load table contents from file
  FILE *f = fopen(table_file, "r");
  int c;

  // Ignore header lines
  char buffer[64];
  for (int i = 0; i < 4; i++) fgets(buffer, 64, f);

  // Table properties
  c = fscanf(f, "%d %d", &mat->num_rho, &mat->num_T);
  if (c != 2) {
    error("Failed to read EOS table");
  }

  // Allocate table memory
  mat->table_log_rho = (float *)malloc(mat->num_rho * sizeof(float));
  mat->table_log_u_rho_T = (float *)malloc(mat->num_rho * mat->num_T *
                                           sizeof(float));
  mat->table_P_rho_T = (float *)malloc(mat->num_rho * mat->num_T *
                                       sizeof(float));
  mat->table_c_rho_T = (float *)malloc(mat->num_rho * mat->num_T *
                                       sizeof(float));

  // Densities (not log yet)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    c = fscanf(f, "%f", &mat->table_log_rho[i_rho]);
    if (c != 1) {
      error("Failed to read EOS table");
    }
  }

  // Sp. int. energies (not log yet), pressures, and sound speeds
  for (int i_T = 0; i_T < mat->num_T; i_T++) {
    for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
      c = fscanf(f, "%f %f %f",
                 &mat->table_log_u_rho_T[i_rho*mat->num_T + i_T],
                 &mat->table_P_rho_T[i_rho*mat->num_T + i_T],
                 &mat->table_c_rho_T[i_rho*mat->num_T + i_T]);
      if (c != 3) {
        error("Failed to read EOS table");
      }
    }
  }

  fclose(f);
}

// Misc. modifications
INLINE static void prepare_table_SESAME(struct SESAME_params *mat,
                                        const struct unit_system *us) {

  // Convert densities to log(density)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    mat->table_log_rho[i_rho] = logf(mat->table_log_rho[i_rho]);
  }

  // Convert sp. int. energies to log(sp. int. energy)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = 0; i_T < mat->num_T; i_T++) {
      // If not positive then set very small for the log
      if (mat->table_log_u_rho_T[i_rho*mat->num_T + i_T] <= 0) {
        mat->table_log_u_rho_T[i_rho*mat->num_T + i_T] = 1.f;
      }

      mat->table_log_u_rho_T[i_rho*mat->num_T + i_T] =
        logf(mat->table_log_u_rho_T[i_rho*mat->num_T + i_T]);
    }
  }

  // Enforce that the 1D arrays of u (at each rho) are monotonic
  // This is necessary because, for some high-density u slices at very low T,
  // u decreases (very slightly) with T, which makes the interpolation fail
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = mat->num_T-1; i_T > 0; i_T--) {

      // If the one-lower-T u is greater than this u
      if (mat->table_log_u_rho_T[i_rho*mat->num_T + i_T] <
          mat->table_log_u_rho_T[i_rho*mat->num_T + i_T - 1]) {

        // Replace it and all elements below it with that value
        for (int j_u = 0; j_u < i_T; j_u++) {
          mat->table_log_u_rho_T[i_rho*mat->num_T + j_u] =
              mat->table_log_u_rho_T[i_rho*mat->num_T + i_T];
        }
        break;
      }
    }
  }
}

// Convert to internal units
INLINE static void convert_units_SESAME(struct SESAME_params *mat,
                                        const struct unit_system *us) {

  const float kg_m3_to_g_cm3 = 1e-3f;   // Convert kg/m^3 to g/cm^3
  const float Pa_to_Ba = 1e1f;          // Convert Pascals to Barye
  const float J_kg_to_erg_g = 1e4f;     // Convert J/kg to erg/g
  const float m_s_to_cm_s = 1e2f;       // Convert m/s to cm/s

  // All table values in SI
  // Densities (log)
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
      mat->table_log_rho[i_rho] +=
          logf(kg_m3_to_g_cm3 /
               units_cgs_conversion_factor(us, UNIT_CONV_DENSITY));
  }

  // Sp. Int. Energies (log), pressures, and sound speeds
  for (int i_rho = 0; i_rho < mat->num_rho; i_rho++) {
    for (int i_T = 0; i_T < mat->num_T; i_T++) {
      mat->table_log_u_rho_T[i_rho*mat->num_T + i_T] +=
          logf(J_kg_to_erg_g /
               units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS));
      mat->table_P_rho_T[i_rho*mat->num_T + i_T] *=
          Pa_to_Ba / units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
      mat->table_c_rho_T[i_rho*mat->num_T + i_T] *=
          m_s_to_cm_s / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
    }
  }
}

// gas_internal_energy_from_entropy
INLINE static float SESAME_internal_energy_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_pressure_from_entropy
INLINE static float SESAME_pressure_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_pressure
INLINE static float SESAME_entropy_from_pressure(
    float density, float pressure, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_entropy
INLINE static float SESAME_soundspeed_from_entropy(
    float density, float entropy, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_internal_energy
INLINE static float SESAME_entropy_from_internal_energy(
    float density, float u, const struct SESAME_params *mat) {

  return 0;
}

// gas_pressure_from_internal_energy
INLINE static float SESAME_pressure_from_internal_energy(
    float density, float u, const struct SESAME_params *mat) {

  float P;

  if (u <= 0.f) {
    return 0.f;
  }

  int rho_idx, u_idx_1, u_idx_2;
  float intp_rho, intp_u_1, intp_u_2;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u)) to find P(rho, u)
  // Density index
  find_value_in_monotonic_array(log_rho, mat->table_log_rho, mat->num_rho,
                                &rho_idx);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  find_value_in_monotonic_array(log_u,
                                mat->table_log_u_rho_T + rho_idx*mat->num_T,
                                mat->num_T, &u_idx_1);
  find_value_in_monotonic_array(log_u,
                                mat->table_log_u_rho_T + (rho_idx+1)*mat->num_T,
                                mat->num_T, &u_idx_2);

  intp_rho = (log_rho - mat->table_log_rho[rho_idx]) /
             (mat->table_log_rho[rho_idx+1] - mat->table_log_rho[rho_idx]);
  intp_u_1 = (log_u - mat->table_log_u_rho_T[rho_idx*mat->num_T + u_idx_1]) /
             (mat->table_log_u_rho_T[rho_idx*mat->num_T + (u_idx_1+1)] -
              mat->table_log_u_rho_T[rho_idx*mat->num_T + u_idx_1]);
  intp_u_2 = (log_u - mat->table_log_u_rho_T[(rho_idx+1)*mat->num_T + u_idx_2]) /
             (mat->table_log_u_rho_T[(rho_idx+1)*mat->num_T + (u_idx_2+1)] -
              mat->table_log_u_rho_T[(rho_idx+1)*mat->num_T + u_idx_2]);

  // Return zero if outside the table
  if (rho_idx < 0) {                            // Too-low rho
    P = 0.f;
    if ((u_idx_1 < 0) || (u_idx_2 < 0)) {       // and too-low u
      P = 0.f;
    }
  } else if ((u_idx_1 < 0) || (u_idx_2 < 0)) {  // Too-low u
    P = 0.f;
  }
  else if (rho_idx >= mat->num_rho - 1) {       // Too-high rho
    if ((u_idx_1 >= mat->num_T - 1) ||
        (u_idx_2 >= mat->num_T - 1)) {          // and too-high u
      P = 0.f;
    } else {
      P = 0.f;
    }
  } else if ((u_idx_1 >= mat->num_T - 1) ||
             (u_idx_2 >= mat->num_T - 1)) {     // Too-high u
    P = 0.f;
  }
  // Normal interpolation within the table
  else {
    P = (1.f - intp_rho) * (
            (1.f - intp_u_1) *
                mat->table_P_rho_T[rho_idx*mat->num_T + u_idx_1] +
            intp_u_1 *
                mat->table_P_rho_T[rho_idx*mat->num_T + u_idx_1 + 1]) +
        intp_rho * (
            (1.f - intp_u_2) *
                mat->table_P_rho_T[(rho_idx + 1)*mat->num_T + u_idx_2] +
            intp_u_2 *
                mat->table_P_rho_T[(rho_idx + 1)*mat->num_T + u_idx_2 + 1]);
  }

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float SESAME_internal_energy_from_pressure(
    float density, float P, const struct SESAME_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_internal_energy
INLINE static float SESAME_soundspeed_from_internal_energy(
    float density, float u, const struct SESAME_params *mat) {

  float c;

  if (u <= 0.f) {
    return 0.f;
  }

  int rho_idx, u_idx_1, u_idx_2;
  float intp_rho, intp_u_1, intp_u_2;
  const float log_rho = logf(density);
  const float log_u = logf(u);

  // 2D interpolation (bilinear with log(rho), log(u)) to find c(rho, u)
  // Density index
  find_value_in_monotonic_array(log_rho, mat->table_log_rho, mat->num_rho,
                                &rho_idx);

  // Sp. int. energy at this and the next density (in relevant slice of u array)
  find_value_in_monotonic_array(log_u,
                                mat->table_log_u_rho_T + rho_idx*mat->num_T,
                                mat->num_T, &u_idx_1);
  find_value_in_monotonic_array(log_u,
                                mat->table_log_u_rho_T + (rho_idx+1)*mat->num_T,
                                mat->num_T, &u_idx_2);

  intp_rho = (log_rho - mat->table_log_rho[rho_idx]) /
             (mat->table_log_rho[rho_idx+1] - mat->table_log_rho[rho_idx]);
  intp_u_1 = (log_u - mat->table_log_u_rho_T[rho_idx*mat->num_T + u_idx_1]) /
             (mat->table_log_u_rho_T[rho_idx*mat->num_T + (u_idx_1+1)] -
              mat->table_log_u_rho_T[rho_idx*mat->num_T + u_idx_1]);
  intp_u_2 = (log_u - mat->table_log_u_rho_T[(rho_idx+1)*mat->num_T + u_idx_2]) /
             (mat->table_log_u_rho_T[(rho_idx+1)*mat->num_T + (u_idx_2+1)] -
              mat->table_log_u_rho_T[(rho_idx+1)*mat->num_T + u_idx_2]);

  // Return zero if outside the table
  if (rho_idx < 0) {                            // Too-low rho
    c = 0.f;
    if ((u_idx_1 < 0) || (u_idx_2 < 0)) {       // and too-low u
      c = 0.f;
    }
  } else if ((u_idx_1 < 0) || (u_idx_2 < 0)) {  // Too-low u
    c = 0.f;
  }
  else if (rho_idx >= mat->num_rho - 1) {       // Too-high rho
    if ((u_idx_1 >= mat->num_T - 1) ||
        (u_idx_2 >= mat->num_T - 1)) {          // and too-high u
      c = 0.f;
    } else {
      c = 0.f;
    }
  } else if ((u_idx_1 >= mat->num_T - 1) ||
             (u_idx_2 >= mat->num_T - 1)) {     // Too-high u
    c = 0.f;
  }
  // Normal interpolation within the table
  else {
    c = (1.f - intp_rho) * (
            (1.f - intp_u_1) *
                mat->table_c_rho_T[rho_idx*mat->num_T + u_idx_1] +
            intp_u_1 *
                mat->table_c_rho_T[rho_idx*mat->num_T + u_idx_1 + 1]) +
        intp_rho * (
            (1.f - intp_u_2) *
                mat->table_c_rho_T[(rho_idx + 1)*mat->num_T + u_idx_2] +
            intp_u_2 *
                mat->table_c_rho_T[(rho_idx + 1)*mat->num_T + u_idx_2 + 1]);
  }

  return c;
}

// gas_soundspeed_from_pressure
INLINE static float SESAME_soundspeed_from_pressure(
    float density, float P, const struct SESAME_params *mat) {

  float c;

  int rho_idx, P_idx_1, P_idx_2;
  float intp_rho, intp_P_1, intp_P_2;
  const float log_rho = logf(density);

  // 2D interpolation (bilinear with log(rho), P) to find c(rho, P)
  // Density index
  find_value_in_monotonic_array(log_rho, mat->table_log_rho, mat->num_rho,
                                &rho_idx);

  // Pressure at this and the next density (in relevant slice of P array)
  find_value_in_monotonic_array(P,
                                mat->table_P_rho_T + rho_idx*mat->num_T,
                                mat->num_T, &P_idx_1);
  find_value_in_monotonic_array(P,
                                mat->table_P_rho_T + (rho_idx+1)*mat->num_T,
                                mat->num_T, &P_idx_2);

  intp_rho = (log_rho - mat->table_log_rho[rho_idx]) /
             (mat->table_log_rho[rho_idx+1] - mat->table_log_rho[rho_idx]);
  intp_P_1 = (P - mat->table_P_rho_T[rho_idx*mat->num_T + P_idx_1]) /
             (mat->table_P_rho_T[rho_idx*mat->num_T + (P_idx_1+1)] -
              mat->table_P_rho_T[rho_idx*mat->num_T + P_idx_1]);
  intp_P_2 = (P - mat->table_P_rho_T[(rho_idx+1)*mat->num_T + P_idx_2]) /
             (mat->table_P_rho_T[(rho_idx+1)*mat->num_T + (P_idx_2+1)] -
              mat->table_P_rho_T[(rho_idx+1)*mat->num_T + P_idx_2]);

  // Return zero if outside the table
  if (rho_idx < 0) {                            // Too-low rho
    c = 0.f;
    if ((P_idx_1 < 0) || (P_idx_2 < 0)) {       // and too-low P
      c = 0.f;
    }
  } else if ((P_idx_1 < 0) || (P_idx_2 < 0)) {  // Too-low P
    c = 0.f;
  }
  else if (rho_idx >= mat->num_rho - 1) {       // Too-high rho
    if ((P_idx_1 >= mat->num_T - 1) ||
        (P_idx_2 >= mat->num_T - 1)) {          // and too-high P
      c = 0.f;
    } else {
      c = 0.f;
    }
  } else if ((P_idx_1 >= mat->num_T - 1) ||
             (P_idx_2 >= mat->num_T - 1)) {     // Too-high P
    c = 0.f;
  }
  // Normal interpolation within the table
  else {
    c = (1.f - intp_rho) * (
            (1.f - intp_P_1) *
                mat->table_c_rho_T[rho_idx*mat->num_T + P_idx_1] +
            intp_P_1 *
                mat->table_c_rho_T[rho_idx*mat->num_T + P_idx_1 + 1]) +
        intp_rho * (
            (1.f - intp_P_2) *
                mat->table_c_rho_T[(rho_idx + 1)*mat->num_T + P_idx_2] +
            intp_P_2 *
                mat->table_c_rho_T[(rho_idx + 1)*mat->num_T + P_idx_2 + 1]);
  }

  return c;
}

#endif /* SWIFT_SESAME_EQUATION_OF_STATE_H */
