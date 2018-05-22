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
#ifndef SWIFT_TILLOTSON_EQUATION_OF_STATE_H
#define SWIFT_TILLOTSON_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/tillotson.h
 *
 * Contains the Tillotson EOS functions for
 * equation_of_state/planetary/equation_of_state.h
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

// Tillotson parameters
struct Til_params {
  float rho_0, a, b, A, B, u_0, u_iv, u_cv, alpha, beta, eta_min, P_min;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material (cgs units)
INLINE static void set_Til_iron(struct Til_params *mat,
                                enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 7.800f;
  mat->a = 0.5f;
  mat->b = 1.5f;
  mat->A = 1.28e12f;
  mat->B = 1.05e12f;
  mat->u_0 = 9.5e10f;
  mat->u_iv = 2.4e10f;
  mat->u_cv = 8.67e10f;
  mat->alpha = 5.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.0f;
  mat->P_min = 0.0f;
}
INLINE static void set_Til_granite(struct Til_params *mat,
                                   enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 2.680f;
  mat->a = 0.5f;
  mat->b = 1.3f;
  mat->A = 1.8e11f;
  mat->B = 1.8e11f;
  mat->u_0 = 1.6e11f;
  mat->u_iv = 3.5e10f;
  mat->u_cv = 1.8e11f;
  mat->alpha = 5.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.0f;
  mat->P_min = 0.0f;
}
INLINE static void set_Til_water(struct Til_params *mat,
                                 enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 0.998f;
  mat->a = 0.7f;
  mat->b = 0.15f;
  mat->A = 2.18e10f;
  mat->B = 1.325e11f;
  mat->u_0 = 7.0e10f;
  mat->u_iv = 4.19e9f;
  mat->u_cv = 2.69e10f;
  mat->alpha = 10.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.915f;
  mat->P_min = 0.0f;
}

// Convert from cgs to internal units
INLINE static void convert_units_Til(struct Til_params *mat,
                                     const struct unit_system *us) {

  mat->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  mat->A /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->B /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->u_0 /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->u_iv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->u_cv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->P_min /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
}

// gas_internal_energy_from_entropy
INLINE static float Til_internal_energy_from_entropy(
    float density, float entropy, const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_pressure_from_entropy
INLINE static float Til_pressure_from_entropy(float density, float entropy,
                                              const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_pressure
INLINE static float Til_entropy_from_pressure(float density, float pressure,
                                              const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_entropy
INLINE static float Til_soundspeed_from_entropy(float density, float entropy,
                                                const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_entropy_from_internal_energy
INLINE static float Til_entropy_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  return 0;
}

// gas_pressure_from_internal_energy
INLINE static float Til_pressure_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  const float eta = density / mat->rho_0;
  const float eta_sq = eta*eta;
  const float mu = eta - 1.f;
  const float nu = 1.f / eta - 1.f;
  const float w = u / (mat->u_0 * eta_sq) + 1.f;
  const float w_inv = 1.f / w_inv;
  float P_c, P_e, P;

  // Condensed or cold
  if (eta < mat->eta_min) {
    P_c = 0.f;
  } else {
    P_c = (mat->a + mat->b*w_inv) * density * u + mat->A * mu + mat->B * mu*mu;
  }
  // Expanded and hot
  P_e = mat->a * density * u +
        (mat->b * density * u * w_inv +
         mat->A * mu * expf(-mat->beta * nu)) * expf(-mat->alpha * nu*nu);

  // Condensed or cold state
  if ((1.f < eta) || (u < mat->u_iv)) {
    P = P_c;
  }
  // Expanded and hot state
  else if ((eta < 1.f) && (mat->u_cv < u)) {
    P = P_e;
  }
  // Hybrid state
  else {
    P = ((u - mat->u_iv) * P_e + (mat->u_cv - u) * P_c) /
        (mat->u_cv - mat->u_iv);
  }

  // Minimum pressure
  if (P < mat->P_min) {
    P = mat->P_min;
  }

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float Til_internal_energy_from_pressure(
    float density, float P, const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0;
}

// gas_soundspeed_from_internal_energy
INLINE static float Til_soundspeed_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  const float eta = density / mat->rho_0;
  const float rho_inv = 1.f / density;
  const float eta_sq = eta*eta;
  const float mu = eta - 1.f;
  const float nu = 1.f / eta - 1.f;
  const float w = u / (mat->u_0 * eta_sq) + 1.f;
  const float w_inv = 1.f / w_inv;
  const float w_inv_sq = w_inv*w_inv;
  float P_c, P_e, P, c_sq_c, c_sq_e, c_sq;

  // Condensed or cold
  if (eta < mat->eta_min) {
      P_c = 0.f;
  }
  else {
    P_c = (mat->a + mat->b*w_inv) * density * u + mat->A * mu + mat->B * mu*mu;
  }
  c_sq_c = (1.f - mat->a - mat->b * w_inv) * P_c * rho_inv
           + mat->b * (w - 1.f) * w_inv_sq * (2.f*u + P_c * rho_inv)
           + (mat->A + mat->B * (eta_sq - 1.f)) * rho_inv;

  c_sq_c = max(c_sq_c, mat->A / mat->rho_0);

  // Expanded and hot
  P_e = mat->a * density * u +
        (mat->b * density * u * w_inv +
         mat->A * mu * expf(-mat->beta * nu)) * expf(-mat->alpha * nu*nu);

  c_sq_e = 1;///###wilo

  // Condensed or cold state
  if ((1.f < eta) || (u < mat->u_iv)) {
      c_sq = c_sq_c;
  }
  // Expanded and hot state
  else if ((eta < 1.f) && (mat->u_cv < u)) {
      c_sq = c_sq_e;
  }
  // Hybrid state
  else {
    c_sq = ((u - mat->u_iv)*c_sq_e + (mat->u_cv - u)*c_sq_c) /
          (mat->u_cv - mat->u_iv);

    c_sq = max(c_sq_c, mat->A / mat->rho_0);
  }

  return sqrtf(c_sq);
}

// gas_soundspeed_from_pressure
INLINE static float Til_soundspeed_from_pressure(float density, float P,
                                                 const struct Til_params *mat) {

  float c = sqrtf(mat->A / mat->rho_0);

  return c;
}

#endif /* SWIFT_TILLOTSON_EQUATION_OF_STATE_H */
