/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* This file's header */
#include "stars_properties.h"

/* Local includes */
#include "common_io.h"
#include "hydro_properties.h"
#include "kernel_hydro.h"
#include "restart.h"
#include "units.h"

/**
 * @brief Initialize the global properties of the stars scheme.
 *
 * By default, takes the values provided by the hydro.
 *
 * @param sp The #stars_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param p The already read-in properties of the hydro scheme.
 */
void stars_props_init(struct stars_props *sp,
                      const struct phys_const *phys_const,
                      const struct unit_system *us, struct swift_params *params,
                      const struct hydro_props *p) {

  /* Kernel properties */
  sp->eta_neighbours = parser_get_opt_param_float(
      params, "Stars:resolution_eta", p->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  sp->h_tolerance =
      parser_get_opt_param_float(params, "Stars:h_tolerance", p->h_tolerance);

  /* Get derived properties */
  sp->target_neighbours = pow_dimension(sp->eta_neighbours) * kernel_norm;
  const float delta_eta = sp->eta_neighbours * (1.f + sp->h_tolerance);
  sp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(sp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  sp->max_smoothing_iterations = parser_get_opt_param_int(
      params, "Stars:max_ghost_iterations", p->max_smoothing_iterations);

  /* Properties of the EAGLE feedback model */
  sp->feedback.delta_T =
      parser_get_param_double(params, "EAGLEFeedback:SNII_DeltaT_K");
  sp->feedback.E_SNe_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_erg");
  sp->feedback.SNII_min_mass_Msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Min_mass_Msun");
  sp->feedback.SNII_max_mass_Msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Max_mass_Msun");

  /* Convert the relevant ones to internal units */
  sp->feedback.E_SNe = sp->feedback.E_SNe_cgs /
                       units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
}

/**
 * @brief Print the global properties of the stars scheme.
 *
 * @param sp The #stars_props.
 */
void stars_props_print(const struct stars_props *sp) {

  /* Kernel considerations */
  message("Stars kernel: %s with eta=%f (%.2f neighbours).", kernel_name,
          sp->eta_neighbours, sp->target_neighbours);

  message("Stars relative tolerance in h: %.5f (+/- %.4f neighbours).",
          sp->h_tolerance, sp->delta_neighbours);

  message("Maximal iterations in ghost task set to %d",
          sp->max_smoothing_iterations);

  /* Feedback properties */
  message(
      "Feedback: SNII energy injection of %e erg with a temperature change of "
      "%e K.",
      sp->feedback.E_SNe_cgs, sp->feedback.delta_T);
  message(
      "Feedback: Stars between %.3f and %.3f solar masses are considered to be "
      "SNII.",
      sp->feedback.SNII_min_mass_Msun, sp->feedback.SNII_max_mass_Msun);
}

#if defined(HAVE_HDF5)
void stars_props_print_snapshot(hid_t h_grpstars,
                                const struct stars_props *sp) {

  io_write_attribute_s(h_grpstars, "Kernel function", kernel_name);
  io_write_attribute_f(h_grpstars, "Kernel target N_ngb",
                       sp->target_neighbours);
  io_write_attribute_f(h_grpstars, "Kernel delta N_ngb", sp->delta_neighbours);
  io_write_attribute_f(h_grpstars, "Kernel eta", sp->eta_neighbours);
  io_write_attribute_f(h_grpstars, "Smoothing length tolerance",
                       sp->h_tolerance);
  io_write_attribute_i(h_grpstars, "Max ghost iterations",
                       sp->max_smoothing_iterations);
}

/**
 * @brief Write a #stars_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void stars_props_struct_dump(const struct stars_props *p, FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct stars_props), 1, stream,
                       "starsprops", "stars props");
}

/**
 * @brief Restore a stars_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void stars_props_struct_restore(const struct stars_props *p, FILE *stream) {
  restart_read_blocks((void *)p, sizeof(struct stars_props), 1, stream, NULL,
                      "stars props");
}

#endif
