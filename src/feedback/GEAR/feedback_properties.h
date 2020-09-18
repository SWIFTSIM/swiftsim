/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_FEEDBACK_PROPERTIES_H
#define SWIFT_GEAR_FEEDBACK_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Properties of the GEAR feedback model.
 */
struct feedback_props {
  /*! Energy per supernovae */
  float energy_per_supernovae;

  /*! The stellar model */
  struct stellar_model stellar_model;

  /*! The stellar model for first stars */
  struct stellar_model stellar_model_first_stars;

  /* Metallicity limits for the first stars */
  float metallicity_max_first_stars;
};

/**
 * @brief Print the feedback model.
 *
 * @param feedback_props The #feedback_props
 */
__attribute__((always_inline)) INLINE static void feedback_props_print(
    const struct feedback_props* feedback_props) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  /* Print the name of the elements */
  char txt[GEAR_CHEMISTRY_ELEMENT_COUNT * (GEAR_LABELS_SIZE + 2)] = "";
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (i != 0) {
      strcat(txt, ", ");
    }
    strcat(txt, stellar_evolution_get_element_name(
        &feedback_props->stellar_model, i));
  }

  if (engine_rank == 0) {
    message("Chemistry elements: %s", txt);
  }

  /* Print the feedback properties */
  message("Energy per supernovae = %.2g",
          feedback_props->energy_per_supernovae);
  message("Yields table = %s", feedback_props->stellar_model.yields_table);

  /* Print the stellar model */
  stellar_model_print(&feedback_props->stellar_model);

  /* Print the first stars */
  if (feedback_props->metallicity_max_first_stars != -1) {
    message("Yields table first stars = %s",
            feedback_props->stellar_model_first_stars.yields_table);
    stellar_model_print(&feedback_props->stellar_model_first_stars);
    message("Metallicity max for the first stars (in mass fraction) = %g",
            feedback_props->metallicity_max_first_stars);
  }
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void feedback_props_init(
    struct feedback_props* fp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo) {

  /* Supernovae energy */
  double e_feedback =
      parser_get_param_double(params, "GEARFeedback:supernovae_energy_erg");
  e_feedback /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  fp->energy_per_supernovae = e_feedback;

  /* filename of the chemistry tables. */
  parser_get_param_string(params, "GEARFeedback:yields_table", fp->stellar_model.yields_table);

  /* Initialize the stellar models. */
  stellar_evolution_props_init(&fp->stellar_model, phys_const, us, params,
                               cosmo);

  /* Now the same for first stars. */
  fp->metallicity_max_first_stars = parser_get_opt_param_float(
      params, "GEARFeedback:metallicity_max_first_stars", -1);

  if (fp->metallicity_max_first_stars == -1) {
    message("First stars are disabled.");
  }
  else {
    if (fp->metallicity_max_first_stars < 0) {
      error("The metallicity threshold for the first stars is in mass fraction. "
            "It cannot be lower than 0.");
    }
    message("Reading the stellar model for the first stars");
    parser_get_param_string(params, "GEARFeedback:yields_table_first_stars",
                            fp->stellar_model_first_stars.yields_table);
    stellar_evolution_props_init(&fp->stellar_model_first_stars,
                                 phys_const, us, params, cosmo);
  }

  /* Print the stellar properties */
  feedback_props_print(fp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Stellar feedback initialized");
  }
}

#endif /* SWIFT_GEAR_FEEDBACK_PROPERTIES_H */
