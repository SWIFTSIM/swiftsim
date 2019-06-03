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
#include "stellar_evolution_io.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Properties of the GEAR feedback model.
 */
struct feedback_props {
  /*! Energy per supernovae */
  float energy_per_supernovae;

  /*! Thermal time */
  float thermal_time;

  /*! filename of the chemistry table */
  char filename[PARSER_MAX_LINE_SIZE];

  /*! The stellar model */
  struct stellar_model stellar_model;
};

__attribute__((always_inline)) INLINE static void feedback_read_tables(
    struct feedback_props* fp, const struct phys_const* phys_const,
    const struct unit_system* us);

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * By default, takes the values provided by the hydro.
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
  double e_feedback = parser_get_param_double(params, "GEARFeedback:SupernovaeEnergy_erg");
  e_feedback /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  fp->energy_per_supernovae = e_feedback;

  /* Thermal time */
  fp->thermal_time =
    parser_get_param_float(params, "GEARFeedback:ThermalTime_Myr");
  fp->thermal_time *= phys_const->const_year * 1e6;

  /* filename of the chemistry table */
  parser_get_param_string(params, "GEARFeedback:ChemistryTable", fp->filename);

  /* Read tables */
  // feedback_read_tables(fp, phys_const, us);

  /* Initialize the stellar model */
  stellar_evolution_props_init(&fp->stellar_model, phys_const,
			       us, params, cosmo);

  /* Print a final message. */
  message("initialized stellar feedback");
}

#ifdef HAVE_HDF5
/**
 * @brief Initialize the initial mass function.
 *
 * @param fp The #feedback_props.
 * @param file_id The HDF5 file to use.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE static void feedback_initialize_initial_mass_function(
    struct feedback_props *fp, hid_t file_id,
    const struct phys_const *phys_const, const struct unit_system *us) {
}

/**
 * @brief Initialize the stellar lifetime.
 *
 * @param fp The #feedback_props.
 * @param file_id The HDF5 file to use.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE static void feedback_initialize_lifetime(
    struct feedback_props *fp, hid_t file_id,
    const struct phys_const *phys_const, const struct unit_system *us) {
}

/**
 * @brief Initialize the SNIa.
 *
 * @param fp The #feedback_props.
 * @param file_id The HDF5 file to use.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE static void feedback_initialize_supernovae_ia(
    struct feedback_props *fp, hid_t file_id,
    const struct phys_const *phys_const, const struct unit_system *us) {
}

/**
 * @brief Initialize the SNII.
 *
 * @param fp The #feedback_props.
 * @param file_id The HDF5 file to use.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE static void feedback_initialize_supernovae_ii(
    struct feedback_props *fp, hid_t file_id,
    const struct phys_const *phys_const, const struct unit_system *us) {
}

/**
 * @brief Reads the chemistry table
 *
 * @param fp The #feedback_props.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 */
__attribute__((always_inline)) INLINE static void feedback_read_tables(
    struct feedback_props* fp, const struct phys_const* phys_const,
    const struct unit_system* us) {

  hid_t status;

  /* Open file */
  hid_t file_id = H5Fopen(fp->filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("Unable to open file %s", fp->filename);

  /* Initialize the IMF */
  feedback_initialize_initial_mass_function(fp, file_id, phys_const, us);

  /* Initialize the lifetime */
  feedback_initialize_lifetime(fp, file_id, phys_const, us);

  /* Initialize the SNIa */
  feedback_initialize_supernovae_ia(fp, file_id, phys_const, us);

  /* Initialize the SNII */
  feedback_initialize_supernovae_ii(fp, file_id, phys_const, us);

  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");

}
#else // HAVE_HDF5
#error Cannot use GEAR feedback without HDF5
#endif


#endif /* SWIFT_GEAR_FEEDBACK_PROPERTIES_H */
