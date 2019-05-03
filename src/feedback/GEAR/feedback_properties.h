/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @brief Properties of the GEAR feedback model.
 */
struct feedback_props {
  /*! Energy per supernovae */
  float energy_per_supernovae;

  /*! Thermal time */
  float thermal_time;

  /*! Mass ejected per supernovae */
  float mass_ejected;
};


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

  /* Mass ejection */
  double mass_ejected =
    parser_get_param_float(params, "GEARFeedback:MassEjected_Msun");
  mass_ejected *= phys_const->const_solar_mass;
  fp->mass_ejected = mass_ejected;

  /* Print a final message. */
  message("initialized stellar feedback");
}

#endif /* SWIFT_GEAR_FEEDBACK_PROPERTIES_H */
