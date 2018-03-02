/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_GRACKLE_IO_H
#define SWIFT_COOLING_GRACKLE_IO_H

#include "../config.h"
#include "cooling_struct.h"
#include "io_properties.h"

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Cooling Model", "Grackle");
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param cooling The #cooling_function_data
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct xpart* xparts, struct io_props* list,
    const struct cooling_function_data *cooling) {

  int num = 0;

  if (cooling->output_mode == 0)
    return num;

#if COOLING_GRACKLE_MODE >= 1
  /* List what we want to write */
  list[0] = io_make_output_field("HI", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HI_frac);

  list[1] = io_make_output_field("HII", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HII_frac);
  
  list[2] = io_make_output_field("HeI", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HeI_frac);

  list[3] = io_make_output_field("HeII", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HeII_frac);

  list[4] = io_make_output_field("HeIII", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HeIII_frac);

  list[5] = io_make_output_field("e", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.e_frac);

  num += 6;
#endif

  if (cooling->output_mode == 1)
    return num;

#if COOLING_GRACKLE_MODE >= 2
  list += num;

  list[0] = io_make_output_field("HM", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HM_frac);
    
  list[1] = io_make_output_field("H2I", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.H2I_frac);

  list[2] = io_make_output_field("H2II", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.H2II_frac);

  num += 3;
#endif

  if (cooling->output_mode == 2)
    return num;

#if COOLING_GRACKLE_MODE >= 3
  list += num;

  list[0] = io_make_output_field("DI", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.DI_frac);
    
  list[1] = io_make_output_field("DII", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.DII_frac);

  list[2] = io_make_output_field("HDI", FLOAT, 1, UNIT_CONV_NO_UNITS,
                                 xparts, cooling_data.HDI_frac);

  num += 3;
#endif

  return num;

}

/**
 * @brief Parser the parameter file and initialize the #cooling_function_data
 * @param parameter_file The parser parameter file
 * @param cooling The cooling properties to initialize
 */
__attribute__((always_inline)) INLINE static void cooling_parse_arguments(
    const struct swift_params* parameter_file,
    struct cooling_function_data* cooling) {

  parser_get_param_string(parameter_file, "GrackleCooling:CloudyTable",
                          cooling->cloudy_table);
  cooling->with_uv_background =
      parser_get_param_int(parameter_file, "GrackleCooling:WithUVbackground");

  cooling->redshift =
      parser_get_param_double(parameter_file, "GrackleCooling:Redshift");

  cooling->with_metal_cooling =
    parser_get_param_int(parameter_file, "GrackleCooling:WithMetalCooling");

  cooling->provide_volumetric_heating_rates =
    parser_get_param_int(parameter_file, "GrackleCooling:ProvideVolumetricHeatingRates");

  cooling->provide_specific_heating_rates =
    parser_get_param_int(parameter_file, "GrackleCooling:ProvideSpecificHeatingRates");

  cooling->self_shielding_method =
    parser_get_param_int(parameter_file, "GrackleCooling:SelfShieldingMethod");

  cooling->output_mode =
    parser_get_param_int(parameter_file, "GrackleCooling:OutputMode");
}


#endif // SWIFT_COOLING_GRACKLE_IO_H
