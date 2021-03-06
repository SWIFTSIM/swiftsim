/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_PROPERTIES_GEAR_H
#define SWIFT_RT_PROPERTIES_GEAR_H

/**
 * @file src/rt/GEAR/rt_properties.h
 * @brief Main header file for the 'GEAR' radiative transfer scheme
 * properties.
 */

/**
 * @brief Properties of the 'GEAR' radiative transfer model
 */
struct rt_props {

  /* Are we running with hydro or star controlled injection?
   * This is added to avoid #ifdef macros as far as possible */
  int hydro_controlled_injection;

  /* Are we using constant stellar emission rates? */
  int use_const_emission_rates;

  /* Frequency bin edges for photon groups
   * Includes 0 as leftmost edge, doesn't include infinity as
   * rightmost bin edge*/
  float photon_groups[RT_NGROUPS];

  /* Global constant stellar emission rates */
  double stellar_const_emission_rates[RT_NGROUPS];

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Do extended tests where we assume that all parts
   * have spart neighbours? */
  /* skip this for GEAR */
  /* int debug_do_all_parts_have_stars_checks; */

  /* radiation emitted by stars this step. This is not really a property,
   * but a placeholder to sum up a global variable */
  int debug_radiation_emitted_this_step;

  /* total radiation emitted by stars. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_emitted_tot;

  /* radiation absorbed by gas this step. This is not really a property,
   * but a placeholder to sum up a global variable */
  int debug_radiation_absorbed_this_step;

  /* total radiation absorbed by gas. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_absorbed_tot;
#endif
};

/**
 * @brief Print the RT model.
 *
 * @param rtp The #rt_props
 */
__attribute__((always_inline)) INLINE static void rt_props_print(
    const struct rt_props* rtp) {

  /* Only the master print */
  if (engine_rank != 0) return;

  message("Radiative transfer scheme: '%s'", RT_IMPLEMENTATION);
  char messagestring[200] = "Using photon frequency bins: [0.";
  char freqstring[20];
  for (int g = 1; g < RT_NGROUPS; g++) {
    sprintf(freqstring, ", %.3g", rtp->photon_groups[g]);
    strcat(messagestring, freqstring);
  }
  strcat(messagestring, "]");
  message("%s", messagestring);

  if (rtp->use_const_emission_rates) {
    strcpy(messagestring, "Using constant stellar emission rates: [ ");
    for (int g = 0; g < RT_NGROUPS; g++) {
      sprintf(freqstring, "%.3g ", rtp->stellar_const_emission_rates[g]);
      strcat(messagestring, freqstring);
    }
    strcat(messagestring, "]");
    message("%s", messagestring);
  }
}

/**
 * @brief Initialize the global properties of the RT scheme.
 *
 * @param rtp The #rt_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void rt_props_init(
    struct rt_props* rtp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    struct cosmology* cosmo) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  rtp->hydro_controlled_injection = 1;
#else
  rtp->hydro_controlled_injection = 0;
#endif

  if (RT_NGROUPS <= 0) {
    error("You need to run GEAR-RT with at least 1 photon group, you have %d",
          RT_NGROUPS);
  } else if (RT_NGROUPS == 1) {
    rtp->photon_groups[0] = 0.f;
  } else {
    /* Read in parameters */
    float frequencies[RT_NGROUPS - 1];
    parser_get_param_float_array(params, "GEARRT:photon_groups_Hz",
                                 RT_NGROUPS - 1, frequencies);
    float Hz_internal = units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);
    float Hz_internal_inv = 1.f / Hz_internal;
    for (int g = 0; g < RT_NGROUPS - 1; g++) {
      rtp->photon_groups[g + 1] = frequencies[g] * Hz_internal_inv;
    }
    rtp->photon_groups[0] = 0.f;
  }

  /* Are we using constant emission rates? */
  rtp->use_const_emission_rates = parser_get_opt_param_float(
      params, "GEARRT:use_const_emission_rates", /* default = */ 0);

  if (rtp->use_const_emission_rates) {
    double emission_rates[RT_NGROUPS];
    parser_get_param_double_array(params, "GEARRT:star_emission_rates_LSol",
                                  RT_NGROUPS, emission_rates);
    const double unit_power = units_cgs_conversion_factor(us, UNIT_CONV_POWER);
    const double unit_power_inv = 1. / unit_power;
    for (int g = 0; g < RT_NGROUPS; g++) {
      rtp->stellar_const_emission_rates[g] = emission_rates[g] * unit_power_inv;
    }
  } else {
    /* kill the run for now */
    error("GEAR-RT can't run without constant stellar emission rates for now.");
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  rtp->debug_radiation_emitted_tot = 0ULL;
  rtp->debug_radiation_absorbed_tot = 0ULL;
#endif

  /* After initialisation, print params to screen */
  rt_props_print(rtp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Radiative transfer initialized");
  }
}

/**
 * @brief Write an RT properties struct to the given FILE as a
 * stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_dump(
    const struct rt_props* props, FILE* stream) {

  restart_write_blocks((void*)props, sizeof(struct rt_props), 1, stream,
                       "RT props", "RT properties struct");
}

/**
 * @brief Restore an RT properties struct from the given FILE as
 * a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_restore(
    struct rt_props* props, FILE* stream) {

  restart_read_blocks((void*)props, sizeof(struct rt_props), 1, stream, NULL,
                      "RT properties struct");
}

#endif /* SWIFT_RT_PROPERTIES_GEAR_H */
