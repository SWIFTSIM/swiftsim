/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_DTD_H
#define SWIFT_DTD_H

/**
 * @file src/snia_dtd.h
 * @brief Branches between the different SNIa delay time distributions recipies.
 */

/* Config parameters. */
#include "../config.h"

#include "DTD_struct.h"
#include "feedback_properties.h"
#include "parser.h"
#include "physical_constants.h"
#include "units.h"

/* Import the right DTD definition */
#if defined(SNIA_DTD_EXP)

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * We follow Foerster et al. 2006, MNRAS, 368
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                        const double t1,
                                        const struct feedback_props* fp) {

  /* Return zero if below delay time */
  if (t1 < fp->dtd_data.delay_time_Gyr) return 0;

  /* Start from the correct time */
  const double tmin = max(t0, fp->dtd_data.delay_time_Gyr);

  /* The calculation is written as the integral between t0 and t1 of
   * eq. 3 of Schaye 2015 paper. */
  const double tau_inv = fp->dtd_data.SNIa_timescale_Gyr_inv;
  const double nu = fp->dtd_data.SNIa_efficiency;
  const double num_SNIa_per_Msun =
      nu * (exp(-tmin * tau_inv) - exp(-t1 * tau_inv));

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp,
                            const struct phys_const* phys_const,
                            const struct unit_system* us,
                            struct swift_params* params) {

  /* Get the SNIa efficiency */
  fp->dtd_data.SNIa_efficiency =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  /* Get the exponential delay time for the DTD */
  fp->dtd_data.SNIa_timescale_Gyr =
      parser_get_param_float(params, "SNIaDTD:SNIa_timescale_Gyr");

  /* Calculate the inverse of the exponential delay time */
  fp->dtd_data.SNIa_timescale_Gyr_inv = 1.f / fp->dtd_data.SNIa_timescale_Gyr;

  /* Set the delay time of the DTD */
  fp->dtd_data.delay_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  /* Properly normalize the exponential DTD */
  fp->dtd_data.SNIa_efficiency *=
      exp(-fp->dtd_data.delay_time_Gyr * fp->dtd_data.SNIa_timescale_Gyr_inv);
}

#elif defined(SNIA_DTD_POWER)

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * This model assumes that the SNIa DTD is given by a power law, this is
 * a common DTD model, Moaz & Mannucci (2012), PASA, 29, 447 gives an
 * overview of the different variables found in the literature
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                        const double t1,
                                        const struct feedback_props* fp) {

  /* Return zero if below delay time */
  if (t1 < fp->dtd_data.delay_time_Gyr) return 0;

  /* Start from the correct time */
  const double tmin = max(t0, fp->dtd_data.delay_time_Gyr);

  /* The calculation is written as the integral between t0 and t1 of
   * a constant DTD given by \nu / \tau */
  const double norm = fp->dtd_data.norm;
  const double power = fp->dtd_data.power;
  const double num_SNIa_per_Msun =
      norm * (pow(t1, 1. - power) - pow(tmin, 1. - power));

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp,
                            const struct phys_const* phys_const,
                            const struct unit_system* us,
                            struct swift_params* params) {

  /* Get the SNIa efficiency */
  fp->dtd_data.SNIa_efficiency =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  /* Get the power of the power law */
  fp->dtd_data.power =
      parser_get_param_double(params, "SNIaDTD:power_law_slope");

  /* Get the normalization time over which the DTD is normalized */
  fp->dtd_data.normalization_timescale_Gyr =
      parser_get_param_double(params, "SNIaDTD:normalization_timescale_Gyr");

  /* Get the delay time */
  fp->dtd_data.delay_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  /* Calculate the normalization of the power-law DTD */
  const double below_frac =
      pow(fp->dtd_data.normalization_timescale_Gyr, 1 - fp->dtd_data.power) -
      pow(fp->dtd_data.delay_time_Gyr, 1 - fp->dtd_data.power);
  fp->dtd_data.norm = fp->dtd_data.SNIa_efficiency / below_frac;
}

#elif defined(SNIA_DTD_POWER_BETA1)

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * This model assumes that the SNIa DTD is a power law with \beta = 1,
 * this model is a special case of the power law model because the
 * integrals have a different functional form than the general power
 * law. A lot of observations seem to be close to a power law like this,
 * see Moaz & Mannucci (2012), PASA, 29, 447 for a review
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                        const double t1,
                                        const struct feedback_props* fp) {

  /* Return zero if below delay time */
  if (t1 < fp->dtd_data.delay_time_Gyr) return 0;

  /* Start from the correct time */
  const double tmin = max(t0, fp->dtd_data.delay_time_Gyr);

  /* The calculation is written as the integral between t0 and t1 of
   * a power law DTD given by ~\nu /t  */
  const double norm = fp->dtd_data.norm;
  const double num_SNIa_per_Msun = norm * log(t1 / tmin);

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp,
                            const struct phys_const* phys_const,
                            const struct unit_system* us,
                            struct swift_params* params) {

  /* Get the SNIa efficiency */
  fp->dtd_data.SNIa_efficiency =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  /* Get the SNIa normalization timescale */
  fp->dtd_data.normalization_timescale_Gyr =
      parser_get_param_double(params, "SNIaDTD:normalization_timescale_Gyr");

  /* Get the delay time of the DTD */
  fp->dtd_data.delay_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  /* If the delay time is zero, exit with an error */
  if (fp->dtd_data.delay_time_Gyr == 0.) error("Delay time cannot be zero");

  /* Calculate the normalization of the power law DTD with beta=1 */
  const double below_frac = log(fp->dtd_data.normalization_timescale_Gyr /
                                fp->dtd_data.delay_time_Gyr);
  fp->dtd_data.norm = fp->dtd_data.SNIa_efficiency / below_frac;
}

#elif defined(SNIA_DTD_GAUSSIAN)

#include <math.h>

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * This model assumes that the SNIa DTD is a Gaussian following for
 * example Dahlen et al. 2004, ApJ, 613, 189.
 * There is the option of using a constant besided the Gaussian this is
 * following the approach adopted by the Fire 2 simulations (Hopkins et al.
 * 2018, 480, 800)
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                        const double t1,
                                        const struct feedback_props* fp) {

  /* Return zero if below delay time */
  if (t1 < fp->dtd_data.delay_time_Gyr) return 0;

  /* Start from the correct time */
  const double tmin = max(t0, fp->dtd_data.delay_time_Gyr);

  /* The calculation is written as the integral between t0 and t1 of
   * a constant DTD given by \nu / \tau */
  
  /* First calculate the number of SNIa for the constant part of the DTD: */
  const double norm = fp->dtd_data.norm_const;
  const double num_SNIa_per_Msun_const = norm * (t1 - tmin);

  /* Calculate the number of SNIa for the gaussian part of the DTD: */
  const double nu_gauss = fp->dtd_data.SNIa_efficiency_gauss;

  /* get one over the standard deviation time */
  const double inv_std = fp->dtd_data.std_characteristic_time_Gyr_inv;

  /* Get the difference with the characteristic time of the Gaussian 
   * for both the minimum and maximum time on this interval: */
  const double tdif0 = tmin - fp->dtd_data.characteristic_time_Gyr;
  const double tdif1 = t1 - fp->dtd_data.characteristic_time_Gyr;
  
  /* calculate the upper term of the integral (not factor 2) */
  const double gauss_up = erf(tdif1 * inv_std);
  
  /* The lower part */
  const double gauss_low = erf(tdif0 * inv_std);

  /* Finish the integral */
  const double integral = .5 * (gauss_up - gauss_low);

  /* The number of SNIa from the Gaussian part */
  const double num_SNIa_per_Msun_gauss = nu_gauss * integral;

  /* calculate the total Gaussian + constant */
  const double num_SNIa_per_Msun = num_SNIa_per_Msun_gauss + num_SNIa_per_Msun_const;

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp,
                            const struct phys_const* phys_const,
                            const struct unit_system* us,
                            struct swift_params* params) {

  /* Set the SNIa efficiency for the constant part of the DTD */
  fp->dtd_data.SNIa_efficiency_const =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_const_p_Msun");

  /* Set the SNIa efficiency for the Gaussian part of the DTD */
  fp->dtd_data.SNIa_efficiency_gauss =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_gauss_p_Msun");

  /* Set the normalization time for the constant part of the DTD */
  fp->dtd_data.normalization_timescale_Gyr =
      parser_get_param_double(params, "SNIaDTD:normalization_timescale_Gyr");

  /* Set the delay time of the DTD */
  fp->dtd_data.delay_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  /* Set the characteristic time of the Gaussian part of the DTD */
  fp->dtd_data.characteristic_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:characteristic_time_Gyr");

  /* Set the standard deviation of the Gaussian part of the DTD */
  fp->dtd_data.std_characteristic_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:STD_characteristic_time_Gyr");

  /* Calculate the inverse of the standard deviation of the Gaussian part
   * divided by sqrt 2 */
  fp->dtd_data.std_characteristic_time_Gyr_inv =
      1. / fp->dtd_data.std_characteristic_time_Gyr / sqrt(2);

  /* Calculate the norm of the constant part of the DTD */
  fp->dtd_data.norm_const = fp->dtd_data.SNIa_efficiency_const /
                            fp->dtd_data.normalization_timescale_Gyr;
}

#elif defined(SNIA_DTD_CONST)

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * This model assumes that the SNIa DTD is constant.
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                        const double t1,
                                        const struct feedback_props* fp) {

  /* Return zero if below delay time */
  if (t1 < fp->dtd_data.delay_time_Gyr) return 0;

  /* Start from the correct time */
  const double tmin = max(t0, fp->dtd_data.delay_time_Gyr);

  /* The calculation is written as the integral between t0 and t1 of
   * a constant DTD given by \nu / \tau */
  const double norm = fp->dtd_data.norm;
  const double num_SNIa_per_Msun = norm * (t1 - tmin);

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp,
                            const struct phys_const* phys_const,
                            const struct unit_system* us,
                            struct swift_params* params) {

  fp->dtd_data.SNIa_efficiency =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  fp->dtd_data.normalization_timescale_Gyr =
      parser_get_param_float(params, "SNIaDTD:normalization_timescale_Gyr");

  fp->dtd_data.norm =
      fp->dtd_data.SNIa_efficiency / fp->dtd_data.normalization_timescale_Gyr;

  /* Set the delay time of the DTD */
  fp->dtd_data.delay_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");
}

#elif defined(SNIA_DTD_BROKEN_POWER_LAW)

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * This model assumes that the SNIa DTD is a broken power law. This shape
 * is common in theoretical models of the DTD in the double degenerate (DD)
 * scenario, in which the SNIa DTD has a shallow slope below a break time
 * and a deeper slope after a break time. For a review see Moaz & Mannucci
 * (2012), PASA, 29, 447
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
static inline double dtd_number_of_SNIa(const struct spart* sp, const double t0,
                                        const double t1,
                                        const struct feedback_props* fp) {

  /* Return zero if below delay time */
  if (t1 < fp->dtd_data.delay_time_Gyr) return 0;

  /* Start from the correct time */
  const double tmin = max(t0, fp->dtd_data.delay_time_Gyr);

  /* The calculation is written as the integral between t0 and t1 of
   * a constant DTD given by \nu / \tau */
  const double tb = fp->dtd_data.break_time_Gyr;

  /* Check if we are in the long time regime */
  if (t0 > tb) {
    const double norm = fp->dtd_data.norm_long;
    const double power = fp->dtd_data.power_long_time;
    const double num_SNIa_per_Msun =
        norm * (pow(t1, 1. - power) - pow(t0, 1. - power));
    return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
  }

  /* Check if we are in the short time regime */
  if (t1 < tb) {
    const double norm = fp->dtd_data.norm_short;
    const double power = fp->dtd_data.power_short_time;
    const double num_SNIa_per_Msun =
        norm * (pow(t1, 1. - power) - pow(tmin, 1. - power));
    return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
  }

  /* Now we are in the intermediate regime that requires more calculations */
  const double power_short = fp->dtd_data.power_short_time;
  const double power_long = fp->dtd_data.power_long_time;
  const double norm_short = fp->dtd_data.norm_short;
  const double norm_long = fp->dtd_data.norm_long;

  const double num_SNIa_per_Msun =
      norm_short * (pow(tb, 1. - power_short) - pow(tmin, 1. - power_short)) +
      norm_long * (pow(tb, 1. - power_long) - pow(t1, 1. - power_long));

  return num_SNIa_per_Msun * sp->mass_init * fp->mass_to_solar_mass;
}

/**
 * @brief initialize the DTD
 *
 * @param fp the properties of the stellar model.
 * @param phys_const the constant
 * @param us the unit system
 * @param params the input parameters
 */
static inline void dtd_init(struct feedback_props* fp,
                            const struct phys_const* phys_const,
                            const struct unit_system* us,
                            struct swift_params* params) {

  /* Get the SNIa efficiency */
  fp->dtd_data.SNIa_efficiency =
      parser_get_param_float(params, "SNIaDTD:SNIa_efficiency_p_Msun");

  /* Get the short time power law slope */
  fp->dtd_data.power_short_time =
      parser_get_param_double(params, "SNIaDTD:power_law_slope_short_time");

  /* Get the long time power law slope */
  fp->dtd_data.power_long_time =
      parser_get_param_double(params, "SNIaDTD:power_law_slope_long_time");

  /* Get the normalization time over which the DTD is normalized */
  fp->dtd_data.normalization_timescale_Gyr =
      parser_get_param_double(params, "SNIaDTD:normalization_timescale_Gyr");

  /* Get the delay time */
  fp->dtd_data.delay_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:SNIa_delay_time_Gyr");

  /* Get the break time */
  fp->dtd_data.break_time_Gyr =
      parser_get_param_double(params, "SNIaDTD:break_time_Gyr");

  /* Calculate the normalization of the power-law DTD */
  
  /* Calculate 1 minus the short time power */
  const double one_minus_power_short = 1. - fp->dtd_data.power_short_time;
  /* Calculate one over this number */
  const double one_minus_power_short_inv = 1./one_minus_power_short;

  /* Get the delay time*/
  const double t_delay = fp->dtd_data.delay_time_Gyr;

  /* Get the break time */
  const  double t_break = fp->dtd_data.break_time_Gyr;

  /* Calculate the first norm */
  const double norm1_inv = one_minus_power_short_inv * (1. - pow(t_delay/t_break, one_minus_power_short));

  /* Calculate 1 minus the long time power */
  const double one_minus_power_long = 1. - fp->dtd_data.power_short_time;
  /* Calculate one over this number */
  const double one_minus_power_long_inv = 1./one_minus_power_short;

  const double t_norm = fp->dtd_data.normalization_timescale_Gyr;

  /* Calculate the second norm */
  const double norm2_inv = one_minus_power_long_inv * (pow(t_norm/t_break, one_minus_power_long) - 1.);

  /* Store the normalization */
  fp->dtd_data.norm_short = fp->dtd_data.SNIa_efficiency /
                            (norm1_inv + norm2_inv) * one_minus_power_short_inv;
  fp->dtd_data.norm_long = fp->dtd_data.SNIa_efficiency /
                           (norm1_inv + norm2_inv) * one_minus_power_long_inv;
}

#else
#error "Invalid choice of the SNIa delay time distribution (DTD)"
#endif

#endif /* SWIFT_DTD_H */
