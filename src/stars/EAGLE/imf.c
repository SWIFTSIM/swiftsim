/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @file src/stars/EAGLE/imf.c
 * @brief EAGLE functions related to the IMF.
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes */
#include "error.h"
#include "imf.h"

/**
 *
 * @param imf_model The #eagle_stars_imf in use.
 * @param age Age of the stars in giga-years.
 * @param Z Metallicity (metal mass fraction) of the stars.
 */
double eagle_stars_dying_mass(const struct eagle_stars_imf *imf_model,
                              const double age, const double Z) {

  switch (imf_model->lifetime_model) {

    case eagle_stars_lifetime_P98: {
      if (age <= 0.) {
        return imf_model->maximal_mass;
      }

      /* Logarithm of age in years */
      const double log_age = log10(age) + 9.;

    } break;

    default:
      error("Invalid life-time model! Not implemented yet!");
  }

  return 0.;
}
