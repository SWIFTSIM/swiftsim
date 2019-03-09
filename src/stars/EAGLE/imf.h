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
#ifndef SWIFT_STARS_EAGLE_IMF_H
#define SWIFT_STARS_EAGLE_IMF_H

/**
 * @brief Stellar inital mass function.
 */
enum eagle_stars_imf_model {
  eagle_stars_imf_Chabrier, /*!< Chabrier 2003 */
  eagle_stars_imf_PowerLaw, /*!< Power-law IMF */
};

/**
 * @brief Model used for the stars' life time.
 */
enum eagle_stars_lifetime_model {

  eagle_stars_lifetime_PM93, /*!< Padovani & Mateucci 1993 */
  eagle_stars_lifetime_MM89, /*!< Maeder & Meynet 1989 */
  eagle_stars_lifetime_P98   /*!< Portinari et al. 1998 */
};

struct eagle_stars_imf {

  /*! Minimal mass considered in solar masses */
  double minimal_mass;

  /*! Maximal mass considered in solar masses */
  double maximal_mass;

  /*! The life-time model (isochrones) */
  enum eagle_stars_lifetime_model lifetime_model;

  /*! The IMF shape */
  enum eagle_stars_imf_model imf;
};

#endif /* SWIFT_STARS_EAGLE_IMF_H */
