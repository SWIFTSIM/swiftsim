/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 * @brief SNIa struct for DTD properties
 */
struct dtd {
  
  /*! Efficiency of the SNIa model */
  double SNIa_efficiency;

  /*! Power of power law */
  double power;

  /*! Normalization time */
  double normalization_timescale_Gyr;

  /*! Delay time */
  double delay_time_Gyr;

  /*! DTD normalization */
  double norm;
};
