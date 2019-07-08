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

  /*! Power short time power law */
  double power_short_time;

  /*! Power long time power law */
  double power_long_time;

  /*! Power law break time */
  double break_time_Gyr;

  /*! Normalization time scale */
  double normalization_timescale_Gyr;

  /*! Delay time */
  double delay_time_Gyr;

  /*! DTD normalization short time */
  double norm_short;

  /*! DTD normalization long time */
  double norm_long;
};
