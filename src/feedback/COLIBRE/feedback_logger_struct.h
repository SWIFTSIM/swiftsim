/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_COLIBRE_FEEDBACK_LOGGER_STRUCT_H
#define SWIFT_COLIBRE_FEEDBACK_LOGGER_STRUCT_H

/* feedback history struct */
struct feedback_history_SNIa {

  /*! Total new SNIa injected energy */
  double SNIa_energy;

  /*! Number of SNIa */
  float N_SNIa;

  /*! Number of heating events */
  int heating;

};

/* Global variable for the SNIa events */
extern struct feedback_history_SNIa log_SNIa;

/* feedback history struct accumulator */
struct feedback_history_accumulator {

  /*! Previous storage step*/
  int step_prev;

  /*! Previous storage time */
  float time_prev;

  /*! Previous storage redshift */
  float z_prev;

  /*! Previous storage scale factor */
  float a_prev;

};

#endif /* SWIFT_COLIBRE_FEEDBACK_LOGGER_STRUCT_H */
