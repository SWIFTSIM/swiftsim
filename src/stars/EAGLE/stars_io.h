/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STARS_IO_H
#define SWIFT_EAGLE_STARS_IO_H

#include "io_properties.h"
#include "stars_part.h"

/**
 * @brief Specifies which s-particle fields to read from a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void stars_read_particles(struct spart *sparts,
                                        struct io_props *list,
                                        int *num_fields) {

  /* Say how much we want to read */
  *num_fields = 5;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, sparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, sparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, sparts, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, sparts, h);
}

/**
 * @brief Specifies which s-particle fields to write to a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void stars_write_particles(const struct spart *sparts,
                                         struct io_props *list,
                                         int *num_fields) {

  /* Say how much we want to write */
  *num_fields = 10;

  /* List what we want to write */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 sparts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, sparts, v);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, sparts, mass);
  list[3] = io_make_output_field("ParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 sparts, id);
  list[4] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 sparts, h);
  list[5] = io_make_output_field("BirthDensity", FLOAT, 1, UNIT_CONV_DENSITY,
                                 sparts, birth_density);
  list[6] = io_make_output_field("InitialMasses", FLOAT, 1, UNIT_CONV_MASS,
                                 sparts, mass_init);
  list[7] = io_make_output_field("BirthTime", FLOAT, 1, UNIT_CONV_TIME, sparts,
                                 birth_time);
  list[8] = io_make_output_field("NeighbourGasDensity", FLOAT, 1,
                                 UNIT_CONV_DENSITY, sparts, rho_gas);
  list[9] = io_make_output_field("NeighbourGasMass", FLOAT, 1, UNIT_CONV_MASS,
                                 sparts, density.neighbour_mass);
}

#endif /* SWIFT_EAGLE_STAR_IO_H */
