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
#ifndef SWIFT_COLIBRE_STAR_FORMATION_STRUCT_H
#define SWIFT_COLIBRE_STAR_FORMATION_STRUCT_H

/* Do we need unique IDs (only useful when spawning
   new particles, conversion gas->stars does not need unique IDs) */
#define star_formation_need_unique_id 0

/**
 * @brief Star-formation-related properties stored in the extended particle
 * data.
 */
struct star_formation_xpart_data {

  /*! Star formation rate */
  float SFR;
};

/**
 * @brief Star-formation-related properties stored in the particle data.
 */
struct star_formation_part_data {

  /*! Velocity dispersion squared in physical internal units */
  float sigma_v2;

  /*! Unweighted gas mass */
  float gas_mass_unweighted;
};

/**
 * @brief Star-formation-related properties stored in the star particle
 * data.
 */
struct star_formation_spart_data {

  /*! The physical birth density */
  float birth_density;

  /*! The birth temperature */
  float birth_temperature;

  /*! The physical subgrid birth density */
  float birth_subgrid_density;

  /*! The subgrid birth temperature */
  float birth_subgrid_temperature;

  /*! The velocity dispersion at birth time in physical internal units */
  float birth_velocity_dispersion;
};

#endif /* SWIFT_COLIBRE_STAR_FORMATION_STRUCT_H */
