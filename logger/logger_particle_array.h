/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
 * @file logger_particle_array.h
 * @brief This file contains a structure that contains the different type of
 * particles.
 */

#ifndef LOGGER_LOGGER_PARTICLE_ARRAY_H
#define LOGGER_LOGGER_PARTICLE_ARRAY_H

/* Include config */
#include "../config.h"

/* Include local files */
#include "logger_particle.h"

/**
 * @brief Structure containing all the particles.
 *
 * This is used for the interface with the user.
 */
struct logger_particle_array {
  /* Hydro */
  struct {
    /* The array of particles */
    struct logger_particle *parts;

    /* The number of particles */
    size_t n;

    /* Size of the allocated buffer */
    size_t allocated_size;

  } hydro;

  /* Gravity */
  struct {
    /* The array of particles */
    struct logger_gparticle *parts;

    /* The number of particles */
    size_t n;

    /* Size of the allocated buffer */
    size_t allocated_size;
  } grav;

  /* Stars */
  struct {
    /* The array of particles */
    struct logger_sparticle *parts;

    /* The number of particles */
    size_t n;

    /* Size of the allocated buffer */
    size_t allocated_size;
  } stars;
};

void logger_particle_array_init(struct logger_particle_array *array);
void logger_particle_array_allocate(struct logger_particle_array *array,
                                    size_t n_part, size_t n_gpart,
                                    size_t n_spart, const int empty);
void logger_particle_array_free(struct logger_particle_array *array);
void logger_particle_array_change_size(struct logger_particle_array *array,
                                       size_t new_n_part, size_t new_n_gpart,
                                       size_t new_n_spart);

void logger_particle_array_add_hydro(struct logger_particle_array *array,
                                     size_t offset);
void logger_particle_array_add_stars(struct logger_particle_array *array,
                                     size_t offset);
void logger_particle_array_add_gravity(struct logger_particle_array *array,
                                       size_t offset);

void logger_particle_array_update(struct logger_particle_array *prev,
                                  struct logger_particle_array *next,
                                  struct logger_particle_array *tmp,
                                  size_t n_deleted_hydro, size_t n_deleted_grav,
                                  size_t n_deleted_stars);

#endif  // LOGGER_LOGGER_PARTICLE_ARRAY_H
