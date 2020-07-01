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

/* Include corresponding header */
#include "logger_particle_array.h"

/**
 * @brief Initialize the #logger_particle_array structure.
 *
 * @param array The array to initialize.
 */
void logger_particle_array_init(struct logger_particle_array *array) {
  /* Hydro */
  array->hydro.n = 0;
  array->hydro.allocated_size = 0;
  array->hydro.parts = NULL;

  /* Dark matter */
  array->grav.n = 0;
  array->grav.allocated_size = 0;
  array->grav.parts = NULL;

  /* Stars */
  array->stars.n = 0;
  array->stars.allocated_size = 0;
  array->stars.parts = NULL;
}

/**
 * @brief The allocate the required memory.
 *
 * @param array The array to allocate.
 * @param n_part The number of #logger_particle.
 * @param n_gpart The number of #logger_gparticle.
 * @param n_spart The number of #logger_sparticle.
 * @param empty Are we considering that the array is empty?
 */
void logger_particle_array_allocate(struct logger_particle_array *array,
                                    size_t n_part, size_t n_gpart,
                                    size_t n_spart, const int empty) {
  /* Hydro */
  array->hydro.n = empty ? 0 : n_part;
  array->hydro.allocated_size = n_part;
  array->hydro.parts =
      (struct logger_particle *)malloc(n_part * sizeof(struct logger_particle));
  if (array->hydro.parts == NULL) {
    error("Failed to allocate the hydro particles");
  }

  /* Gravity */
  array->grav.n = empty ? 0 : n_gpart;
  array->grav.allocated_size = n_gpart;
  array->grav.parts = (struct logger_gparticle *)malloc(
      n_gpart * sizeof(struct logger_gparticle));
  if (array->grav.parts == NULL) {
    error("Failed to allocate the gravity particles");
  }

  /* Stars */
  array->stars.n = empty ? 0 : n_spart;
  array->stars.allocated_size = n_spart;
  array->stars.parts = (struct logger_sparticle *)malloc(
      n_spart * sizeof(struct logger_sparticle));
  if (array->stars.parts == NULL) {
    error("Failed to allocate the stars particles");
  }
}

/**
 * @brief Free the allocated memory.
 *
 * @param array The array to free.
 */
void logger_particle_array_free(struct logger_particle_array *array) {
  free(array->hydro.parts);
  free(array->grav.parts);
  free(array->stars.parts);

  logger_particle_array_init(array);
}

/**
 * @brief Change the size of the allocated memory.
 *
 * @param array The array to increase.
 * @param new_n_part The new number of (hydro) particles.
 * @param new_n_gpart The new number of (gravity) particles.
 * @param new_n_spart The new number of (stars) particles.
 */
void logger_particle_array_change_size(struct logger_particle_array *array,
                                       size_t new_n_part, size_t new_n_gpart,
                                       size_t new_n_spart) {
  /* Hydro */
  /* Check if need to free memory */
  if (new_n_part == 0 && array->hydro.n != 0) {
    free(array->hydro.parts);
    array->hydro.parts = NULL;
  }
  /* Check if need to change the size */
  else if (new_n_part != array->hydro.n) {

    array->hydro.allocated_size = new_n_part;

    /* Allocate the new array */
    struct logger_particle *parts = (struct logger_particle *)malloc(
        new_n_part * sizeof(struct logger_particle));
    if (parts == NULL) {
      error("Failed to allocate the hydro particles.");
    }

    /* Copy the previous particles */
    if (array->hydro.n != 0) {
      memcpy(parts, array->hydro.parts,
             array->hydro.n * sizeof(struct logger_particle));
    }
    free(array->hydro.parts);
    array->hydro.parts = parts;
  }

  array->hydro.n = new_n_part;

  /* Dark matter */
  /* Check if need to free memory */
  if (new_n_gpart == 0 && array->grav.n != 0) {
    free(array->grav.parts);
    array->grav.parts = NULL;
  }
  /* Check if need to change the size */
  else if (new_n_gpart != array->grav.n) {

    array->grav.allocated_size = new_n_gpart;

    /* Allocate the new array */
    struct logger_gparticle *parts = (struct logger_gparticle *)malloc(
        new_n_gpart * sizeof(struct logger_gparticle));
    if (parts == NULL) {
      error("Failed to allocate the dark matter particles.");
    }

    /* Copy the previous particles */
    if (array->grav.n != 0) {
      memcpy(parts, array->grav.parts,
             array->grav.n * sizeof(struct logger_gparticle));
    }
    free(array->grav.parts);
    array->grav.parts = parts;
  }

  array->grav.n = new_n_gpart;

  /* Stars */
  /* Check if need to free memory */
  if (new_n_spart == 0 && array->stars.n != 0) {
    free(array->stars.parts);
    array->stars.parts = NULL;
  }
  /* Check if need to change the size */
  else if (new_n_spart != array->stars.n) {

    array->stars.allocated_size = new_n_spart;

    /* Allocate the new array */
    struct logger_sparticle *parts = (struct logger_sparticle *)malloc(
        new_n_spart * sizeof(struct logger_sparticle));
    if (parts == NULL) {
      error("Failed to allocate the stars particles.");
    }

    /* Copy the previous particles */
    if (array->stars.n != 0) {
      memcpy(parts, array->stars.parts,
             array->stars.n * sizeof(struct logger_sparticle));
    }
    free(array->stars.parts);
    array->stars.parts = parts;
  }

  array->stars.n = new_n_spart;
}

/**
 * @brief Add an hydro particle (save only the offset).
 *
 * @param array The #logger_particle_array.
 * @param offset The offset of the new particle in the logfile.
 */
void logger_particle_array_add_hydro(struct logger_particle_array *array,
                                     size_t offset) {
  /* Save the offset */
  array->hydro.parts[array->hydro.n].offset = offset;

  /* Update the number of particles */
  array->hydro.n++;

  /* Increase the size if required */
  if (array->hydro.n == array->hydro.allocated_size) {
    logger_particle_array_change_size(array, 2 * array->hydro.allocated_size,
                                      array->grav.allocated_size,
                                      array->stars.allocated_size);
  }
}

/**
 * @brief Add a star (save only the offset).
 *
 * @param array The #logger_particle_array.
 * @param offset The offset of the new particle in the logfile.
 */
void logger_particle_array_add_stars(struct logger_particle_array *array,
                                     size_t offset) {
  /* Save the offset */
  array->stars.parts[array->stars.n].offset = offset;

  /* Update the number of particles */
  array->stars.n++;

  /* Increase the size if required */
  if (array->stars.n == array->stars.allocated_size) {
    logger_particle_array_change_size(array, array->hydro.allocated_size,
                                      array->grav.allocated_size,
                                      2 * array->stars.allocated_size);
  }
}

/**
 * @brief Add a gravity particle (save only the offset).
 *
 * @param array The #logger_particle_array.
 * @param offset The offset of the new particle in the logfile.
 */
void logger_particle_array_add_gravity(struct logger_particle_array *array,
                                       size_t offset) {

  /* Save the offset */
  array->grav.parts[array->grav.n].offset = offset;

  /* Update the number of particles */
  array->grav.n++;

  /* Increase the size if required */
  if (array->grav.n == array->grav.allocated_size) {
    logger_particle_array_change_size(array, array->hydro.allocated_size,
                                      2 * array->grav.allocated_size,
                                      array->stars.allocated_size);
  }
}

/**
 * @brief Remove the flagged particles from the arrays and add the particles
 * from the dynamic array.
 *
 * @param prev The array containing the particles before the required time.
 * @param next The array containing the particles after the required time.
 * @param new The new particles.
 * @param n_deleted_hydro The number of flagged particles (hydro).
 * @param n_deleted_grav The number of flagged particles (gravity).
 * @param n_deleted_stars The number of flagged particles (stars).
 */
void logger_particle_array_update(struct logger_particle_array *prev,
                                  struct logger_particle_array *next,
                                  struct logger_particle_array *tmp,
                                  size_t n_deleted_hydro, size_t n_deleted_grav,
                                  size_t n_deleted_stars) {

  if (prev->hydro.n != next->hydro.n) {
    error("The previous and next arrays are not compatible.");
  }
  if (prev->grav.n != next->grav.n) {
    error("The previous and next arrays are not compatible.");
  }
  if (prev->stars.n != next->stars.n) {
    error("The previous and next arrays are not compatible.");
  }

  /* Number of particles without the deleted particles */
  size_t n_hydro = prev->hydro.n;
  size_t n_grav = prev->grav.n;
  size_t n_stars = prev->stars.n;

  /* Remove the hydro particles */
  for (size_t i = 0; i < n_hydro && n_deleted_hydro != 0; i++) {
    if (prev->hydro.parts[i].flag == logger_flag_delete) {
      /* Ensure that we do not have a deleted particle
         at the end of the array. */
      while (prev->hydro.parts[prev->hydro.n - 1].flag == logger_flag_delete) {
        n_hydro--;
        n_deleted_hydro--;
      }

      /* Copy the particle from the end to the current one. */
      prev->hydro.parts[i] = prev->hydro.parts[n_hydro - 1];
      next->hydro.parts[i] = next->hydro.parts[n_hydro - 1];
      n_hydro--;
      n_deleted_hydro--;
    }
  }

  /* Remove the gravity particles */
  for (size_t i = 0; i < n_grav && n_deleted_grav != 0; i++) {
    if (prev->grav.parts[i].flag == logger_flag_delete) {
      /* Ensure that we do not have a deleted particle
         at the end of the array. */
      while (prev->grav.parts[prev->grav.n - 1].flag == logger_flag_delete) {
        n_grav--;
        n_deleted_grav--;
      }

      /* Copy the particle from the end to the current one. */
      prev->grav.parts[i] = prev->grav.parts[n_grav - 1];
      next->grav.parts[i] = next->grav.parts[n_grav - 1];
      n_grav--;
      n_deleted_grav--;
    }
  }

  /* Remove the stars particles */
  for (size_t i = 0; i < n_stars && n_deleted_stars != 0; i++) {
    if (prev->stars.parts[i].flag == logger_flag_delete) {
      /* Ensure that we do not have a deleted particle
         at the end of the array. */
      while (prev->stars.parts[prev->stars.n - 1].flag == logger_flag_delete) {
        n_stars--;
        n_deleted_stars--;
      }

      /* Copy the particle from the end to the current one. */
      prev->stars.parts[i] = prev->stars.parts[n_stars - 1];
      next->stars.parts[i] = next->stars.parts[n_stars - 1];
      n_stars--;
      n_deleted_stars--;
    }
  }

  /* Check that we have deleted all the particles */
  if (n_deleted_hydro != 0) {
    error("Something went wrong during the deletion of hydro particles.");
  }

  if (n_deleted_grav != 0) {
    error("Something went wrong during the deletion of gravity particles.");
  }

  if (n_deleted_stars != 0) {
    error("Something went wrong during the deletion of stars particles.");
  }

  /* Number of particles at the end of this function */
  const size_t new_n_hydro = n_hydro + tmp->hydro.n;
  const size_t new_n_grav = n_grav + tmp->grav.n;
  const size_t new_n_stars = n_stars + tmp->stars.n;

  /* Reallocate the memory */
  logger_particle_array_change_size(prev, new_n_hydro, new_n_grav, new_n_stars);
  logger_particle_array_change_size(next, new_n_hydro, new_n_grav, new_n_stars);

  /* Copy the new particles at the end of the arrays */
  if (tmp->hydro.n != 0) {
    memcpy(prev->hydro.parts + n_hydro, tmp->hydro.parts,
           tmp->hydro.n * sizeof(struct logger_particle));
  }
  if (tmp->grav.n != 0) {
    memcpy(prev->grav.parts + n_grav, tmp->grav.parts,
           tmp->grav.n * sizeof(struct logger_gparticle));
  }
  if (tmp->stars.n != 0) {
    memcpy(prev->stars.parts + n_stars, tmp->stars.parts,
           tmp->stars.n * sizeof(struct logger_sparticle));
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we removed all the required particles */
  for (size_t i = 0; i < prev->hydro.n; i++) {
    if (prev->hydro.parts[i].flag == logger_flag_delete) {
      error("Failed to delete a particle.");
    }
    /* Only the offset is correctly set for the next particles */
    if (i > n_hydro && next->hydro.parts[i].flag == logger_flag_delete) {
      error("Failed to delete a particle.");
    }
  }

  /* Check that we removed all the required particles */
  for (size_t i = 0; i < prev->grav.n; i++) {
    if (prev->grav.parts[i].flag == logger_flag_delete) {
      error("Failed to delete a particle.");
    }
    /* Only the offset is correctly set for the next particles */
    if (i > n_grav && next->grav.parts[i].flag == logger_flag_delete) {
      error("Failed to delete a particle.");
    }
  }

  /* Check that we removed all the required particles */
  for (size_t i = 0; i < prev->stars.n; i++) {
    if (prev->stars.parts[i].flag == logger_flag_delete) {
      error("Failed to delete a particle.");
    }
    /* Only the offset is correctly set for the next particles */
    if (i > n_stars && next->stars.parts[i].flag == logger_flag_delete) {
      error("Failed to delete a particle.");
    }
  }
#endif
}
