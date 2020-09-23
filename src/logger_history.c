/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "logger_history.h"

/* Standard includes */
#include <string.h>

/* Local include */
#include "logger_io.h"
#include "part.h"

#if defined(WITH_LOGGER)

#define LOGGER_HISTORY_INIT_SIZE 1024

/**
 * @brief Initialize the structure for the first time.
 *
 * @param hist The #logger_history.
 */
void logger_history_first_init(struct logger_history *hist) {

  for (int i = 0; i < swift_type_count; i++) {
    /* Set the counters to their initial value */
    hist->size[i] = 0;
    hist->capacity[i] = LOGGER_HISTORY_INIT_SIZE;

    hist->data[i] = (struct logger_index_data *)malloc(
        sizeof(struct logger_index_data) * LOGGER_HISTORY_INIT_SIZE);
    if (hist->data[i] == NULL) {
      error(
          "Failed to allocate memory for the logger_history (particle type %i)",
          i);
    }
  }
}

/**
 * @brief Initialize the structure (for example just after a dump).
 *
 * @param hist The #logger_history.
 * @param params The #swift_params.
 * @param already_allocated Are the data already allocated? (Need to free it?)
 */
void logger_history_init(struct logger_history *hist) {

  for (int i = 0; i < swift_type_count; i++) {
    /* Set the counters to their initial value */
    hist->size[i] = 0;
    hist->capacity[i] = LOGGER_HISTORY_INIT_SIZE;

    /* Allocate the initial size */
    free(hist->data[i]);

    hist->data[i] = (struct logger_index_data *)malloc(
        sizeof(struct logger_index_data) * LOGGER_HISTORY_INIT_SIZE);
    if (hist->data[i] == NULL) {
      error(
          "Failed to allocate memory for the logger_history (particle type %i)",
          i);
    }
  }
}

/**
 * @brief Clean the structure (e.g. just before exiting).
 *
 * @param hist The #logger_history.
 */
void logger_history_clean(struct logger_history *hist) {
  for (int i = 0; i < swift_type_count; i++) {
    /* Set the counters to 0 */
    hist->size[i] = 0;
    hist->capacity[i] = 0;

    /* Free the memory */
    if (hist->data[i] != NULL) {
      free(hist->data[i]);
      hist->data[i] = NULL;
    }
  }
}

/**
 * @brief Log a the particle information into the #logger_history.
 *
 * @param hist The #logger_history.
 * @param data The data from the particle.
 * @param part_type The particle type.
 */
void logger_history_log(struct logger_history *hist,
                        struct logger_index_data data,
                        enum part_type part_type) {

  /* Check if enough space is left */
  if (hist->size[part_type] + 1 >= hist->capacity[part_type]) {
    /* Compute the previous amount of memory */
    const size_t memsize =
        sizeof(struct logger_index_data) * hist->capacity[part_type];

    /* Increase the capacity of the array */
    hist->capacity[part_type] *= 2;

    /* Allocate the new array and copy the content of the previous one */
    struct logger_index_data *tmp =
        (struct logger_index_data *)malloc(2 * memsize);

    memcpy(tmp, hist->data[part_type], memsize);

    /* Free the previous array and switch the pointers */
    free(hist->data[part_type]);
    hist->data[part_type] = tmp;
  }

  /* Save the new particle */
  hist->data[part_type][hist->size[part_type]] = data;

  /* Increase the element counter */
  hist->size[part_type] += 1;
}

/**
 * @brief Log a particle in the history.
 *
 * @param hist The #logger_history.
 * @param p The #part to log.
 * @param xp The #xpart to log.
 */
void logger_history_log_part(struct logger_history *hist, const struct part *p,
                             const struct xpart *xp) {
  /* Save the required data */
  struct logger_index_data data;
  data.id = p->id;
  data.offset = xp->logger_data.last_offset;

  /* Log the data */
  logger_history_log(hist, data, swift_type_gas);
}

/**
 * @brief Log a particle in the history.
 *
 * @param hist The #logger_history.
 * @param sp The #spart to log.
 */
void logger_history_log_spart(struct logger_history *hist,
                              const struct spart *sp) {
  /* Save the required data */
  struct logger_index_data data;
  data.id = sp->id;
  data.offset = sp->logger_data.last_offset;

  /* Log the data */
  logger_history_log(hist, data, swift_type_stars);
}

/**
 * @brief Log a particle in the history.
 *
 * @param hist The #logger_history.
 * @param gp The #gpart to log.
 */
void logger_history_log_gpart(struct logger_history *hist,
                              const struct gpart *gp) {
  /* Save the required data */
  struct logger_index_data data;
  data.id = gp->id_or_neg_offset;
  data.offset = gp->logger_data.last_offset;

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->id_or_neg_offset < 0) {
    error("Trying to log a non gravity particle");
  }
#endif

  /* Log the data */
  logger_history_log(hist, data, swift_type_dark_matter);
}

/**
 * @brief Log a particle in the history.
 *
 * @param hist The #logger_history.
 * @param bp The #bpart to log.
 */
void logger_history_log_bpart(struct logger_history *hist,
                              const struct bpart *bp) {
  error("TODO");
}

/**
 * @brief Write the history into an index file.
 *
 * @param hist The #logger_history.
 * @param e The #engine.
 * @param f The file where to write the history.
 */
void logger_history_write(struct logger_history *hist, struct engine *e,
                          FILE *f) {
  /* Write the number of particles */
  fwrite(hist->size, sizeof(uint64_t), swift_type_count, f);

  /* write the particles */
  for (int i = 0; i < swift_type_count; i++) {
    /* Generate the structures for writing the index file */
    const int num_fields = 2;
    struct io_props list[2];
    list[0] =
        io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS,
                             0.f, hist->data[i], id, "Field not used");
    list[1] = io_make_output_field("Offset", UINT64, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   hist->data[i], offset, "Field not used");

    write_index_array(e, f, list, num_fields, hist->size[i]);
  }

  /* Reset the logger history */
  logger_history_init(hist);
}

/**
 * @brief Compute the number of particles in the history.
 *
 * @param hist The #logger_history.
 *
 * @return The total number of particle stored.
 */
size_t logger_history_get_size(const struct logger_history *hist) {
  size_t number = 0;
  for (int i = 0; i < swift_type_count; i++) {
    number += hist->size[i];
  }

  return number;
}

void logger_history_dump(const struct logger_history *hist, FILE *stream) {
  restart_write_blocks((void *)hist, sizeof(struct logger_history), 1, stream,
                       "logger_history", "logger_history");

  for (int i = 0; i < swift_type_count; i++) {
    if (hist->size[i] != 0)
      restart_write_blocks(
          (void *)hist->data[i], sizeof(struct logger_index_data),
          hist->size[i], stream, "logger_history_data", "logger_history_data");
  }
}

void logger_history_restore(struct logger_history *hist, FILE *stream) {
  restart_read_blocks((void *)hist, sizeof(struct logger_history), 1, stream,
                      NULL, "logger_history");

  for (int i = 0; i < swift_type_count; i++) {
    hist->data[i] =
        malloc(hist->capacity[i] * sizeof(struct logger_index_data));
    if (hist->data[i] == NULL) {
      error("Failed to allocate array for logger history");
    }

    if (hist->size[i] != 0)
      restart_read_blocks((void *)hist->data[i],
                          sizeof(struct logger_index_data), hist->size[i],
                          stream, NULL, "logger_history_data");
  }
}

#endif  // WITH_LOGGER
