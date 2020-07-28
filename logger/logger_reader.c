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

/* Include corresponding header */
#include "logger_reader.h"

/* Include standard library */
#include <sys/sysinfo.h>
#include <unistd.h>

/* Include local headers */
#include "threadpool.h"

#define nr_threads get_nprocs()

/**
 * @brief Initialize the reader.
 *
 * @param reader The #logger_reader.
 * @param basename The basename of the logger files.
 * @param verbose The verbose level.
 */
void logger_reader_init(struct logger_reader *reader, const char *basename,
                        int verbose) {
  if (verbose > 1) message("Initializing the reader.");

  /* Set the variable to the default values */
  reader->time.time = -1.;
  reader->time.int_time = 0;
  reader->time.time_offset = 0;

  /* Copy the base name */
  strcpy(reader->basename, basename);

  /* Initialize the reader variables. */
  reader->verbose = verbose;

  /* Generate the logfile filename */
  char logfile_name[STRING_SIZE];
  sprintf(logfile_name, "%s.dump", basename);

  /* Initialize the log file. */
  logger_logfile_init_from_file(&reader->log, logfile_name, reader,
                                /* only_header */ 0);

  /* Initialize the index files */
  logger_reader_init_index(reader);

  if (verbose > 1) message("Initialization done.");
}

/**
 * @brief Initialize the index part of the reader.
 *
 * @param reader The #logger_reader.
 */
void logger_reader_init_index(struct logger_reader *reader) {
  /* Initialize the logger_index */
  logger_index_init(&reader->index.index, reader);

  /* Count the number of files */
  int count = 0;
  while (1) {
    char filename[STRING_SIZE + 50];
    sprintf(filename, "%s_%04i.index", reader->basename, count);

    /* Check if file exists */
    if (access(filename, F_OK) != -1) {
      count++;
    } else {
      break;
    }
  }

  reader->index.n_files = count;

  /* Initialize the arrays */
  if ((reader->index.times = (double *)malloc(count * sizeof(double))) == NULL) {
    error("Failed to allocate the list of times");
  }
  if ((reader->index.int_times = (integertime_t *)malloc(count * sizeof(integertime_t))) == NULL) {
    error("Failed to allocate the list of times");
  }

  /* Get the information contained in the headers */
  for (int i = 0; i < reader->index.n_files; i++) {
    char filename[STRING_SIZE + 50];
    sprintf(filename, "%s_%04i.index", reader->basename, i);

    /* Read the header */
    logger_index_read_header(&reader->index.index, filename);

    /* Save the required information */
    reader->index.times[i] = reader->index.index.time;
    reader->index.int_times[i] = reader->index.index.integer_time;
  }
}

/**
 * @brief Free the reader.
 *
 * @param reader The #logger_reader.
 */
void logger_reader_free(struct logger_reader *reader) {
  /* Free the log. */
  logger_logfile_free(&reader->log);

  if (reader->time.time != -1.) {
    logger_index_free(&reader->index.index);
  }

  free(reader->index.int_times);
  free(reader->index.times);
}

/**
 * @brief Set the reader to a given time and read the correct index file.
 *
 * @param reader The #logger_reader.
 * @param time The requested time.
 */
void logger_reader_set_time(struct logger_reader *reader, double time) {
  /* Set the time */
  reader->time.time = time;

  /* Find the correct index */
  unsigned int left = 0;
  unsigned int right = reader->index.n_files - 1;

  while (left != right) {
    /* Do a ceil - division */
    unsigned int m = (left + right + 1) / 2;
    if (reader->index.times[m] > time) {
      right = m - 1;
    } else {
      left = m;
    }
  }

  /* Generate the filename */
  char filename[STRING_SIZE + 50];
  sprintf(filename, "%s_%04u.index", reader->basename, left);

  /* Check if the file is already mapped */
  if (reader->index.index.index.map != NULL) {
    logger_index_free(&reader->index.index);
  }

  /* Read the file */
  logger_index_read_header(&reader->index.index, filename);
  logger_index_map_file(&reader->index.index, filename, /* sorted */ 1);

  /* Get the offset of the time chunk */
  size_t ind = time_array_get_index_from_time(&reader->log.times, time);
  /* ind == 0 and ind == 1 are the same time, but when reading we need
     data before the initial time. */
  if (ind == 0) {
    ind = 1;
  }

  /* Check if we requested exactly a time step  */
  if (reader->log.times.records[ind].time != time) {
    /* In order to interpolate, we need to be above and not below the time */
    ind += 1;
  }

  /* Save the values */
  reader->time.index = ind;
  reader->time.int_time = reader->log.times.records[ind].int_time;
  reader->time.time_offset = reader->log.times.records[ind].offset;
}

/**
 * @brief Provides the number of particle (per type) from the index file.
 *
 * @param reader The #logger_reader.
 * @param n_type (output) The number of particle type possible.
 *
 * @return For each type possible, the number of particle.
 */
const uint64_t *logger_reader_get_number_particles(struct logger_reader *reader,
                                                   int *n_type) {
  *n_type = swift_type_count;
  return reader->index.index.nparts;
}

/**
 * @brief Read a single gpart from the logger and interpolate it if needed.
 *
 * The temporary arrays are provided to this function in order to allocate
 * them only once for all the particles.
 *
 * @param reader The #logger_reader.
 * @param time The requested time.
 * @param offset_time The offset of the corresponding time record.
 * @param interp_type The type of interpolation requested.
 * @param part_type The particle type.
 * @param offset_last_full_record The offset of this particle last record
 containing all the fields.
 * @param fields_wanted The fields wanted
 (sorted according to the particle's type reading order).
 * @param n_fields_wanted Size of fields_wanted.
 * @param output Array of pointers to the output array
 (sorted according to the particle's type reading order).
 * @param tmp_output Temporary array for storing the particle's fields.
 * @param time_before Temporary array for storing the time of each field read.
 * @param field_found Temporary array for storing the fields found when
 looking the the records after the requested time.
 */
void logger_reader_read_single_particle(
    struct logger_reader *reader, double time, size_t offset_time,
    enum logger_reader_type interp_type, enum part_type part_type,
    const size_t offset_last_full_record,
    const int *fields_wanted, const int n_fields_wanted, void **output,
    void **tmp_output, double *time_before, int *fields_found) {

  const struct header *h = &reader->log.header;
  size_t offset = offset_last_full_record;

  /* Get the mask_id from the particle type */
  const int *logger_mask_id = NULL;
  switch(part_type) {
    case swift_type_gas:
      logger_mask_id = hydro_logger_mask_id;
      break;
    case swift_type_dark_matter:
      logger_mask_id = gravity_logger_mask_id;
      break;
    case swift_type_stars:
      logger_mask_id = stars_logger_mask_id;
      break;
  default:
      error("Type not implemented yet.");
  }

  /* Find the data for the previous record.
     As we start from a full record,
     no need to check if all the fields are found.
  */
  while(offset < offset_time) {
    /* Read the particle. */
    size_t mask, h_offset;
    switch (part_type) {
      case swift_type_gas:
        logger_particle_read(reader, offset, fields_wanted, n_fields_wanted, tmp_output,
                             &mask, &h_offset);
        break;
      case swift_type_dark_matter:
        logger_gparticle_read(reader, offset, fields_wanted, n_fields_wanted, tmp_output,
                              &mask, &h_offset);
        break;
      case swift_type_stars:
        logger_sparticle_read(reader, offset, fields_wanted, n_fields_wanted, tmp_output,
                              &mask, &h_offset);
        break;
      default:
        error("Type not implemented yet.");
    }

    /* Get the time. */
    double current_time = time_array_get_time(&reader->log.times, offset);

    /* Copy into the output array. */
    for(int i = 0; i < n_fields_wanted; i++) {
      const int field_id = logger_mask_id[fields_wanted[i]];
      /* Check if the field is present in this record. */
      if (mask & h->masks[field_id].mask) {
        time_before[i] = current_time;
        memcpy(output[i], tmp_output[i], h->masks[field_id].size);
      }
    }

    /* Go to the next record. */
    offset += h_offset;
  }

  /* When interpolating, we need to get the next record after
     the requested time.
  */
  if (interp_type != logger_reader_const) {

    /* Initialize the variables that record the fields found. */
    int number_field_still_to_recover = n_fields_wanted;
    for(int i = 0; i < n_fields_wanted; i++) {
      fields_found[i] = 0;
    }

    /* Loop over the records until having all the fields. */
    while (number_field_still_to_recover != 0) {

      /* Read the particle. */
      size_t mask, h_offset;
      logger_gparticle_read(reader, offset, fields_wanted, n_fields_wanted, tmp_output,
                            &mask, &h_offset);

      /* Get the time of the record. */
      double current_time = time_array_get_time(&reader->log.times, offset);

      /* Copy the data into the output array. */
      for(int i = 0; i < n_fields_wanted; i++) {
        const int field_id = logger_mask_id[fields_wanted[i]];
        /* Check if the mask is present in this record and still not found. */
        if (!fields_found[i] && mask & h->masks[field_id].mask) {

          /* Mark the field as being found. */
          number_field_still_to_recover -= 1;
          fields_found[i] = 1;

          /* If the times are the same, no need to interpolate
             (e.g. when requesting a logger_log_all_particles). */
          if (current_time == time_before[i])
            continue;

          /* Interpolate the data. */
          switch(part_type) {
            case swift_type_gas:
              hydro_logger_interpolate_field(
                  tmp_output[i], output[i], output[i], time_before[i],
                  current_time, time, fields_wanted[i]);
              break;
            case swift_type_dark_matter:
              gravity_logger_interpolate_field(
                  tmp_output[i], output[i], output[i], time_before[i],
                  current_time, time, fields_wanted[i]);
              break;
            case swift_type_stars:
              stars_logger_interpolate_field(
                  tmp_output[i], output[i], output[i], time_before[i],
                  current_time, time, fields_wanted[i]);
            break;
            default:
              error("Type not implemented yet.");
          }
        }
      }
    }
  }
}

/**
 * @brief Sort the fields according to the reading order of a given particle type.
 *
 * @param fields_wanted The fields to sort.
 * @param sorted_indices (out) The indices of fields_wanted in the sorted order (need to be allocated).
 * @param sorted_fields_wanted (out) fields_wanted in the sorted order (need to be allocated).
 * @param n_fields_wanted Size of fields_wanted, sorted_indices and sorted_fields_wanted.
 * @param type The type of the particle.
 */
void logger_reader_sort_ids(const int *fields_wanted, int *sorted_indices,
                            int *sorted_fields_wanted,
                            const int n_fields_wanted, enum part_type type) {

  /* Get the correct variables. */
  int n_max = 0;
  int *fields_ids = NULL;
  switch(type) {
    case swift_type_gas:
      n_max = hydro_logger_field_count;
      fields_ids = hydro_logger_mask_id;
      break;
    case swift_type_dark_matter:
      n_max = gravity_logger_field_count;
      fields_ids = gravity_logger_mask_id;
      break;
    case swift_type_stars:
      n_max = stars_logger_field_count;
      fields_ids = stars_logger_mask_id;
      break;
    default:
      error("Particle type not implemented yet.");
  }

  /* No need to have an efficient sort here (max ~10 items). */

  int current = 0;
  /* Go over all the existing fields in the correct order. */
  for(int j = 0; j < n_max; j++) {
    int field_id = fields_ids[j];
    /* Find if the fields is requested. */
    for(int i = 0; i < n_fields_wanted; i++) {
      /* The field is found, go to the next one. */
      if (field_id == fields_wanted[i]) {
        sorted_indices[current] = i;
        sorted_fields_wanted[current] = j;
        current += 1;
        break;
      }
    }
  }

  /* Check if we found all the fields. */
  if (current != n_fields_wanted) {
    error("Failed to find a field for %s", part_type_names[type]);
  }
}

/**
 * @brief Read all the particles of a given type.
 *
 * @param part_type The name of the particle type following part_type.h (e.g. gas, dark_matter, stars).
 * @param type The type of functions (e.g. hydro, gravity, stars).
 */
#define READ_ALL_PARTICLES(part_type, type)     \
  ({                                                                    \
  /* Do the hydro. */                                                   \
  struct index_data *data = logger_index_get_data(&reader->index.index, swift_type_##part_type); \
                                                                        \
  /* Sort the fields in order to read the correct bits. */              \
  logger_reader_sort_ids(fields_wanted, sorted_indices,                 \
                         sorted_fields_wanted,                          \
                         n_fields_wanted, swift_type_##part_type);      \
                                                                        \
  /* Allocate the memory for the temporary output. */                   \
  for(int i = 0; i < n_fields_wanted; i++) {                            \
    const int field_id = type##_logger_mask_id[sorted_fields_wanted[i]]; \
    if ((tmp_output_single[i] = malloc(h->masks[field_id].size)) == NULL) { \
      error("Failed to allocate the temporary output.");                \
    }                                                                   \
  }                                                                     \
                                                                        \
  /* Read the dark matter particles */                                  \
  for (size_t i = 0; i < n_part[swift_type_##part_type]; i++) {         \
    /* Get the offset */                                                \
    size_t offset = data[i].offset;                                     \
                                                                        \
    /* Sort the output into output_single. */                           \
    for(int field = 0; field < n_fields_wanted; field++) {              \
      const int field_id = fields_wanted[sorted_indices[field]];        \
      output_single[field] = output[sorted_indices[field]] + i * h->masks[field_id].size; \
    }                                                                   \
                                                                        \
    /* Read the particle. */                                            \
    logger_reader_read_single_particle(                                 \
                                       reader, time, reader->time.time_offset, interp_type, swift_type_##part_type, offset, \
                                       sorted_fields_wanted, n_fields_wanted, output_single, tmp_output_single, \
                                       time_output, fields_found);      \
  }                                                                     \
                                                                        \
  /* Free the allocated memory for this type of particle. */            \
  for(int i = 0; i < n_fields_wanted; i++) {                            \
    free(tmp_output_single[i]);                                         \
  }                                                                     \
})

/**
 * @brief Read all the particles from the index file.
 *
 * @param reader The #logger_reader.
 * @param time The requested time for the particle.
 * @param interp_type The type of interpolation.
 * @param fields_wanted The fields requested (index of the header->masks).
 * @param n_fields_wanted Number of field requested.
 * @param output Pointer to the output array. Size: (n_fields_wanted, sum(n_part)).
 * @param n_part Number of particles of each type.
 */
void logger_reader_read_all_particles(struct logger_reader *reader, double time,
                                      enum logger_reader_type interp_type,
                                      const int *fields_wanted, const int n_fields_wanted,
                                      void **output, const uint64_t *n_part) {

  const struct header *h = &reader->log.header;

  /* Allocate temporary memory. */
  /* The index of fields_wanted sorted according to the fields order. */
  int *sorted_indices = (int *) malloc(sizeof(int) * n_fields_wanted);
  if (sorted_indices == NULL) {
    error("Failed to allocate the array of sorted indices.");
  }
  /* fields_wanted sorted according to the fields order. */
  int *sorted_fields_wanted = (int *) malloc(sizeof(int) * n_fields_wanted);
  if (sorted_fields_wanted == NULL) {
    error("Failed to allocate the array of sorted fields.");
  }
  /* Pointer to the output array in the correct reading order. */
  void **output_single = malloc(sizeof(void*) * n_fields_wanted);
  if (output_single == NULL) {
    error("Failed to allocate the output.");
  }
  /* Temporary array for storing a single particle read. */
  void **tmp_output_single = malloc(sizeof(void*) * n_fields_wanted);
  if (tmp_output_single == NULL) {
    error("Failed to allocate the temporary array for a single particle.");
  }

  /* Output time of the fields (used internally by the single particle readers). */
  double *time_output = (double *) malloc(sizeof(double) * n_fields_wanted);
  if (time_output == NULL) {
    error("Failed to allocate the array of the times.");
  }
  /* Fields found when moving after the requested time (used internally by the single particle readers) */
  int *fields_found = (int *) malloc(sizeof(int) * n_fields_wanted);
  if (fields_found == NULL) {
    error("Failed to allocate the array of fields found.");
  }

  /* Do the hydro. */
  if (n_part[swift_type_gas] != 0)
    READ_ALL_PARTICLES(gas, hydro);

  /* Do the dark matter now. */
  if (n_part[swift_type_dark_matter] != 0)
    READ_ALL_PARTICLES(dark_matter, gravity);

  /* Do the stars now. */
  if (n_part[swift_type_stars] != 0)
    READ_ALL_PARTICLES(stars, stars);

  /* Free the memory. */
  free(output_single);
  free(tmp_output_single);
  free(sorted_indices);
  free(sorted_fields_wanted);
  free(time_output);
  free(fields_found);
}

/**
 * @brief Get the simulation initial time.
 *
 * @param reader The #logger_reader.
 *
 * @return The initial time
 */
double logger_reader_get_time_begin(struct logger_reader *reader) {
  return reader->log.times.records[0].time;
}

/**
 * @brief Get the simulation final time.
 *
 * @param reader The #logger_reader.
 *
 * @return The final time
 */
double logger_reader_get_time_end(struct logger_reader *reader) {
  const size_t ind = reader->log.times.size;
  return reader->log.times.records[ind - 1].time;
}

/**
 * @brief Get the offset of the last timestamp before a given time.
 *
 * @param reader The #logger_reader.
 * @param time The requested time.
 *
 * @return The offset of the timestamp.
 */
size_t logger_reader_get_next_offset_from_time(struct logger_reader *reader,
                                               double time) {
  size_t ind = time_array_get_index_from_time(&reader->log.times, time);
  /* We do not want to have the sentiel */
  if (reader->log.times.size - 2 == ind) {
    ind -= 1;
  }
  return reader->log.times.records[ind + 1].offset;
}

/**
 * @brief Read a record without knowing if it is a particle or a timestamp.
 *
 * WARNING This function asssumes that all the particles are hydro particles.
 * Thus it should be used only for testing the code.
 *
 * @param reader The #logger_reader.
 * @param output The already allocated buffer containing all the fields possible for
 * an hydro particle. (out) The particle if the record is a particle
 * @param time (out) The time if the record is a timestamp.
 * @param is_particle (out) 1 if the record is a particle 0 otherwise.
 * @param offset The offset of the record to read.
 *
 * @return The offset after the record.
 */
size_t logger_reader_read_record(struct logger_reader *reader, void **output,
                                 double *time, int *is_particle, size_t offset) {

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  size_t mask = 0;
  size_t h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, (char *)map + offset, &mask, &h_offset);

  *is_particle = !(mask & h->timestamp_mask);
  /* The record is a particle. */
  if (*is_particle) {
    /* We request all the fields */
    int *required_fields = malloc(hydro_logger_field_count * sizeof(int));
    if (required_fields == NULL) {
      error("Failed to allocate the required fields.");
    }
    for(int i = 0; i < hydro_logger_field_count; i++) {
      required_fields[i] = i;
    }

    offset = logger_particle_read(reader, offset, required_fields, hydro_logger_field_count,
                                  output, &mask, &h_offset);

    free(required_fields);
  }
  /* The record is a timestamp. */
  else {
    integertime_t not_used = 0;
    offset = time_read(&not_used, time, reader, offset);
  }

  return offset;
}
