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
  reader->index.times = (double *)malloc(count * sizeof(double));
  reader->index.int_times =
      (integertime_t *)malloc(count * sizeof(integertime_t));

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
}

/**
 * @brief Read a record (timestamp or particle).
 *
 * WARNING THIS FUNCTION WORKS ONLY WITH HYDRO PARTICLES.
 *
 * @param reader The #logger_reader.
 * @param lp (out) The #logger_particle (if the record is a particle).
 * @param time (out) The time read (if the record is a timestamp).
 * @param is_particle Is the record a particle (or a timestamp)?
 * @param offset The offset in the file.
 *
 * @return The offset after this record.
 */
size_t logger_reader_read_record(struct logger_reader *reader,
                                 struct logger_particle *lp, double *time,
                                 int *is_particle, size_t offset) {

  struct logger_logfile *log = &reader->log;

  /* Read mask to find out if timestamp or particle. */
  size_t mask = 0;
  logger_loader_io_read_mask(&log->header, (char *)log->log.map + offset, &mask,
                             NULL);

  /* Check if timestamp or not. */
  int ind = header_get_field_index(&log->header, "Timestamp");
  if (ind == -1) {
    error("File header does not contain a mask for time.");
  }
  if (log->header.masks[ind].mask == mask) {
    *is_particle = 0;
    integertime_t int_time = 0;
    offset = time_read(&int_time, time, reader, offset);
  } else {
    *is_particle = 1;
    offset =
        logger_particle_read(lp, reader, offset, *time, logger_reader_const);
  }

  return offset;
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
void logger_reader_read_single_gparticle(
    struct logger_reader *reader, double time, size_t offset_time,
    enum logger_reader_type interp_type, const size_t offset_last_full_record,
    const int *fields_wanted, const int n_fields_wanted, void **output,
    void **tmp_output, double *time_before, int *fields_found) {

  const struct header *h = &reader->log.header;
  size_t offset = offset_last_full_record;


  /* Find the data for the previous record.
     As we start from a full record,
     no need to check if all the fields are found.
  */
  while(offset < offset_time) {

    /* Read the particle. */
    size_t mask, h_offset;
    logger_gparticle_read(reader, offset, fields_wanted, n_fields_wanted, tmp_output,
                          &mask, &h_offset);

    /* Get the time. */
    double current_time = time_array_get_time(&reader->log.times, offset);

    /* Copy the into the output array. */
    for(int i = 0; i < n_fields_wanted; i++) {
      const int field_id = gravity_logger_mask_id[fields_wanted[i]];
      /* Check if the field is present in this record. */
      if (mask & h->masks[field_id].mask) {
        time_before[i] = current_time;
        memcpy(output[i], tmp_output[i], h->masks[gravity_logger_mask_id[i]].size);
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
        const int field_id = gravity_logger_mask_id[fields_wanted[i]];
        /* Check if the mask is present in this record and still not found. */
        if (!fields_found[i] && mask & h->masks[field_id].mask) {
          /* Interpolate the data. */
          logger_gparticle_interpolate_field(
              tmp_output[i], output[i], output[i], time_before[i],
              current_time, time, fields_wanted[i]);
          /* Mark the field as being found. */
          number_field_still_to_recover -= 1;
          fields_found[i] = 1;
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
    case swift_type_dark_matter:
      n_max = gravity_logger_field_count;
      fields_ids = gravity_logger_mask_id;
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
  /* fields_wanted sorted according to the fields order. */
  int *sorted_fields_wanted = (int *) malloc(sizeof(int) * n_fields_wanted);
  /* Pointer to the output array in the correct reading order. */
  void **output_single = malloc(sizeof(void) * n_fields_wanted);
  /* Temporary array for storing a single particle read. */
  void **tmp_output_single = malloc(sizeof(void) * n_fields_wanted);
  /* Output time of the fields (used internally by the single particle readers). */
  double *time_output = (double *) malloc(sizeof(double) * n_fields_wanted);
  /* Fields found when moving after the requested time (used internally by the single particle readers) */
  int *fields_found = (int *) malloc(sizeof(int) * n_fields_wanted);

  /* Do the dark matter now. */
  struct index_data *data = logger_index_get_data(&reader->index.index, swift_type_dark_matter);

  /* Sort the fields in order to read the correct bits. */
  logger_reader_sort_ids(fields_wanted, sorted_indices,
                         sorted_fields_wanted,
                         n_fields_wanted, swift_type_dark_matter);

  /* Allocate the memory for the temporary output. */
  for(int i = 0; i < n_fields_wanted; i++) {
    const int field_id = gravity_logger_mask_id[sorted_fields_wanted[i]];
    tmp_output_single[i] = malloc(h->masks[field_id].size);
  }

  /* Read the dark matter particles */
  for (size_t i = 0; i < n_part[swift_type_dark_matter]; i++) {
    /* Get the offset */
    size_t offset = data[i].offset;

    /* Sort the output into output_single. */
    for(int field = 0; field < n_fields_wanted; field++) {
      const int field_id = fields_wanted[sorted_indices[field]];
      output_single[field] = output[sorted_indices[field]] + i * h->masks[field_id].size;
    }

    /* Read the particle. */
    logger_reader_read_single_gparticle(
      reader, time, reader->time.time_offset, interp_type, offset,
      sorted_fields_wanted, n_fields_wanted, output_single, tmp_output_single,
      time_output, fields_found);
  }

  /* Free the allocated memory for this type of particle. */
  for(int i = 0; i < n_fields_wanted; i++) {
    free(tmp_output_single[i]);
  }

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
