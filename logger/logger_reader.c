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
 * @brief Read a record (timestamp or particle)
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
size_t reader_read_record(struct logger_reader *reader,
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

  message("Found time");
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
 * @brief Read all the particles from the index file.
 *
 * @param reader The #logger_reader.
 * @param time The requested time for the particle.
 * @param interp_type The type of interpolation.
 * @param parts The array of particles to use.
 * @param n_tot The total number of particles
 */
void logger_reader_read_all_particles(struct logger_reader *reader, double time,
                                      enum logger_reader_type interp_type,
                                      struct logger_particle_array *array,
                                      size_t n_tot) {

  /* Shortcut to some structures */
  struct logger_index *index = &reader->index.index;
  const int spec_flag_ind =
    header_get_field_index(&reader->log.header, "special flags");

  /* Get the correct index file */
  logger_reader_set_time(reader, time);

  /* Read the hydro */
  struct index_data *data = logger_index_get_data(index, /* Particle type */ 0);

  /* Read the particles */
  for (size_t i = 0; i < array->hydro.n; i++) {
    /* Get the offset */
    size_t prev_offset = data[i].offset;
    size_t next_offset = prev_offset;

#ifdef SWIFT_DEBUG_CHECKS
    /* check with the offset of the next timestamp.
     * (the sentinel protects against overflow)
     */
    const size_t ind = reader->time.index + 1;
    if (prev_offset >= reader->log.times.records[ind].offset) {
      error("An offset is out of range (%zi > %zi).", prev_offset,
            reader->log.times.records[ind].offset);
    }
#endif

    while (next_offset < reader->time.time_offset) {
      // TODO use a single call for the offset and mask
      size_t mask = 0;
      logger_loader_io_read_mask(
        &reader->log.header, reader->log.log.map + prev_offset,
        &mask, /* diff_offset */ NULL);

      if (mask & reader->log.header.masks[spec_flag_ind].mask) {
        error("TODO");
      }

      int test = tools_get_next_record(&reader->log.header, reader->log.log.map,
                                       &next_offset, reader->log.log.mmap_size);

      if (test == -1) {
        mask = 0;
        logger_loader_io_read_mask(&reader->log.header,
                                   (char *)reader->log.log.map + prev_offset,
                                   &mask, &next_offset);
        error(
            "Trying to get a particle without next record (mask: %zi, diff "
            "offset: %zi)",
            mask, next_offset);
      }
    }

    /* Read the particle */
    logger_particle_read(&array->hydro.parts[i], reader, prev_offset, reader->time.time,
                         interp_type);
  }
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
 * @brief Get the two particle records around the requested time.
 *
 * @param reader The #logger_reader.
 * @param prev (in) A record before the requested time. (out) The last record
 * before the time.
 * @param next (out) The first record after the requested time.
 * @param time_offset The offset of the requested time.
 */
void logger_reader_move_forward(struct logger_reader *reader,
                                struct logger_particle_array *prev,
                                struct logger_particle_array *next,
                                size_t time_offset) {

  /* Ensure to have enough place in the next */
  logger_particle_array_change_size(next, prev->hydro.n, prev->dark_matter.n, prev->stars.n);

  for(size_t i = 0; i < prev->hydro.n; i++) {
    enum logger_reader_event event = logger_reader_get_next_particle(
      reader, &prev->hydro.parts[i], &next->hydro.parts[i], time_offset);

    if (event != logger_reader_event_null) {
      error("TODO");
    }
  }

  for(size_t i = 0; i < prev->dark_matter.n; i++) {
    enum logger_reader_event event = logger_reader_get_next_gparticle(
      reader, &prev->dark_matter.parts[i], &next->dark_matter.parts[i], time_offset);

    if (event != logger_reader_event_null) {
      error("TODO");
    }
  }

  for(size_t i = 0; i < prev->stars.n; i++) {
    enum logger_reader_event event = logger_reader_get_next_sparticle(
      reader, &prev->stars.parts[i], &next->stars.parts[i], time_offset);

    if (event != logger_reader_event_null) {
      error("TODO");
    }
  }

}

/**
 * @brief Get the two particle records around the requested time.
 *
 * @param reader The #logger_reader.
 * @param prev (in) A record before the requested time. (out) The last record
 * before the time.
 * @param next (out) The first record after the requested time.
 * @param time_offset The offset of the requested time.
 *
 * @return The event encountered (update also the offset of the previous particle).
 */
enum logger_reader_event logger_reader_get_next_particle(
    struct logger_reader *reader,
    struct logger_particle *prev,
    struct logger_particle *next,
    size_t time_offset) {

  void *map = reader->log.log.map;
  size_t prev_offset = prev->offset;
  size_t next_offset = 0;

  /* Get the mask index of the special flags */
  const int spec_flag_ind =
      header_get_field_index(&reader->log.header, "special flags");
  if (spec_flag_ind < -1) {
    error("The logfile does not contain the special flags field.");
  }

  while (1) {
    /* Read the offset to the next particle */
    size_t mask = 0;
    logger_loader_io_read_mask(&reader->log.header, (char *)map + prev_offset,
                               &mask, &next_offset);

    /* Check if something special happened */
    if (mask & reader->log.header.masks[spec_flag_ind].mask) {
      /* Read the special flag's data */
      struct logger_particle tmp;
      logger_particle_read(&tmp, reader, prev_offset, /* Time */ -1,
                           logger_reader_const);

      /* Get the data from the flag */
      int data = 0;
      enum logger_special_flags flag = logger_unpack_flags_and_data(
        tmp.flag, &data);

      /* Save the offset information */
      prev->offset = prev_offset;

      /* Deal with the different cases */
      switch(flag) {

        case logger_flag_change_type:
          if (data == swift_type_stars) {
            return logger_reader_event_stars;
          }
          else {
            error("Not implemented yet.");
          }

        case logger_flag_mpi_exit:
        case logger_flag_delete:
          return logger_reader_event_deleted;

        default:
          error("Failed");
      }
    }

    /* Are we at the end of the file? */
    if (next_offset == 0) {
      time_array_print(&reader->log.times);
      error(
          "End of file for particle %lli offset %zi when requesting time %g "
          "with offset %zi",
          prev->id, prev_offset,
          time_array_get_time(&reader->log.times, time_offset), time_offset);
    }

    next_offset += prev_offset;

    /* Have we found the next particle? */
    if (next_offset > time_offset) {
      break;
    }

    /* Update the previous offset */
    prev_offset = next_offset;
  }

  /* Read the previous offset if required */
  if (prev_offset != prev->offset) {
    logger_particle_read(prev, reader, prev_offset, /* Time */ 0,
                         logger_reader_const);
  }

  /* Read the next particle */
  logger_particle_read(next, reader, next_offset, /* Time */ 0,
                       logger_reader_const);

  return logger_reader_event_null;
}

/**
 * @brief Get the two particle records around the requested time.
 *
 * @param reader The #logger_reader.
 * @param prev (in) A record before the requested time. (out) The last record
 * before the time.
 * @param next (out) The first record after the requested time.
 * @param time_offset The offset of the requested time.
 *
 * @return The event encountered (update also the offset of the previous particle).
 */
enum logger_reader_event logger_reader_get_next_gparticle(
                                                         struct logger_reader *reader,
                                                         struct logger_gparticle *prev,
                                                         struct logger_gparticle *next,
                                                         size_t time_offset) {
  error("TODO");
}

/**
 * @brief Get the two particle records around the requested time.
 *
 * @param reader The #logger_reader.
 * @param prev (in) A record before the requested time. (out) The last record
 * before the time.
 * @param next (out) The first record after the requested time.
 * @param time_offset The offset of the requested time.
 *
 * @return The event encountered (update also the offset of the previous particle).
 */
enum logger_reader_event logger_reader_get_next_sparticle(
                                                         struct logger_reader *reader,
                                                         struct logger_sparticle *prev,
                                                         struct logger_sparticle *next,
                                                         size_t time_offset) {

  error("TODO");
}
