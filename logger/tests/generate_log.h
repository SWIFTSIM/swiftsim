/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include config */
#include "../../config.h"

/* Local headers */
#include "hydro.h"
#include "logger.h"

/* Not all the fields are written at every step.
 * Here we define how often a few fields are written.
 */
#define period_rho 2
#define period_h 4

/**
 * @brief Generate the data of a bunch of particles.
 *
 * @param parts The list of particles.
 * @param xparts The list of extra particles.
 * @param nparts The number of particles.
 */
void generate_particles(struct part *parts, struct xpart *xparts, size_t nparts) {
  struct hydro_space hs;

  for (size_t i = 0; i < nparts; i++) {
    /* Set internal energy. */
    hydro_set_init_internal_energy(&parts[i], 100);

    /* Initialize particle. */
    hydro_first_init_part(&parts[i], &xparts[i]);
    hydro_init_part(&parts[i], &hs);

    for (int j = 0; j < 3; j++) {
      parts[i].x[j] = i;
      parts[i].v[j] = (j == 0) ? -1 : 0;
      parts[i].a_hydro[j] = (j == 1) ? 1e-2 : 0;
    }
    parts[i].h = 15;
    parts[i].rho = 50;
    parts[i].id = i;
    hydro_set_mass(&parts[i], 1.5);
    xparts[i].logger_data.last_offset = 0;

    /* Add time bin in order to skip particles. */
    parts[i].time_bin = (i % 10) + 1;
  }

}

/** Provides a integer time given the step number.*/
integertime_t get_integer_time(int step) { return step; }

/** Provides a double time given the step number. */
double get_double_time(int step) {
  const double time_base = 1e-4;
  return step * time_base;
}

/**
 * @brief Write a few particles during multiple time steps.
 *
 * As only the logger is tested, there is no need to really
 * evolve the particles.
 *
 * @param log The #logger_writer.
 * @param parts The list of particles.
 * @param xparts The list of x-particles.
 * @param nparts The number of particles.
 */
void write_particles(struct logger_writer *log, struct part *parts,
                     struct xpart *xparts, size_t nparts) {

  const int number_steps = 100;

  /* Loop over all the steps. */
  for (int i = 0; i < number_steps; i++) {
    integertime_t ti_int = get_integer_time(i);
    double ti_double = get_double_time(i);

    /* Mark the current time step in the particle logger file. */
    logger_log_timestamp(log, ti_int, ti_double, &log->timestamp_offset);
    /* Make sure that we have enough space in the particle logger file
     * to store the particles in current time step. */
    logger_ensure_size(log, nparts, /* number gpart */ 0, 0);

    /* Loop over all the particles. */
    for (size_t j = 0; j < nparts; j++) {

      /* Skip some particles. */
      if (i % parts[j].time_bin != 0) continue;

      /* Write a time information to check that the correct particle is read. */
      parts[j].x[0] = i;

      /* Write this particle. */
      unsigned int mask =
          logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
          logger_mask_data[logger_a].mask | logger_mask_data[logger_u].mask |
          logger_mask_data[logger_consts].mask;

      int number_particle_step = i / parts[j].time_bin;

      if (number_particle_step % period_h == 0)
        mask |= logger_mask_data[logger_h].mask;
      if (number_particle_step % period_rho == 0)
        mask |= logger_mask_data[logger_rho].mask;

      logger_log_part(log, &parts[j], mask, &xparts[j].logger_data.last_offset);
    }

    // TODO write index files.
  }

  /* Mark the current time step in the particle logger file. */
  integertime_t ti_int = get_integer_time(number_steps);
  double ti_double = get_double_time(number_steps);
  logger_log_timestamp(log, ti_int, ti_double, &log->timestamp_offset);
}

void generate_log(struct swift_params *params, struct part *parts,
                  struct xpart *xparts, size_t nparts) {
  /* Initialize the particles */
  generate_particles(parts, xparts, nparts);

  /* Initialize the writer */
  struct logger_writer log;
  logger_init(&log, params);

  /* Write file header */
  logger_write_file_header(&log);

  /* Write particles */
  write_particles(&log, parts, xparts, nparts);

  /* Cleanup the memory */
  logger_free(&log);
}
