/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Config parameters. */
#include "../config.h"

#if defined(WITH_LOGGER) && defined(HAVE_HDF5) && !defined(WITH_MPI)

/* Some standard headers. */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "logger_io.h"

/* Local includes. */
#include "chemistry_io.h"
#include "common_io.h"
#include "cooling.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "gravity_io.h"
#include "gravity_properties.h"
#include "hydro_io.h"
#include "hydro_properties.h"
#include "io_properties.h"
#include "kernel_hydro.h"
#include "parallel_io.h"
#include "part.h"
#include "serial_io.h"
#include "single_io.h"
#include "stars_io.h"
#include "units.h"
#include "xmf.h"

/**
 * @brief Writes a logger index file
 *
 * @param log The #logger.
 * @param e The engine containing all the system.
 *
 * Creates an output file and writes the offset and id of particles
 * contained in the engine. If such a file already exists, it is erased and
 * replaced by the new one.
 *
 * Calls #error() if an error occurs.
 *
 */
void logger_write_index_file(struct logger *log, struct engine* e) {

  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Ntot = e->s->nr_gparts;

  struct part* parts = e->s->parts;
  struct xpart* xparts = e->s->xparts;
  // struct gpart* gparts = e->s->gparts;
  // struct gpart* dmparts = NULL;
  // struct spart* sparts = e->s->sparts;
  static int outputCount = 0;

  /* Number of unassociated gparts */
  const size_t Ndm = Ntot > 0 ? Ntot - (Ngas + Nstars) : 0;

  long long N_total[swift_type_count] = {Ngas, Ndm, 0, 0, Nstars, 0};

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%.100s_%04i.index", e->logger->base_name,
           outputCount);

  /* Open file */
  FILE *f = NULL;
  f = fopen(fileName, "wb");

  if (f == NULL) {
    error("Failed to open file %s", fileName);
  }

    
  /* Open header to write simulation properties */
  fprintf(f, "# Time=%g IntergerTime=%lli\n", e->time, e->ti_current);
  fprintf(f, "# Number of particles\n# ");
  for(int i = 0; i < swift_type_count; i++) {
    fprintf(f, "%lli ", N_total[i]);
  }
  fprintf(f, "\n");

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N_total[ptype] == 0) continue;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas:
        hydro_write_index(parts, xparts, Ngas, f);
        break;

      case swift_type_dark_matter:
        error("TODO");
        break;

      case swift_type_stars:
        error("TODO");
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

  }

  /* message("Done writing particles..."); */

  /* Close file */
  fclose(f);

  ++outputCount;
}

#endif /* WITH_LOGGER && HAVE_HDF5 && !WITH_MPI */
