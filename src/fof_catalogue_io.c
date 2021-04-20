/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020  Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Local headers */
#include "engine.h"
#include "fof.h"
#include "hydro_io.h"

void write_fof_hdf5_header(hid_t h_file, const struct engine *e,
                           const struct fof_props *LOS_params) {

  /* Open header to write simulation properties */
  hid_t h_grp =
      H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_TIME);
  const double factor_length = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute_d(h_grp, "Time", dblTime);
  io_write_attribute_d(h_grp, "Dimension", (int)hydro_dimension);
  io_write_attribute_d(h_grp, "Redshift", e->cosmology->z);
  io_write_attribute_d(h_grp, "Scale-factor", e->cosmology->a);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);

  /* Write out the particle types */
  io_write_part_type_names(h_grp);

  /* Write out the time-base */
  if (e->policy | engine_policy_cosmology) {
    io_write_attribute_d(h_grp, "TimeBase_dloga", e->time_base);
    const double delta_t = cosmology_get_timebase(e->cosmology, e->ti_current);
    io_write_attribute_d(h_grp, "TimeBase_dt", delta_t);
  } else {
    io_write_attribute_d(h_grp, "TimeBase_dloga", 0);
    io_write_attribute_d(h_grp, "TimeBase_dt", e->time_base);
  }

  /* Store the time at which the snapshot was written */
  time_t tm = time(NULL);
  struct tm *timeinfo = localtime(&tm);
  char snapshot_date[64];
  strftime(snapshot_date, 64, "%T %F %Z", timeinfo);
  io_write_attribute_s(h_grp, "Snapshot date", snapshot_date);

  /* GADGET-2 legacy values */
  /* Number of particles of each type */
  long long N_total[swift_type_count] = {0};
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};
  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);
  }
  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, N_total,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total", UINT, numParticles,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total_HighWord", UINT,
                     numParticlesHighWord, swift_type_count);
  double MassTable[swift_type_count] = {0};
  io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
  io_write_attribute(h_grp, "InitialMassTable", DOUBLE,
                     e->s->initial_mean_mass_particles, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                     swift_type_count);
  io_write_attribute_i(h_grp, "NumFilesPerSnapshot", e->nr_nodes);
  io_write_attribute_i(h_grp, "ThisFile", e->nodeID);
  io_write_attribute_s(h_grp, "OutputType", "FOF");

  /* Close group */
  H5Gclose(h_grp);

  io_write_meta_data(h_file, e, e->internal_units, e->snapshot_units);
}

void write_fof_hdf5_catalogue(const struct fof_props *props,
                              const struct engine *e) {}
