/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Some standard headers. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <hdf5.h>

/* This object's header. */
#include "lightcone.h"

/* Local headers */
#include "common_io.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "lightcone_io.h"
#include "lightcone_replications.h"
#include "lock.h"
#include "parser.h"
#include "particle_buffer.h"
#include "periodic.h"
#include "restart.h"
#include "space.h"
#include "timeline.h"

/* Whether to dump the replication list */
//#define DUMP_REPLICATIONS
#ifdef DUMP_REPLICATIONS
static int output_nr = 0;
#endif

/* MPI rank for diagnostic messages */
extern int engine_rank;


static void lightcone_allocate_buffers(struct lightcone_props *props) {

  /* Initialize particle output buffers */
  const size_t elements_per_block = (size_t) props->buffer_chunk_size;

  if(props->use_type[swift_type_gas]) {
    particle_buffer_init(&props->buffer[swift_type_gas],
                         sizeof(struct lightcone_gas_data),
                         elements_per_block, "lightcone_gas");
  }
  
  if(props->use_type[swift_type_dark_matter]) {
    particle_buffer_init(&props->buffer[swift_type_dark_matter],
                         sizeof(struct lightcone_dark_matter_data),
                         elements_per_block, "lightcone_dm");
  }

  if(props->use_type[swift_type_dark_matter_background]) {  
    particle_buffer_init(&props->buffer[swift_type_dark_matter_background],
                         sizeof(struct lightcone_dark_matter_data),
                         elements_per_block, "lightcone_dm_bg");
  }

  if(props->use_type[swift_type_stars]) {  
    particle_buffer_init(&props->buffer[swift_type_stars],
                         sizeof(struct lightcone_stars_data),
                         elements_per_block, "lightcone_stars");
  }

  if(props->use_type[swift_type_black_hole]) {
  particle_buffer_init(&props->buffer[swift_type_black_hole],
                       sizeof(struct lightcone_black_hole_data),
                       elements_per_block, "lightcone_bh");
  }

  if(props->use_type[swift_type_neutrino]) {
    particle_buffer_init(&props->buffer[swift_type_neutrino],
                         sizeof(struct lightcone_neutrino_data),
                         elements_per_block, "lightcone_neutrino");
  }  
}


/**
 * @brief Dump lightcone_props struct to the output stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to write to.
 */
void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  /* Don't dump the replication list - will regenerate it as needed */
  struct lightcone_props tmp = *props;
  tmp.replication_list.nrep = 0;
  tmp.replication_list.replication = NULL;
  tmp.have_replication_list = 0;

  /* Don't write out particle buffers - must flush before dumping restart. */
  memset(tmp.buffer, 0, sizeof(struct particle_buffer)*swift_type_count);

  restart_write_blocks((void *) &tmp, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");
}


/**
 * @brief Restore lightcone_props struct from the input stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to read from.
 */
void lightcone_struct_restore(struct lightcone_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                      NULL, "lightcone_props");
  lightcone_allocate_buffers(props);
}


/**
 * @brief Initialise the properties of the lightcone code.
 *
 * @param props the #lightcone_props structure to fill.
 * @param params the parameter file parser.
 */
void lightcone_init(struct lightcone_props *props,
                    const struct space *s,
                    const struct cosmology *cosmo,
                    struct swift_params *params) {
  
  /* Whether we generate lightcone output */
  props->enabled = 1;

  /* Which particle types to use */
  for(int i=0; i<swift_type_count; i+=1)
    props->use_type[i] = 0;
  props->use_type[swift_type_gas] = parser_get_param_int(params, "Lightcone:use_gas");
  props->use_type[swift_type_dark_matter] = parser_get_param_int(params, "Lightcone:use_dm");
  props->use_type[swift_type_dark_matter_background] = parser_get_param_int(params, "Lightcone:use_dm_background");
  props->use_type[swift_type_stars] = parser_get_param_int(params, "Lightcone:use_stars");
  props->use_type[swift_type_black_hole] = parser_get_param_int(params, "Lightcone:use_black_hole");
  props->use_type[swift_type_neutrino] = parser_get_param_int(params, "Lightcone:use_neutrino");

  /* Base name for output files */
  parser_get_param_string(params, "Lightcone:basename", props->basename);

  /* Redshift range for the lightcone */
  props->z_min = parser_get_param_double(params, "Lightcone:z_min");
  props->z_max = parser_get_param_double(params, "Lightcone:z_max");

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, "Lightcone:observer_position", 3,
                                props->observer_position);

  /* Write particles to disk if this many or more are in the buffer */
  props->max_particles_buffered = parser_get_opt_param_int(params, "Lightcone:max_particles_buffered", 100000);

  /* Chunk size for particles buffered in memory  */
  props->buffer_chunk_size = parser_get_opt_param_int(params, "Lightcone:buffer_chunk_size", 20000);

  /* Chunk size for particles buffered in memory  */
  props->hdf5_chunk_size = parser_get_opt_param_int(params, "Lightcone:hdf5_chunk_size", 16384);

  /* Get the size of the simulation box */
  props->boxsize = s->dim[0];
  if(s->dim[1] != s->dim[0] || s->dim[2] != s->dim[0])
    error("Lightcones require a cubic simulation box.");

  /* Get top level cell size */
  props->cell_width = s->width[0];
  if(s->width[1] != s->width[0] || s->width[2] != s->width[0])
    error("Lightcones require cubic top level cells.");

  /* Store expansion factors at z_min, z_max */
  props->a_at_z_min = 1.0/(1.0+props->z_min);
  props->a_at_z_max = 1.0/(1.0+props->z_max);

  /* Find corresponding distance squared */
  props->r2_max = pow(cosmology_get_comoving_distance(cosmo, props->a_at_z_max), 2.0);
  props->r2_min = pow(cosmology_get_comoving_distance(cosmo, props->a_at_z_min), 2.0);

  /* Estimate number of particles which will be output.

     Assumptions:
     - flat cosmology (haven't implemented comoving volume calculation for non-flat)
     - uniform box
  */
  const long long nr_gparts = s->nr_gparts;
  long long total_nr_gparts;
#ifdef WITH_MPI
  MPI_Reduce(&nr_gparts, &total_nr_gparts, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  total_nr_gparts = nr_gparts;
#endif
  if(engine_rank==0) {
    const double lightcone_rmax = cosmology_get_comoving_distance(cosmo, props->a_at_z_max);
    const double lightcone_rmin = cosmology_get_comoving_distance(cosmo, props->a_at_z_min);
    const double volume = 4./3.*M_PI*(pow(lightcone_rmax, 3.)-pow(lightcone_rmin, 3.));
    const long long est_nr_output = total_nr_gparts / pow(props->boxsize, 3.0) * volume;
    message("gparts in lightcone (if uniform box+flat cosmology): %lld", est_nr_output);
  }

  /* Initially have no replication list */
  props->have_replication_list = 0;
  props->ti_old = 0;
  props->ti_current = 0;

  /* Initialize various counters */
  for(int i=0; i<swift_type_count; i+=1) {
    props->tot_num_particles_written[i] = 0;
    props->num_particles_written_to_file[i] = 0;
  }
  props->current_file = -1;
  
  /* Always start a new file initially */
  props->start_new_file = 1;

  /* Allocate lightcone output buffers */
  lightcone_allocate_buffers(props);

}


/**
 * @brief Flush any buffers which exceed the specified size.
 */
void lightcone_flush_buffers(struct lightcone_props *props,
                             int flush_all, int end_file) {

  /* Will flush any buffers with more particles than this */
  size_t max_to_buffer = (size_t) props->max_particles_buffered;
  if(flush_all)max_to_buffer = 0;

  /* Count how many types have data to write out */
  int types_to_flush = 0;
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    if(props->use_type[ptype]) {
      const size_t num_to_write = particle_buffer_num_elements(&props->buffer[ptype]);
      if(num_to_write >= max_to_buffer && num_to_write > 0)types_to_flush += 1;
    }
  }
  
  /* Check if there's anything to do */
  if(types_to_flush>0) {
    
    /* We have data to flush, so open or create the output file */
    hid_t file_id;
    char fname[FILENAME_BUFFER_SIZE];
    if(props->start_new_file) {

      /* Get the name of the next file */
      props->current_file += 1;
      if(snprintf(fname, FILENAME_BUFFER_SIZE, "%s_%04d.%d.hdf5",
                  props->basename, props->current_file, engine_rank) < 0)
        error("Lightcone output filename truncated");

      /* Create the file */
      file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if(file_id < 0)error("Unable to create new lightcone file: %s", fname);

      /* We have now written no particles to the current file */
      for(int ptype=0; ptype<swift_type_count; ptype+=1)
        props->num_particles_written_to_file[ptype] = 0;    

      /* We no longer need to create a new file */
      props->start_new_file = 0;

    } else {

      /* Re-open an existing file */
      if(snprintf(fname, FILENAME_BUFFER_SIZE, "%s_%04d.%d.hdf5",
                  props->basename, props->current_file, engine_rank) < 0)
        error("Lightcone output filename truncated");
      file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
      if(file_id < 0)error("Unable to open current lightcone file: %s", fname);

    }

    /* Loop over particle types */
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      if(props->use_type[ptype]) {
        const size_t num_to_write = particle_buffer_num_elements(&props->buffer[ptype]);
        if(num_to_write >= max_to_buffer && num_to_write > 0) {
          switch(ptype) {
          case swift_type_gas:
            lightcone_write_gas(props, file_id, ptype);
            break;
          case swift_type_dark_matter:
          case swift_type_dark_matter_background:
            lightcone_write_dark_matter(props, file_id, ptype);
            break;
          case swift_type_stars:
            lightcone_write_stars(props, file_id, ptype);
            break;
          case swift_type_black_hole:
            lightcone_write_black_hole(props, file_id, ptype);
            break;
          case swift_type_neutrino:
            lightcone_write_neutrino(props, file_id, ptype);
            break;
          default:
            /* Don't support this type in lightcones */
            error("Trying to write out unsupported lightcone particle type");
          }
          particle_buffer_empty(&props->buffer[ptype]);
          props->num_particles_written_to_file[ptype] += num_to_write;          
        }
      }
    } /* Next particle type */

    /* We're done updating the output file */
    H5Fclose(file_id);
  }

  /* If we need to start a new file next time, record this */
  if(end_file)props->start_new_file = 1;

}


/**
 * @brief Flush all buffers and deallocate lightcone data.
 */
void lightcone_clean(struct lightcone_props *props) {

  /* Flush particle buffers */
  lightcone_flush_buffers(props, /* flush_all = */ 1, /* end_file = */ 1);

  /* Deallocate particle buffers */
  for(int i=0; i<swift_type_count; i+=1) {
    if(props->use_type[i])
      particle_buffer_free(&props->buffer[i]);
  }

  /* Free replication list, if we have one */
  if(props->have_replication_list)replication_list_clean(&props->replication_list);

}


/**
 * @brief Determine periodic copies of the simulation box which could
 * contribute to the lightcone.
 *
 *                     \
 *           \          \
 *            |         |
 * Obs      A |    B    | C
 *            |         |
 *           /          /
 *          R1         /
 *                    R0
 *
 * Consider a single particle being drifted. Here R0 is the comoving
 * distance to the time the particle is drifted FROM. R1 is the comoving
 * distance to the time the particle is drifted TO on this step.
 *
 * Particles which are beyond the lightcone surface at the start of
 * their drift (C) cannot cross the lightcone on this step if v < c.
 * Particles between the lightcone surfaces at the start and end of
 * their drift (B) may cross the lightcone (and certainly will if they
 * have zero velocity).
 *
 * Particles just within the lightcone surface at the start of their
 * drift (A) may be able to cross the lightcone due to their velocity so
 * we need to allow a boundary layer on the inside edge of the shell.
 * If we assume v < c, then we can use a layer of thickness R0-R1.
 *
 * Here we compute the earliest and latest times particles may be drifted
 * between, find the corresponding comoving distances R0 and R1, reduce
 * the inner distance by R0-R1, and find all periodic copies of the 
 * simulation box which overlap this spherical shell.
 *
 * Later we use this list to know which periodic copies to check when
 * particles are drifted.
 *
 * @param props The #lightcone_props structure
 * @param cosmo The #cosmology structure
 * @param s The #space structure
 * @param time_min Start of the time step
 * @param time_max End of the time step
 */
void lightcone_init_replication_list(struct lightcone_props *props,
                                     const struct cosmology *cosmo,
                                     const integertime_t ti_old,
                                     const integertime_t ti_current,
                                     const double dt_max) {

  /* Deallocate the old list, if there is one */
  if(props->have_replication_list)replication_list_clean(&props->replication_list);

  /* Get the size of the simulation box */
  const double boxsize = props->boxsize;

  /* Get a lower limit on earliest time particle may be drifted from */
  float dt = cosmo->time_end - cosmo->time_begin;
  while (dt > dt_max) dt /= 2.f;
  timebin_t bin = get_time_bin(dt*cosmo->time_base_inv);
  integertime_t ti_lim = get_integer_time_begin(ti_old, bin);

  /* Get expansion factor at earliest and latest times particles might be drifted between */
  double a_current = cosmo->a_begin * exp(ti_current * cosmo->time_base);
  double a_old = cosmo->a_begin * exp(ti_lim * cosmo->time_base);
  if(a_old < cosmo->a_begin)a_old = cosmo->a_begin;

  /* Convert redshift range to a distance range */
  double lightcone_rmin = cosmology_get_comoving_distance(cosmo, a_current);
  double lightcone_rmax = cosmology_get_comoving_distance(cosmo, a_old);
  if(lightcone_rmin > lightcone_rmax)
    error("Lightcone has rmin > rmax");

  /* Allow inner boundary layer, assuming all particles have v < c.
     This is to account for particles moving during the time step. */
  double boundary = lightcone_rmax-lightcone_rmin;
  lightcone_rmin -= boundary;
  if(lightcone_rmin < 0)lightcone_rmin = 0;

  if(a_current < props->a_at_z_max || a_old > props->a_at_z_min) {
    /* Timestep does not overlap the lightcone redshift range */
    replication_list_init_empty(&props->replication_list);
  } else {
    /* Timestep may contribute particles to the lightcone */
    replication_list_init(&props->replication_list, boxsize,
                          props->cell_width,
                          props->observer_position,
                          lightcone_rmin, lightcone_rmax);
  }

  /* Record that we made the list */
  props->have_replication_list = 1;

  /* Store times we used to make the list, for consistency check later */
  props->ti_old = ti_lim;
  props->ti_current = ti_current;

  /* Report the size of the list */
#ifdef DUMP_REPLICATIONS
  if(engine_rank==0) {
    message("no. of replications to check: %d", props->replication_list.nrep);
    message("shell to search inner radius=%e, outer radius=%e", lightcone_rmin,
            lightcone_rmax);
  }
#endif

  /* Write out the list, if required */
#ifdef DUMP_REPLICATIONS
  if(engine_rank==0) {
    char fname[500];
    sprintf(fname, "replication_list.%d.txt", output_nr);
    FILE *fd_rep = fopen(fname, "w");
    fprintf(fd_rep, "# Observer x, y, z\n");
    fprintf(fd_rep, "%e, %e, %e\n", props->observer_position[0],
            props->observer_position[1], props->observer_position[2]); 
    fprintf(fd_rep, "# Box size, inner radius, outer radius\n");
    fprintf(fd_rep, "%e, %e, %e\n", boxsize, lightcone_rmin-boundary, lightcone_rmax);
    fprintf(fd_rep, "# x, y, z, rmin2, rmax2\n");
    replication_list_write(&props->replication_list, fd_rep);
    fclose(fd_rep);
    output_nr += 1;
  }
#endif
}
