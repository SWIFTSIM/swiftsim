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
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common_io.h"

/* define a few functions for the yaml file */
#define hid_t FILE *
#define io_write_attribute_s(file, params, value) {            \
    fprintf(file, "    %s: %s\n", params, value);              \
  }
#define io_write_attribute_d(file, params, value) {            \
    fprintf(file, "    %s: %g\n", params, value);              \
  }
#define io_write_attribute_f(file, params, value) {            \
    fprintf(file, "    %s: %g\n", params, value);              \
  }
#define io_write_attribute_i(file, params, value) {            \
    fprintf(file, "    %s: %i\n", params, value);              \
  }
#define io_write_attribute_l(file, params, value) {            \
    fprintf(file, "    %s: %li\n", params, value);             \
  }

#define io_write_attribute(file, params, type, value, dim) {	       \
    fprintf(file, "    %s: [", params);                                \
    for(int i = 0; i < dim-1; i++) {				       \
      fprintf(file, type ", ", value[i]);			       \
    }								       \
    fprintf(file, type "]\n", value[dim-1]);			       \
  }

/* This object's header. */
#include "logger_io.h"

/* Local includes. */
#include "chemistry_io.h"
#include "common_io.h"
#include "cooling_io.h"
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
#include "tracers_io.h"
#include "units.h"
#include "version.h"
#include "xmf.h"



/**
 * @brief Mapper function to copy #part or #gpart fields into a buffer.
 * WARNING Assumes two io_props in extra_data.
 */
void logger_io_copy_mapper(void* restrict temp, int N, void* restrict extra_data) {


  /* Get the io_props */
  const struct io_props *props = (const struct io_props*)(extra_data);
  const struct io_props props1 = props[0];
  const struct io_props props2 = props[1];

  /* Get the sizes */
  const size_t typeSize1 = io_sizeof_type(props1.type);
  const size_t copySize1 = typeSize1 * props1.dimension;
  const size_t typeSize2 = io_sizeof_type(props2.type);
  const size_t copySize2 = typeSize2 * props2.dimension;
  const size_t copySize = copySize1 + copySize2;

  /* How far are we with this chunk? */
  char* restrict temp_c = (char*)temp;
  const ptrdiff_t delta = (temp_c - props1.start_temp_c) / copySize;

  /* Copy the memory to the buffer */
  for (int k = 0; k < N; k++) {
    memcpy(&temp_c[k * copySize], props1.field + (delta + k) * props1.partSize,
           copySize1);
    memcpy(&temp_c[k * copySize + copySize1], props2.field + (delta + k) * props2.partSize,
           copySize2);
  }
}
/**
 * @brief Writes the data array in the index file.
 *
 * @param e The #engine we are writing from.
 * @param f The file to use.
 * @param props The #io_props array.
 * @param n_props The number of element in @props.
 * @param N The number of particles to write.
 */
void writeIndexArray(const struct engine* e, FILE *f,
		     struct io_props *props, size_t n_props, size_t N) {

  /* Check that the assumptions are corrects */
  if (n_props != 2)
    error("Not implemented: The index file can only write two props.");
  
  if (props[0].dimension != 1 || props[1].dimension != 1)
    error("Not implemented: cannot use multidimensional data");

  /* Get a few variables */
  const size_t typeSize = io_sizeof_type(props[0].type) +
    io_sizeof_type(props[1].type);

  const size_t num_elements = N;

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (posix_memalign((void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy the particle data to the temporary buffer */
  /* Set initial buffer position */
  props[0].start_temp_c = temp;
  props[1].start_temp_c = temp;

  /* Copy the whole thing into a buffer */
  threadpool_map((struct threadpool*)&e->threadpool, logger_io_copy_mapper, temp,
		 N, typeSize, 0, props);

  /* Write data to file */
  fwrite(temp, typeSize, num_elements, f);

  /* Free everything */
  free(temp);
}

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
 * An index file is constructed by writing first a few variables (e.g. time, number of particles,
 * if the file is sorted, ...) and then an array of index and offset for each particle type.
 *
 * Calls #error() if an error occurs.
 *
 */
void logger_write_index_file(struct logger *log, struct engine* e) {

  struct part* parts = e->s->parts;
  struct xpart* xparts = e->s->xparts;
  struct gpart* gparts = e->s->gparts;
  // struct spart* sparts = e->s->sparts;
  static int outputCount = 0;

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  /* const size_t Nstars = e->s->nr_sparts; */
  /* const size_t Nblackholes = e->s->nr_bparts; */

  /* Number of particles that we will write */
  const size_t Ntot_written =
      e->s->nr_gparts - e->s->nr_inhibited_gparts - e->s->nr_extra_gparts;
  const size_t Ngas_written =
      e->s->nr_parts - e->s->nr_inhibited_parts - e->s->nr_extra_parts;
  const size_t Nstars_written =
      e->s->nr_sparts - e->s->nr_inhibited_sparts - e->s->nr_extra_sparts;
  const size_t Nblackholes_written =
      e->s->nr_bparts - e->s->nr_inhibited_bparts - e->s->nr_extra_bparts;
  const size_t Nbaryons_written =
      Ngas_written + Nstars_written + Nblackholes_written;
  const size_t Ndm_written =
      Ntot_written > 0 ? Ntot_written - Nbaryons_written : 0;

  /* Format things in a Gadget-friendly array */
  long long N_total[swift_type_count] = {
      (long long)Ngas_written,   (long long)Ndm_written,        0, 0,
      (long long)Nstars_written, (long long)Nblackholes_written};

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

    
  /* Write double time */
  fwrite(&e->time, sizeof(double), 1, f);

  /* Write integer time */
  fwrite(&e->ti_current, sizeof(integertime_t), 1, f);

  /* Write number of particles */
  fwrite(N_total, sizeof(long long), swift_type_count, f);

  /* Write if the file is sorted */
  const char sorted = 0;
  fwrite(&sorted, sizeof(char), 1, f);

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N_total[ptype] == 0) continue;

    /* Number of properties (the code cannot deal with more than two props
       per particle type) */
    size_t N = 0;
    int num_fields = 0;
    struct io_props list[2];

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;
    struct spart* sparts_written = NULL;
    struct bpart* bparts_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas:
        if (Ngas == Ngas_written) {

          /* No inhibted particles: easy case */
          N = Ngas;
          num_fields += hydro_write_index(parts, xparts, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ngas_written;

          /* Allocate temporary arrays */
          if (swift_memalign("parts_written", (void**)&parts_written,
                             part_align,
                             Ngas_written * sizeof(struct part)) != 0)
            error("Error while allocating temporary memory for parts");
          if (swift_memalign("xparts_written", (void**)&xparts_written,
                             xpart_align,
                             Ngas_written * sizeof(struct xpart)) != 0)
            error("Error while allocating temporary memory for xparts");

          /* Collect the particles we want to write */
          io_collect_parts_to_write(parts, xparts, parts_written,
                                    xparts_written, Ngas, Ngas_written);

          /* Select the fields to write */
          num_fields += hydro_write_index(parts, xparts, list);
        } break;

      case swift_type_dark_matter:
        if (Ntot == Ndm_written) {

          /* This is a DM-only run without inhibited particles */
          N = Ntot;
          num_fields += darkmatter_write_index(gparts, list);
        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ndm_written;

          /* Allocate temporary array */
          if (swift_memalign("gparts_written", (void**)&gparts_written,
                             gpart_align,
                             Ndm_written * sizeof(struct gpart)) != 0)
            error("Error while allocating temporary memory for gparts");

          /* Collect the non-inhibited DM particles from gpart */
	  const int with_stf = 0;
          io_collect_gparts_to_write(gparts, e->s->gpart_group_data,
                                     gparts_written, gpart_group_data_written,
                                     Ntot, Ndm_written, with_stf);

          /* Select the fields to write */
          num_fields += darkmatter_write_index(gparts, list);
	} break;

      case swift_type_stars:
        error("TODO");
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    if (num_fields != 2) {
      error("The code expects only two fields per particle type for the logger");
    }

    /* Write ids */
    writeIndexArray(e, f, list, num_fields, N);

    /* Free temporary arrays */
    if (parts_written) swift_free("parts_written", parts_written);
    if (xparts_written) swift_free("xparts_written", xparts_written);
    if (gparts_written) swift_free("gparts_written", gparts_written);
    if (gpart_group_data_written)
      swift_free("gpart_group_written", gpart_group_data_written);
    if (sparts_written) swift_free("sparts_written", sparts_written);
    if (bparts_written) swift_free("bparts_written", bparts_written);
  }

  /* Close file */
  fclose(f);

  ++outputCount;
}

/**
 * @brief Writes the current Unit System
 * @param file The (opened) params file in which to write
 * @param us The unit_system to dump
 * @param groupName The name of the section to write to
 */
void logger_io_write_unit_system(FILE *file, const struct unit_system* us,
				 const char* groupName) {

  fprintf(file, "%s:\n", groupName);
  
  io_write_attribute_d(file, "Unit mass in cgs (U_M)",
		    units_get_base_unit(us, UNIT_MASS));
  io_write_attribute_d(file, "Unit length in cgs (U_L)",
		    units_get_base_unit(us, UNIT_LENGTH));
  io_write_attribute_d(file, "Unit time in cgs (U_t)",
		    units_get_base_unit(us, UNIT_TIME));
  io_write_attribute_d(file, "Unit current in cgs (U_I)",
		    units_get_base_unit(us, UNIT_CURRENT));
  io_write_attribute_d(file, "Unit temperature in cgs (U_T)",
		    units_get_base_unit(us, UNIT_TEMPERATURE));

  fprintf(file, "\n");
}


/**
 * @brief Writes the code version to the file
 * @param file The (opened) params file in which to write
 */
void logger_io_write_code_description(FILE *file) {

  fprintf(file, "Code:\n");

  io_write_attribute_s(file, "Code", "SWIFT");
  io_write_attribute_s(file, "Code Version", package_version());
  io_write_attribute_s(file, "Compiler Name", compiler_name());
  io_write_attribute_s(file, "Compiler Version", compiler_version());
  io_write_attribute_s(file, "Git Branch", git_branch());
  io_write_attribute_s(file, "Git Revision", git_revision());
  io_write_attribute_s(file, "Git Date", git_date());
  io_write_attribute_s(file, "Configuration options",
                    configuration_options());
  io_write_attribute_s(file, "CFLAGS", compilation_cflags());
  io_write_attribute_s(file, "HDF5 library version", hdf5_version());
  io_write_attribute_s(file, "Thread barriers", thread_barrier_version());
  io_write_attribute_s(file, "Allocators", allocator_version());
#ifdef HAVE_FFTW
  io_write_attribute_s(file, "FFTW library version", fftw3_version());
#endif
#ifdef HAVE_LIBGSL
  io_write_attribute_s(file, "GSL library version", libgsl_version());
#endif
#ifdef WITH_MPI
  io_write_attribute_s(file, "MPI library", mpi_version());
#ifdef HAVE_METIS
  io_write_attribute_s(file, "METIS library version", metis_version());
#endif
#ifdef HAVE_PARMETIS
  io_write_attribute_s(file, "ParMETIS library version",
                       parmetis_version());
#endif
#else
  io_write_attribute_s(file, "MPI library", "Non-MPI version of SWIFT");
#endif

  fprintf(file, "\n");
}


/**
 * @brief Write the #engine policy to the file.
 * @param file File to write to.
 * @param e The #engine to read the policy from.
 */
void logger_io_write_engine_policy(FILE *file, const struct engine* e) {

  fprintf(file, "Policy:\n");

  for (int i = 1; i < engine_maxpolicy; ++i) {
    if (e->policy & (1 << i)) {
      io_write_attribute_i(file, engine_policy_names[i + 1], 1);
    }
    else {
      io_write_attribute_i(file, engine_policy_names[i + 1], 0);
    }
  }

  fprintf(file, "\n");
}

/**
 * @brief Write the #hydro_props to the file.
 * @param file File to write to.
 * @param p The #hydro_props to write.
 */
void logger_hydro_props_print_snapshot(FILE *file, const struct hydro_props *p) {

  fprintf(file, "HydroScheme:\n");

  eos_print_snapshot(file, &eos);

  io_write_attribute_i(file, "Dimension", (int)hydro_dimension);
  io_write_attribute_s(file, "Scheme", SPH_IMPLEMENTATION);
  io_write_attribute_s(file, "Kernel function", kernel_name);
  io_write_attribute_f(file, "Kernel target N_ngb", p->target_neighbours);
  io_write_attribute_f(file, "Kernel delta N_ngb", p->delta_neighbours);
  io_write_attribute_f(file, "Kernel eta", p->eta_neighbours);
  io_write_attribute_f(file, "Smoothing length tolerance", p->h_tolerance);
  io_write_attribute_f(file, "Maximal smoothing length [internal units]",
                       p->h_max);
  io_write_attribute_f(file, "CFL parameter", p->CFL_condition);
  io_write_attribute_f(file, "Volume log(max(delta h))",
                       p->log_max_h_change);
  io_write_attribute_f(file, "Volume max change time-step",
                       pow_dimension(expf(p->log_max_h_change)));
  io_write_attribute_i(file, "Max ghost iterations",
                       p->max_smoothing_iterations);
  io_write_attribute_f(file, "Minimal temperature", p->minimal_temperature);
  io_write_attribute_f(file,
                       "Minimal energy per unit mass [internal units]",
                       p->minimal_internal_energy);
  io_write_attribute_f(file, "Initial temperature", p->initial_temperature);
  io_write_attribute_f(file,
                       "Initial energy per unit mass [internal units]",
                       p->initial_internal_energy);
  io_write_attribute_f(file, "Hydrogen mass fraction",
                       p->hydrogen_mass_fraction);
  io_write_attribute_f(file, "Hydrogen ionization transition temperature",
                       p->hydrogen_ionization_temperature);
  io_write_attribute_f(file, "Max v_sig ratio (limiter)",
                       const_limiter_max_v_sig_ratio);

  /* Write out the implementation-dependent viscosity parameters
   * (see hydro/SCHEME/hydro_parameters.h for this implementation) */
  viscosity_print_snapshot(file, &(p->viscosity));

  /* Same for the diffusion */
  diffusion_print_snapshot(file, &(p->diffusion));

}

/**
 * @brief Write the #gravity_props to the file.
 * @param file File to write to.
 * @param p The #gravity_props to write.
 */
void logger_gravity_props_print_snapshot(FILE *file,
                                  const struct gravity_props *p) {

  io_write_attribute_f(file, "Time integration eta", p->eta);
  io_write_attribute_s(file, "Softening style",
                    kernel_gravity_softening_name);
  io_write_attribute_f(
   file, "Comoving softening length [internal units]",
   p->epsilon_comoving * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(
   file,
   "Comoving Softening length (Plummer equivalent)  [internal units]",
   p->epsilon_comoving);
  io_write_attribute_f(
   file, "Maximal physical softening length  [internal units]",
   p->epsilon_max_physical * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(file,
                    "Maximal physical softening length (Plummer equivalent) "
                    " [internal units]",
                    p->epsilon_max_physical);
  io_write_attribute_f(file, "Opening angle", p->theta_crit);
  io_write_attribute_s(file, "Scheme", GRAVITY_IMPLEMENTATION);
  io_write_attribute_i(file, "MM order", SELF_GRAVITY_MULTIPOLE_ORDER);
  io_write_attribute_f(file, "Mesh a_smooth", p->a_smooth);
  io_write_attribute_f(file, "Mesh r_cut_max ratio", p->r_cut_max_ratio);
  io_write_attribute_f(file, "Mesh r_cut_min ratio", p->r_cut_min_ratio);
  io_write_attribute_f(file, "Tree update frequency",
                    p->rebuild_frequency);
  io_write_attribute_s(file, "Mesh truncation function",
                    kernel_long_gravity_truncation_name);
}


/**
 * @brief Write the #cosmology to the file.
 * @param file File to write to.
 * @param c The #cosmology to write.
 */
void logger_cosmology_write_model(FILE *file, const struct cosmology *c) {

  io_write_attribute_d(file, "a_beg", c->a_begin);
  io_write_attribute_d(file, "a_end", c->a_end);
  io_write_attribute_d(file, "time_beg [internal units]", c->time_begin);
  io_write_attribute_d(file, "time_end [internal units]", c->time_end);
  io_write_attribute_d(file, "Universe age [internal units]", c->time);
  io_write_attribute_d(file, "h", c->h);
  io_write_attribute_d(file, "H0 [internal units]", c->H0);
  io_write_attribute_d(file, "H [internal units]", c->H);
  io_write_attribute_d(file, "Hubble time [internal units]", c->Hubble_time);
  io_write_attribute_d(file, "Omega_m", c->Omega_m);
  io_write_attribute_d(file, "Omega_r", c->Omega_r);
  io_write_attribute_d(file, "Omega_b", c->Omega_b);
  io_write_attribute_d(file, "Omega_k", c->Omega_k);
  io_write_attribute_d(file, "Omega_lambda", c->Omega_lambda);
  io_write_attribute_d(file, "w_0", c->w_0);
  io_write_attribute_d(file, "w_a", c->w_a);
  io_write_attribute_d(file, "w", c->w);
  io_write_attribute_d(file, "Critical density [internal units]",
		    c->critical_density);
}


/**
 * @brief Write the contents of the parameter structure to a file
 *
 * @param params Structure that holds the parameters
 * @param file The file.
 * @param write_used Write used fields or unused fields.
 */
void logger_parser_write_params_to_hdf5(const struct swift_params *params,
					FILE *file, int write_used) {

  for (int i = 0; i < params->paramCount; i++) {
    if (write_used && !params->data[i].used)
      continue;
    else if (!write_used && params->data[i].used)
      continue;
    io_write_attribute_s(file, params->data[i].name, params->data[i].value);
  }
}

/**
 * @brief Write the parameters into a yaml file.
 *
 * @params log The #logger.
 * @params e The #engine.
 */
void logger_write_description(struct logger *log, struct engine* e) {
  const struct unit_system *internal_units = e->internal_units;
  const struct unit_system *snapshot_units = e->snapshot_units;

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%.100s.yml", e->logger->base_name);

  /* Open file */
  FILE *f = NULL;
  f = fopen(fileName, "wb");

  if (f == NULL) {
    error("Failed to open file %s", fileName);
  }

  /* Convert basic output information to snapshot units */
  const double factor_length =
    units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  /* Print the relevant information and print status */
  fprintf(f, "Header:\n");
  io_write_attribute(f, "BoxSize", "%g", dim, 3);
  const int dimension = (int)hydro_dimension;
  io_write_attribute_i(f, "Dimension", dimension);
  io_write_attribute_s(f, "Code", "SWIFT");
  time_t tm = time(NULL);
  io_write_attribute_s(f, "Snapshot date", ctime(&tm));

  double MassTable[swift_type_count] = {0};
  io_write_attribute(f, "MassTable", "%g", MassTable, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(f, "Flag_Entropy_ICs", "%u", flagEntropy,
		  swift_type_count);

  fprintf(f, "\n");

  /* Print the code version */
  logger_io_write_code_description(f);

  /* Print the run's policy */
  logger_io_write_engine_policy(f, e);

  /* Print the SPH parameters */
  if (e->policy & engine_policy_hydro) {
    logger_hydro_props_print_snapshot(f, e->hydro_properties);
    hydro_write_flavour(f);

    fprintf(f, "\n");
  }

  /* Print the subgrid parameters */
  fprintf(f, "SubgridScheme:\n");

  entropy_floor_write_flavour(f);
  cooling_write_flavour(f, e->cooling_func);
  chemistry_write_flavour(f);
  tracers_write_flavour(f);
  fprintf(f, "\n");

  /* Print the gravity parameters */
  if (e->policy & engine_policy_self_gravity) {
    fprintf(f, "GravityScheme:\n");
    logger_gravity_props_print_snapshot(f, e->gravity_properties);

    fprintf(f, "\n");
  }

  /* Print the stellar parameters */
  if (e->policy & engine_policy_stars) {
    fprintf(f, "StarsScheme:\n");
    stars_props_print_snapshot(f, e->stars_properties);
    fprintf(f, "\n");
  }

  /* Print the cosmological model  */
  fprintf(f, "Cosmology:\n");
  if (e->policy & engine_policy_cosmology) {
    io_write_attribute_i(f, "Cosmological run", 1);
  }
  else {
    io_write_attribute_i(f, "Cosmological run", 0);
  }
  logger_cosmology_write_model(f, e->cosmology);
  fprintf(f, "\n");

  /* Print the runtime parameters */
  fprintf(f, "Parameters:\n");
  logger_parser_write_params_to_hdf5(e->parameter_file, f, 1);
  fprintf(f, "\n");

  /* Print the runtime unused parameters */
  fprintf(f, "UnusedParameters:\n");
  logger_parser_write_params_to_hdf5(e->parameter_file, f, 0);
  fprintf(f, "\n");

  /* Print the system of Units used in the spashot */
  // TODO Improve output
  logger_io_write_unit_system(f, snapshot_units, "Units");

  /* Print the system of Units used internally */
  logger_io_write_unit_system(f, internal_units, "InternalCodeUnits");

  /* TODO implement cells */

  /* Close file */
  fclose(f);
  
}

#endif /* HAVE_HDF5 */
