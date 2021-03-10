/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_MULTI_SOFTENING_GRAVITY_PART_H
#define SWIFT_MULTI_SOFTENING_GRAVITY_PART_H

#include "fof_struct.h"
#include "logger.h"
#include "mpi_part.h"

/* Gravity particle. */
struct gpart {

  /*! Particle ID. If negative, it is the negative offset of the #part with
     which this gpart is linked. */
  long long id_or_neg_offset;

  /*! Particle position. */
  double x[3];

  /*! Particle velocity. */
  float v_full[3];

  /*! Particle acceleration from the tree. */
  float a_grav[3];

  /*! Particle acceleration from the mesh. */
  float a_grav_mesh[3];

#ifndef SWIFT_GRAVITY_NO_POTENTIAL

  /*! Gravitational potential from the tree.*/
  float potential;

  /*! Gravitational potential from the mesh.*/
  float potential_mesh;

#endif

  /*! Particle mass. */
  float mass;

  /*! Norm of the acceleration at the previous step. */
  float old_a_grav_norm;

  /*! Current co-moving spline softening of the particle */
  float epsilon;

  /*! Particle FoF properties (group ID, group size, ...) */
  struct fof_gpart_data fof_data;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Type of the #gpart (DM, gas, star, ...) */
  enum part_type type;

#ifdef WITH_LOGGER
  /* Additional data for the particle logger */
  struct logger_part_data logger_data;
#endif

#ifdef HAVE_VELOCIRAPTOR_ORPHANS
  /* Flag to indicate this particle should be output at subsequent VR
     invocations because it was the most bound in a group at some point */
  char has_been_most_bound;
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

  /* Time of the last mesh kick */
  integertime_t ti_kick_mesh;

  /* Has this particle been initialised? */
  int initialised;

  /* Total number of gparts this gpart interacted with */
  long long num_interacted;
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  /* Acceleration taken from each component of the tree */
  float a_grav_p2p[3];
  float a_grav_m2p[3];
  float a_grav_m2l[3];

  /* Brute-force particle accelerations */
  double a_grav_exact[3];
  double a_grav_exact_long[3];
  double a_grav_exact_short[3];

  /* Brute-force particle potential. */
  double potential_exact;

  /* Type specific interaction counters */
  long long num_interacted_m2l;
  long long num_interacted_m2p;
  long long num_interacted_p2p;
  long long num_interacted_pm;
#endif
};

/**
 * @brief Gravity particle data when on a foreign node.
 *
 * We only need a small fraction of the information on foreign
 * nodes as the particles don't get updated there.
 */
struct gpart_foreign {

  /*! Particle position. */
  double x[3];

  /*! Particle mass. */
  float mass;

  /*! Particle softening */
  float epsilon;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Type of the #gpart (DM, gas, star, ...) */
  enum part_type type;
};

#ifdef WITH_MPI

/**
 * @brief Register the MPI types to send a subset of a #gpart and receive it
 * as a #gpart_foreign.
 *
 * We here list the fields of a #gpart to be sent and received over MPI.
 *
 * @param gpart_send_mpi_type The MPI type for the sending.
 * @param gpart_send_mpi_type The MPI type for the receiving.
 */
INLINE static void gravity_create_MPI_types(MPI_Datatype *gpart_send_mpi_type,
                                            MPI_Datatype *gpart_recv_mpi_type) {
  create_indexed_mpi_type(struct gpart, *gpart_send_mpi_type, x, mass, epsilon,
                          time_bin, type);
  create_indexed_mpi_type(struct gpart_foreign, *gpart_recv_mpi_type, x, mass,
                          epsilon, time_bin, type);
}
#endif /* WITH_MPI */

#endif /* SWIFT_MULTI_SOFTENING_GRAVITY_PART_H */
