/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
 *               2019 Joel Pfeffer (j.l.pfeffer@ljmu.ac.uk)
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
#ifndef SWIFT_MOSAICS_STAR_PART_H
#define SWIFT_MOSAICS_STAR_PART_H

/* Some standard headers. */
#include <stdlib.h>

/* Read additional aubgrid models */
#include "chemistry_struct.h"
#include "feedback_struct.h"
#include "tracers_struct.h"

/**
 * @brief Particle fields for the star particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct spart {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff_sort[3];

  /*! Particle velocity. */
  float v[3];

  /*! Star mass */
  float mass;

  /*! Particle smoothing length. */
  float h;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  /*! Union for the birth time and birth scale factor */
  union {

    /*! Birth time */
    float birth_time;

    /*! Birth scale factor */
    float birth_scale_factor;
  };

  /*! Scale-factor / time at which this particle last did enrichment */
  float last_enrichment_time;

  /*! Initial star mass */
  float mass_init;

  /*! Birth density */
  float birth_density;

  /*! Birth temperature */
  float birth_temperature;

  /*! SNII Feedback energy fraction */
  float SNII_f_E;

  /*! SNIa Feedback energy fraction */
  float SNIa_f_E;

  /*! last time an HII region was built (age of star in Myr) */
  float HIIregion_last_rebuild;

  /*! current time-step length of star particle */
  float star_timestep;

  /*! HII mass available for ionization (current) */
  float HIIregion_mass_to_ionize;

  /*! mass in kernel when HII region is built (for debugging) */
  float HIIregion_mass_in_kernel;

  /*! Feedback structure */
  struct feedback_spart_data feedback_data;

  /*! Tracer structure */
  struct tracers_xpart_data tracers_data;

  /*! Chemistry structure */
  struct chemistry_part_data chemistry_data;

  /*! Particle time bin */
  timebin_t time_bin;

  /*! Number of time-steps since the last enrichment step */
  char count_since_last_enrichment;

  /* -------------- Now the MOSAICS data ---------------- */

  /*! Second derivative of gravitational potential */
  /* upper symmetric 3*3 matrix:
   * tt[0] == xx
   * tt[1] == xy == yx
   * tt[2] == xz == zx
   * tt[3] == yy
   * tt[4] == yz == zy
   * tt[5] == zz
   */
  float tidal_tensor[3][6];

  // TODO just temporary for testing
  /*! Gravitational potential */
  float potential;

  /*! Birth pressure */
  float birth_pressure;

  /*! Birth subgrid density */
  float birth_subgrid_dens;

  /*! Birth subgrid sound speed */
  float sound_speed_subgrid;

  /*! Cluster formation efficiency */
  float CFE;

  /*! Exponential truncation to mass function */
  float Mcstar;

  /*! The local gas velocity dispersion at formation */
  float gas_vel_disp;

  /*! The local stellar velocity dispersion at formation */
  float star_vel_disp;

  /*! Gas fraction within that same region */
  float fgas;

  /*! Epicyclic frequency at formation */
  float kappa;

  /*! Circular frequency at formation */
  float Omega;

  /*! Local Toomre mass */
  float Toomre_mass;

  /*! Fraction of Mtoomre that may collapse to a GMC */
  float frac_collapse;

  /*! Flag denoting whether we need to run cluster formation  */
  int new_star;

  /*! Flag denoting whether particle has GCs */
  int gcflag;

  /*! In case we want tensors for particles without GCs */
  int calc_tensor;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_STARS

  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_STARS];

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS_STARS];
#endif

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Contains all the constants and parameters of the stars scheme
 */
struct stars_props {

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weighted number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /*! Are we overwriting the stars' birth time read from the ICs? */
  int overwrite_birth_time;

  /*! Value to set birth time of stars read from ICs */
  float spart_first_init_birth_time;

  /* ---------------- MOSAICS parameters ---------------- */

  /*! Flag to force tensor calculation for all stars */
  int calc_all_star_tensors;

  /*! King '66 density profile parameter */
  float W0;

  /*! Star formation efficiency for Mcstar */
  float SFE;

  /*! Use the subgrid turbulent velocity dispersion for CFE */
  int subgrid_gas_vel_disp;

  /*! Sound speed of cold ISM (m/s) */
  float Fixedcs;

  /* ----- Cluster formation efficiency parameters ------ */

  /*! star formation law */
  int sflaw;

  /*! GMC virial parameter */
  float qvir;

  /*! time of first SN (Myr) */
  float tsn;

  /*! time of determining the CFE (Myr) */
  float tview;

  /*! GMC surface density (Msun/pc^2) */
  float surfGMC;

  /*! maximum (protostellar core) SFE */
  float ecore;

  /*! turbulent/magnetic pressure ratio */
  float beta0;

  /*! SN/radiative feedback mode */
  int radfb;

  /* --------------- Conversion factors ----------------- */

  /*! Conversion factor from internal density unit to cgs */
  double density_to_kgm3;

  /*! Conversion factor from internal velocity unit to m/s */
  double velocity_to_ms;

  /*! Conversion factor from internal time unit to s */
  double time_to_cgs;

  /* TODO not sure if needed? */
  /*! Conversion factor from internal mass unit to solar mass */
  // double mass_to_solar_mass;
};

#endif /* SWIFT_MOSAICS_STAR_PART_H */
