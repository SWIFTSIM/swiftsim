/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MOSAICS_STARS_IO_H
#define SWIFT_MOSAICS_STARS_IO_H

#include "io_properties.h"
#include "stars_part.h"

/**
 * @brief Specifies which s-particle fields to read from a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void stars_read_particles(struct spart *sparts,
                                        struct io_props *list,
                                        int *num_fields) {

  /* Say how much we want to read */
  *num_fields = 7;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, sparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, sparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, sparts, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, sparts, h);
  list[5] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass_init);
  list[6] = io_make_input_field("StellarFormationTime", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_NO_UNITS, sparts, birth_time);
}

INLINE static void convert_spart_pos(const struct engine *e,
                                     const struct spart *sp, double *ret) {

  const struct space *s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(sp->x[0] - s->pos_dithering[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(sp->x[1] - s->pos_dithering[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(sp->x[2] - s->pos_dithering[2], 0.0, s->dim[2]);
  } else {
    ret[0] = sp->x[0];
    ret[1] = sp->x[1];
    ret[2] = sp->x[2];
  }
}

INLINE static void convert_spart_vel(const struct engine *e,
                                     const struct spart *sp, float *ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology *cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, sp->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, sp->time_bin);

  /* Get time-step since the last kick */
  float dt_kick_grav;
  if (with_cosmology) {
    dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_grav -=
        cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
  } else {
    dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
  }

  /* Extrapolate the velocites to the current time */
  const struct gpart *gp = sp->gpart;
  ret[0] = gp->v_full[0] + gp->a_grav[0] * dt_kick_grav;
  ret[1] = gp->v_full[1] + gp->a_grav[1] * dt_kick_grav;
  ret[2] = gp->v_full[2] + gp->a_grav[2] * dt_kick_grav;

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

/**
 * @brief Specifies which s-particle fields to write to a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * @param with_cosmology Are we running a cosmological simulation?
 */
INLINE static void stars_write_particles(const struct spart *sparts,
                                         struct io_props *list, int *num_fields,
                                         const int with_cosmology) {

  /* Say how much we want to write */
  *num_fields = 38;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_spart(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, sparts,
      convert_spart_pos, "Co-moving position of the particles");

  list[1] = io_make_output_field_convert_spart(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, sparts, convert_spart_vel,
      "Peculiar velocities of the particles. This is a * dx/dt where x is the "
      "co-moving position of the particles.");

  list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 sparts, mass,
                                 "Masses of the particles at the current point "
                                 "in time (i.e. after stellar losses)");

  list[3] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           sparts, id, "Unique ID of the particles");

  list[4] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sparts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[5] = io_make_output_field("InitialMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 sparts, mass_init,
                                 "Masses of the star particles at birth time");

  if (with_cosmology) {
    list[6] = io_make_output_field(
        "BirthScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        birth_scale_factor, "Scale-factors at which the stars were born");
  } else {
    list[6] = io_make_output_field("BirthTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f,
                                   sparts, birth_time,
                                   "Times at which the stars were born");
  }

  list[7] = io_make_output_field(
      "FeedbackEnergyFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      SNII_f_E,
      "Fractions of the canonical feedback energy that was used for the stars' "
      "SNII feedback events");

  list[8] = io_make_output_field(
      "HIIregions_last_rebuild", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      HIIregion_last_rebuild,
      "Age of star in Myr when HII region was last rebuilt");

  list[9] = io_make_output_field(
      "HIIregions_mass_to_ionize", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      HIIregion_mass_to_ionize,
      "Masses of the HII regions at the current point "
      "in time");

  list[10] = io_make_output_field("Timestep", FLOAT, 1, UNIT_CONV_TIME, 0.f,
                                  sparts, star_timestep,
                                  "Current timestep of the star particle");

  list[11] = io_make_output_field(
      "HIIregions_mass_in_kernel", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      HIIregion_mass_in_kernel,
      "Masses in kernels at time of HII region formation");

  //list[12] = io_make_output_field(
  //    "GCs_IDs", INT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_NO_UNITS, 0.f, sparts,
  //    clusters.id,
  //    "IDs of clusters tagged to particle (only unique to star particle)");

  list[12] = io_make_output_field(
      "GCs_Masses", FLOAT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_MASS, 0.f, sparts,
      clusters.mass, "Masses of clusters tagged to particle");

  list[13] = io_make_output_field(
      "GCs_InitialMasses", FLOAT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_MASS, 0.f,
      sparts, clusters.initial_mass,
      "Initial masses of clusters tagged to particle");

  if (with_cosmology) {
    list[14] = io_make_output_field(
        "GCs_DisruptionScaleFactors", FLOAT, MOSAICS_MAX_CLUSTERS, 
        UNIT_CONV_NO_UNITS, 0.f, sparts, clusters.disruption_time,
        "Scale-factors at which clusters became disrupted (-1 if surviving)");
  } else {
    list[14] = io_make_output_field(
        "GCs_DisruptionTimes", FLOAT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_TIME, 0.f,
        sparts, clusters.disruption_time,
        "Times at which clusters became disrupted (-1 if surviving)");
  }

  //list[15] = io_make_output_field(
  //    "GCs_MassLossEvap", FLOAT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_MASS, 0.f, sparts,
  //    clusters.dmevap, "Cluster mass lost to evaporation");

  list[15] = io_make_output_field(
      "GCs_MassLossShocks", FLOAT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_MASS, 0.f, sparts,
      clusters.dmshock, "Cluster mass lost to tidal shocks");

  list[16] = io_make_output_field(
      "NumberOfClusters", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      num_clusters, "Surviving number of subgrid star clusters");

  list[17] = io_make_output_field(
      "InitialNumberOfClusters", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      initial_num_clusters, "Number of star clusters tried to form in total");

  list[18] = io_make_output_field(
      "InitialNumberOfClustersEvoLim", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      initial_num_clusters_evo, 
      "Number of star clusters formed above evolution mass limit");

  list[19] = io_make_output_field(
      "InitialClusterMassTotal", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      initial_cluster_mass_total, "Sum of initial star cluster masses");

  list[20] = io_make_output_field(
      "InitialClusterMassEvoLim", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      initial_cluster_mass_evo,
      "Sum of initial star cluster masses above evolution mass limit");

  list[21] = io_make_output_field(
      "ClusterMassTotal", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      cluster_mass_total, "Sum of surviving star cluster masses");

  list[22] = io_make_output_field(
      "SubgridFieldMass", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      field_mass, "Subgrid field mass component of star particle");

  list[23] = io_make_output_field(
      "TidalTensors", FLOAT, 18, UNIT_CONV_TIDAL_TENSOR, -3.f, sparts,
      tidal_tensor,
      "Second derivative of gravitational potential at the last 3 timesteps. "
      "We store the upper components of the symmeteric matrix (i.e. 6 values "
      "per tensor: xx, xy, xz, yy, yz, zz)");

  list[24] = io_make_output_field(
      "ClusterFormationEfficiency", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      CFE, "Fraction of star formation which occurs in bound star clusters");

  list[25] = io_make_output_field(
      "ClusterTruncationMass", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts, Mcstar,
      "Exponential truncation of star cluster initial mass function (Mcstar)");

  list[26] = io_make_output_field(
      "BirthStarVelocityDispersions", FLOAT, 1, UNIT_CONV_VELOCITY_SQUARED, 
      0.f, sparts, stars_sigma_v2,
      "Local velocity dispersions (3D) squared of stars at the time of star "
      "formation within the gas smoothing length");

  list[27] = io_make_output_field(
      "BirthGasFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts, fgas,
      "Local gas fraction (within gas smoothing length) at time of star "
      "formation");

  list[28] = io_make_output_field(
      "BirthKernelStarCount", INT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts, scount,
      "Number of stars within gas smoothing length at time of star formation");

  //TODO Birth properties should not have a factors
  list[29] = io_make_output_field(
      "BirthEpicyclicFrequency", FLOAT, 1, UNIT_CONV_FREQUENCY, -1.5f, sparts, 
      kappa_birth, "Epicyclic frequency at time of star formation");

  //TODO Birth properties should not have a factors
  list[30] = io_make_output_field(
      "BirthCircularFrequency", FLOAT, 1, UNIT_CONV_FREQUENCY, -1.5f, sparts, 
      Omega_birth, "Circular frequency at time of star formation");

  list[31] = io_make_output_field(
      "BirthToomreMass", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      Toomre_mass, "Local Toomre mass at formation");

  list[32] = io_make_output_field(
      "ToomreCollapseFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      frac_collapse, "Fraction of Toomre mass which can collapse to a GMC");

  // TODO just temporary
  list[33] = io_make_output_field(
      "BirthPressures", FLOAT, 1, UNIT_CONV_PRESSURE, 0.f, sparts,
      birth_pressure,
      "Physical pressures at the time of birth of the gas particles that "
      "turned into stars (note that "
      "we store the physical pressure at the birth redshift, no conversion is "
      "needed)");

  // TODO temporary?
  list[34] = io_make_output_field(
      "BirthSoundSpeed", FLOAT, 1, UNIT_CONV_SPEED, 0.f, sparts,
      sound_speed_subgrid,
      "Hydro or subgrid sound speed, depending on cooling model");

  // TODO just temporary
  list[35] = io_make_output_field(
      "Potentials", FLOAT, 1, UNIT_CONV_POTENTIAL, -1.f, sparts, potential,
      "Co-moving gravitational potential at position of the particles");

  list[36] = io_make_output_field(
      "GCs_MassLossEvap", FLOAT, MOSAICS_MAX_CLUSTERS, UNIT_CONV_MASS, 0.f, sparts,
      clusters.dmevap, "Cluster mass lost to evaporation");

  list[37] = io_make_output_field(
      "BirthToomreLength", FLOAT, 1, UNIT_CONV_LENGTH, 0.f, sparts, 
      birth_toomre_length, "Toomre length at time of star formation");
}

/**
 * @brief Initialize the global properties of the stars scheme.
 *
 * By default, takes the values provided by the hydro.
 *
 * @param sp The #stars_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param p The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void stars_props_init(struct stars_props *sp,
                                    const struct phys_const *phys_const,
                                    const struct unit_system *us,
                                    struct swift_params *params,
                                    const struct hydro_props *p,
                                    const struct cosmology *cosmo) {

  /* Kernel properties */
  sp->eta_neighbours = parser_get_opt_param_float(
      params, "Stars:resolution_eta", p->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  sp->h_tolerance =
      parser_get_opt_param_float(params, "Stars:h_tolerance", p->h_tolerance);

  /* Get derived properties */
  sp->target_neighbours = pow_dimension(sp->eta_neighbours) * kernel_norm;
  const float delta_eta = sp->eta_neighbours * (1.f + sp->h_tolerance);
  sp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(sp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  sp->max_smoothing_iterations = parser_get_opt_param_int(
      params, "Stars:max_ghost_iterations", p->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "Stars:max_volume_change", -1);
  if (max_volume_change == -1)
    sp->log_max_h_change = p->log_max_h_change;
  else
    sp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Do we want to overwrite the stars' birth time? */
  sp->overwrite_birth_time =
      parser_get_opt_param_int(params, "Stars:overwrite_birth_time", 0);

  /* Read birth time to set all stars in ICs */
  if (sp->overwrite_birth_time) {
    sp->spart_first_init_birth_time =
        parser_get_param_float(params, "Stars:birth_time");
  }

  /* Some useful conversion values ------------------------------------------ */

  //TODO I'm not sure these should be here because of scale factor dependencies?
  sp->density_to_kgm3 =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) * 1e3;

  sp->velocity_to_ms =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) * 0.01;

  sp->time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  sp->tidal_tensor_to_cgs = 
      units_cgs_conversion_factor(us, UNIT_CONV_TIDAL_TENSOR);

  const double Msun_cgs = phys_const->const_solar_mass *
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  sp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  sp->solar_mass_to_mass = 1. / sp->mass_to_solar_mass;

  /* MOSAICS parameters ----------------------------------------------------- */

  /* Use the subgrid turbulent velocity dispersion for CFE */
  sp->subgrid_gas_vel_disp = parser_get_opt_param_int(params, 
      "Stars:use_subgrid_velocity_dispersion", 0);

  /* Use a power-law mass function (default Schechter) */
  sp->power_law_clMF =
      parser_get_opt_param_int(params, "Stars:power_law_clMF", 0);

  /* King parameter */
  sp->W0 = parser_get_opt_param_float(params, "Stars:clusters_King_W0", 5.0);

  /* Half-mass radii */
  sp->rh = parser_get_opt_param_float(params, "Stars:clusters_rh", 4.0);

  /* Disruption timescale */
  sp->fixed_t0evap = 
      parser_get_opt_param_float(params, "Stars:fixed_evap_t0", -1.f);

  /* Evaporation on by default */
  sp->cluster_evap = 
      parser_get_opt_param_int(params, "Stars:clusters_evaporation", 1);

  /* Tidal shocks on by default */
  sp->cluster_shocks = 
      parser_get_opt_param_int(params, "Stars:clusters_tidal_shocks", 1);

  /* Tidal shocks on by default */
  sp->spitzer_evap_term = 
      parser_get_opt_param_int(params, "Stars:Spitzer_evap_term", 1);

  /* Calculate smoothed Omega and kappa over stellar neighbours? */
  sp->smoothed_orbit_frequencies = 
      parser_get_opt_param_int(params, "Stars:smoothed_orbit_frequencies", 0);

  /* Use Omega^2 = -lambda_2, otherwise sum(-lambda/3) */
  sp->Omega_is_lambda2 = 
      parser_get_opt_param_int(params, "Stars:Omega_is_lambda2", 0);

  /* Maximum Toomre length for tensor softening, in kpc */
  sp->max_toomre_length = parser_get_opt_param_float(params,
      "Stars:max_toomre_length", 2.);

  /* kpc in internal units */
  sp->max_toomre_length *= 1000. * phys_const->const_parsec;

  /* Parameters of initial cluster mass function ---------------------------- */

  /* Use a power-law mass function (default Schechter) */
  sp->power_law_clMF =
      parser_get_opt_param_int(params, "Stars:power_law_clMF", 0);

  /* Integrated star formation efficiency in collapse of GMC (for Mcstar) */
  sp->SFE = parser_get_opt_param_float(params, "Stars:GMC_SFE", 0.1);

  /* Cluster mass function minimum (Msun) */
  sp->clMF_min = parser_get_opt_param_float(params, "Stars:clMF_min", 100.0);
  sp->clMF_min *= sp->solar_mass_to_mass;

  /* Initial lowest cluster mass to evolve (Msun) */
  sp->clMF_min_evolve = 
      parser_get_opt_param_float(params, "Stars:clMF_min_evolve", 5000.0);
  sp->clMF_min_evolve *= sp->solar_mass_to_mass;

  /* Cluster mass function maximum (Msun) */
  sp->clMF_max = parser_get_opt_param_float(params, "Stars:clMF_max", 1.0e9);
  sp->clMF_max *= sp->solar_mass_to_mass;

  /* Cluster mass function power-law index M^(-alpha) */
  sp->clMF_slope = parser_get_opt_param_float(params, "Stars:clMF_slope", 2.0);

  /* Use the subgrid turbulent velocity dispersion for CFE */
  sp->subgrid_gas_vel_disp = 
      parser_get_opt_param_int(params, "Stars:use_subgrid_velocity_dispersion", 0);

  /* Value of fixed cluster formation efficiency */
  sp->FixedCFE = parser_get_opt_param_float(params, "Stars:fixed_CFE", -1.0);

  /* Use the subgrid turbulent velocity dispersion for CFE */
  sp->subgrid_gas_vel_disp = 
      parser_get_opt_param_int(params, "Stars:use_subgrid_velocity_dispersion", 0);

  /* Sound speed of cold ISM (m/s) */
  sp->Fixedcs = parser_get_opt_param_float(params, "Stars:fixed_sound_speed", -1.0);

  /* Value of fixed cluster formation efficiency */
  sp->FixedCFE = parser_get_opt_param_float(params, "Stars:fixed_CFE", -1.0);

  /* Use the subgrid turbulent velocity dispersion for CFE */
  sp->subgrid_gas_vel_disp = 
      parser_get_opt_param_int(params, "Stars:use_subgrid_velocity_dispersion", 0);

  /* Sound speed of cold ISM (m/s) */
  sp->Fixedcs = parser_get_opt_param_float(params, "Stars:fixed_sound_speed", -1.0);

  /* star formation law. 0: Elmegreen 02; 1: Krumholz & McKee 05 */
  sp->sflaw = parser_get_opt_param_int(params, "Stars:CFE_sflaw", 0);

  /* GMC virial parameter */
  sp->qvir = parser_get_opt_param_float(params, "Stars:CFE_qvir", 1.3);

  /* time of first SN (Myr) */
  sp->tsn = parser_get_opt_param_float(params, "Stars:CFE_tsn", 3.0);

  /* time of determining the CFE (Myr) */
  sp->tview = parser_get_opt_param_float(params, "Stars:CFE_tview", 10.0);

  /* GMC surface density (Msun/pc^2) */
  sp->surfGMC = parser_get_opt_param_float(params, "Stars:CFE_surfGMC", 100.0);

  /* maximum (protostellar core) SFE */
  sp->ecore = parser_get_opt_param_float(params, "Stars:CFE_ecore", 0.5);

  /* turbulent/magnetic pressure ratio */
  sp->beta0 = parser_get_opt_param_float(params, "Stars:CFE_beta0", 1.0e10);

  /* SN/radiative feedback mode */
  sp->radfb = parser_get_opt_param_int(params, "Stars:CFE_radfb", 0);
}

/**
 * @brief Print the global properties of the stars scheme.
 *
 * @param sp The #stars_props.
 */
INLINE static void stars_props_print(const struct stars_props *sp) {

  /* Now stars */
  message("Stars kernel: %s with eta=%f (%.2f neighbours).", kernel_name,
          sp->eta_neighbours, sp->target_neighbours);

  message("Stars relative tolerance in h: %.5f (+/- %.4f neighbours).",
          sp->h_tolerance, sp->delta_neighbours);

  message(
      "Stars integration: Max change of volume: %.2f "
      "(max|dlog(h)/dt|=%f).",
      pow_dimension(expf(sp->log_max_h_change)), sp->log_max_h_change);

  message("Maximal iterations in ghost task set to %d",
          sp->max_smoothing_iterations);

  if (sp->overwrite_birth_time)
    message("Stars' birth time read from the ICs will be overwritten to %f",
            sp->spart_first_init_birth_time);

  message("MOSAICS maximum clusters: %d", MOSAICS_MAX_CLUSTERS);

  if (sp->FixedCFE > 0)
    message("MOSAICS using fixed CFE: %f", sp->FixedCFE);

  if (sp->power_law_clMF)
    message("MOSAICS using power-law mass function");
}

#if defined(HAVE_HDF5)
INLINE static void stars_props_print_snapshot(hid_t h_grpstars,
                                              const struct stars_props *sp) {

  io_write_attribute_s(h_grpstars, "Kernel function", kernel_name);
  io_write_attribute_f(h_grpstars, "Kernel target N_ngb",
                       sp->target_neighbours);
  io_write_attribute_f(h_grpstars, "Kernel delta N_ngb", sp->delta_neighbours);
  io_write_attribute_f(h_grpstars, "Kernel eta", sp->eta_neighbours);
  io_write_attribute_f(h_grpstars, "Smoothing length tolerance",
                       sp->h_tolerance);
  io_write_attribute_f(h_grpstars, "Volume log(max(delta h))",
                       sp->log_max_h_change);
  io_write_attribute_f(h_grpstars, "Volume max change time-step",
                       pow_dimension(expf(sp->log_max_h_change)));
  io_write_attribute_i(h_grpstars, "Max ghost iterations",
                       sp->max_smoothing_iterations);
  io_write_attribute_i(h_grpstars, "MOSAICS_MAX_CLUSTERS",
                       MOSAICS_MAX_CLUSTERS);
}
#endif

/**
 * @brief Write a #stars_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
INLINE static void stars_props_struct_dump(const struct stars_props *p,
                                           FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct stars_props), 1, stream,
                       "starsprops", "stars props");
}

/**
 * @brief Restore a stars_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
INLINE static void stars_props_struct_restore(const struct stars_props *p,
                                              FILE *stream) {
  restart_read_blocks((void *)p, sizeof(struct stars_props), 1, stream, NULL,
                      "stars props");
}

#endif /* SWIFT_MOSAICS_STAR_IO_H */
