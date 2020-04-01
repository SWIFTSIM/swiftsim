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
  *num_fields = 28;

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

  list[5] = io_make_output_field(
      "BirthDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts, birth_density,
      "Physical densities at the time of birth of the gas particles that "
      "turned into stars (note that "
      "we store the physical density at the birth redshift, no conversion is "
      "needed)");

  list[6] = io_make_output_field("InitialMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 sparts, mass_init,
                                 "Masses of the star particles at birth time");

  if (with_cosmology) {
    list[7] = io_make_output_field(
        "BirthScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        birth_scale_factor, "Scale-factors at which the stars were born");
  } else {
    list[7] = io_make_output_field("BirthTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f,
                                   sparts, birth_time,
                                   "Times at which the stars were born");
  }

  list[8] = io_make_output_field(
      "FeedbackEnergyFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      SNII_f_E,
      "Fractions of the canonical feedback energy that was used for the stars' "
      "SNII feedback events");

  list[9] =
      io_make_output_field("BirthTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE,
                           0.f, sparts, birth_temperature,
                           "Temperatures at the time of birth of the gas "
                           "particles that turned into stars");

  list[10] = io_make_output_field(
      "HIIregions_last_rebuild", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      HIIregion_last_rebuild,
      "Age of star in Myr when HII region was last rebuilt");

  list[11] = io_make_output_field(
      "HIIregions_mass_to_ionize", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      HIIregion_mass_to_ionize,
      "Masses of the HII regions at the current point "
      "in time");

  list[12] = io_make_output_field("Timestep", FLOAT, 1, UNIT_CONV_TIME, 0.f,
                                  sparts, star_timestep,
                                  "Current timestep of the star particle");

  list[13] = io_make_output_field(
      "HIIregions_mass_in_kernel", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      HIIregion_mass_in_kernel,
      "Masses in kernels at time of HII region formation");

  list[14] = io_make_output_field(
      "GCFlag", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, sparts, gcflag,
      "Flag denoting whether the star particle has star surviving"
      "clusters");

  list[15] = io_make_output_field(
      "TidalTensor", FLOAT, 18, UNIT_CONV_TIDAL_TENSOR, 0.f, sparts,
      tidal_tensor,
      "Second derivative of gravitational potential at the last 3 snapshots. We"
      "store the upper components of the symmeteric matrix (i.e. 6 values per"
      "tensor: xx, xy, xz, yy, yz, zz)");

  list[16] = io_make_output_field(
      "BirthPressures", FLOAT, 1, UNIT_CONV_PRESSURE, 0.f, sparts,
      birth_pressure,
      "Physical pressures at the time of birth of the gas particles that "
      "turned into stars (note that "
      "we store the physical pressure at the birth redshift, no conversion is "
      "needed)");

  list[17] = io_make_output_field(
      "ClusterFormationEfficiency", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      CFE, "Fraction of star formation which occurs in bound clusters");

  list[18] = io_make_output_field(
      "Mcstar", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts, Mcstar,
      "Exponential truncation of star cluster mass function");

  list[19] = io_make_output_field(
      "BirthSoundSpeed", FLOAT, 1, UNIT_CONV_SPEED, 0.f, sparts,
      sound_speed_subgrid,
      "Hydro or subgrid sound speed, depending on cooling model");

  list[20] = io_make_output_field(
      "GasVelocityDispersion", FLOAT, 1, UNIT_CONV_SPEED, 0.f, sparts,
      gas_vel_disp,
      "Local velocity dispersion of gas at the time of star formation");

  list[21] = io_make_output_field(
      "StellarVelocityDispersion", FLOAT, 1, UNIT_CONV_SPEED, 0.f, sparts,
      star_vel_disp,
      "Local velocity dispersion of stars at the time of star formation");

  list[22] = io_make_output_field(
      "GasFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts, fgas,
      "Local gas fraction at time of star formation");

  list[23] = io_make_output_field(
      "EpicyclicFrequency", FLOAT, 1, UNIT_CONV_FREQUENCY, 0.f, sparts, kappa,
      "Epicyclic frequency at formation");

  list[24] = io_make_output_field(
      "CircularFrequency", FLOAT, 1, UNIT_CONV_FREQUENCY, 0.f, sparts, Omega,
      "Circular frequency at formation");

  list[25] = io_make_output_field(
      "ToomreMass", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      Toomre_mass, "Local Toomre mass at formation");

  list[26] = io_make_output_field(
      "ToomreCollapseFraction", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      frac_collapse, "Fraction of Toomre mass which can collapse to a GMC");

  // TODO just temporary
  list[27] = io_make_output_field(
      "Potentials", FLOAT, 1, UNIT_CONV_POTENTIAL, -1.f, sparts, potential,
      "Co-moving gravitational potential at position of the particles");
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

  /* MOSAICS parameters ----------------------------------------------------- */

  /* Flag for forcing the tensor calculation */
  sp->calc_all_star_tensors =
      parser_get_opt_param_int(params, "Stars:calculate_all_star_tensors", 0);

  /* King parameter */
  sp->W0 = parser_get_opt_param_float(params, "Stars:King_W0", 5.0);

  /* Integrated star formation efficiency in collapse of GMC (for Mcstar) */
  sp->SFE = parser_get_opt_param_float(params, "Stars:SFE", -1.0);

  /* Use the subgrid turbulent velocity dispersion for CFE */
  sp->subgrid_gas_vel_disp = 
      parser_get_opt_param_int(params, "Stars:use_subgrid_velocity_dispersion", 0);

  /* Sound speed of cold ISM (m/s) */
  sp->Fixedcs = parser_get_opt_param_float(params, "Stars:Fixedcs", -1.0);

  /* Cluster formation efficiency parameters -------------------------------- */

  /* star formation law. 0: Elmegreen 02; 1: Krumholz & McKee 05 */
  sp->sflaw = parser_get_opt_param_int(params, "Stars:sflaw", 0);

  /* GMC virial parameter */
  sp->qvir = parser_get_opt_param_float(params, "Stars:qvir", 1.3);

  /* time of first SN (Myr) */
  sp->tsn = parser_get_opt_param_float(params, "Stars:tsn", 3.0);

  /* time of determining the CFE (Myr) */
  sp->tview = parser_get_opt_param_float(params, "Stars:tview", 10.0);

  /* GMC surface density (Msun/pc^2) */
  sp->surfGMC = parser_get_opt_param_float(params, "Stars:surfGMC", 100.0);

  /* maximum (protostellar core) SFE */
  sp->ecore = parser_get_opt_param_float(params, "Stars:ecore", 0.5);

  /* turbulent/magnetic pressure ratio */
  sp->beta0 = parser_get_opt_param_float(params, "Stars:beta0", 1.0e10);

  /* SN/radiative feedback mode */
  sp->radfb = parser_get_opt_param_int(params, "Stars:radfb", 0);

  /* Some useful conversion values ------------------------------------------ */

  sp->density_to_kgm3 =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY) * 1e3;

  sp->velocity_to_ms =
      units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY) * 0.01;

  sp->time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* TODO not sure if needed */
  // const double Msun_cgs = phys_const->const_solar_mass *
  //                        units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  // const double unit_mass_cgs = units_cgs_conversion_factor(us,
  // UNIT_CONV_MASS);  sp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
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

  message("Calculate tensors for all stars: %d", sp->calc_all_star_tensors);
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
