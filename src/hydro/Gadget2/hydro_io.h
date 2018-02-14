/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_GADGET2_HYDRO_IO_H
#define SWIFT_GADGET2_HYDRO_IO_H

#include "adiabatic_index.h"
#include "hydro.h"
#include "io_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void hydro_read_particles(struct part* parts,
                                        struct io_props* list,
                                        int* num_fields) {

  *num_fields = 8;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                parts, mass);
  list[3] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[4] = io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts, entropy);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, rho);
}

INLINE static void convert_part_u(const struct engine* e, const struct part* p,
                                  const struct xpart* xp, float* ret) {

  ret[0] = hydro_get_comoving_internal_energy(p);
}

INLINE static void convert_part_P(const struct engine* e, const struct part* p,
                                  const struct xpart* xp, float* ret) {

  ret[0] = hydro_get_comoving_pressure(p);
}

INLINE static void convert_part_pos(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, double* ret) {

  if (e->s->periodic) {
    ret[0] = box_wrap(p->x[0], 0.0, e->s->dim[0]);
    ret[1] = box_wrap(p->x[1], 0.0, e->s->dim[1]);
    ret[2] = box_wrap(p->x[2], 0.0, e->s->dim[2]);
  } else {
    ret[0] = p->x[0];
    ret[1] = p->x[1];
    ret[2] = p->x[2];
  }
}

INLINE static void convert_part_vel(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, p->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, p->time_bin);

  /* Get time-step since the last kick */
  float dt_kick_grav, dt_kick_hydro;
  if (with_cosmology) {
    dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_grav -=
        cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
    dt_kick_hydro = cosmology_get_hydro_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_hydro -=
        cosmology_get_hydro_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
  } else {
    dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
    dt_kick_hydro = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
  }

  /* Extrapolate the velocites to the current time */
  hydro_get_drifted_velocities(p, xp, dt_kick_hydro, dt_kick_grav, ret);

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

INLINE static void convert_part_potential(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, float* ret) {

  if (p->gpart != NULL)
    ret[0] = gravity_get_comoving_potential(p->gpart);
  else
    ret[0] = 0.f;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void hydro_write_particles(const struct part* parts,
                                         const struct xpart* xparts,
                                         struct io_props* list,
                                         int* num_fields) {

  *num_fields = 10;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_part("Coordinates", DOUBLE, 3,
                                              UNIT_CONV_LENGTH, parts, xparts,
                                              convert_part_pos);
  list[1] = io_make_output_field_convert_part(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, parts, xparts, convert_part_vel);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, parts, mass);
  list[3] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 parts, h);
  list[4] = io_make_output_field(
      "Entropy", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS, parts, entropy);
  list[5] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);
  list[6] =
      io_make_output_field("Density", FLOAT, 1, UNIT_CONV_DENSITY, parts, rho);
  list[7] = io_make_output_field_convert_part("InternalEnergy", FLOAT, 1,
                                              UNIT_CONV_ENERGY_PER_UNIT_MASS,
                                              parts, xparts, convert_part_u);
  list[8] = io_make_output_field_convert_part(
      "Pressure", FLOAT, 1, UNIT_CONV_PRESSURE, parts, xparts, convert_part_P);

  list[9] = io_make_output_field_convert_part("Potential", FLOAT, 1,
                                              UNIT_CONV_POTENTIAL, parts,
                                              xparts, convert_part_potential);

#ifdef DEBUG_INTERACTIONS_SPH

  list += *num_fields;
  *num_fields += 4;

  list[0] = io_make_output_field("Num_ngb_density", INT, 1, UNIT_CONV_NO_UNITS,
                                 parts, num_ngb_density);
  list[1] = io_make_output_field("Num_ngb_force", INT, 1, UNIT_CONV_NO_UNITS,
                                 parts, num_ngb_force);
  list[2] =
      io_make_output_field("Ids_ngb_density", LONGLONG, MAX_NUM_OF_NEIGHBOURS,
                           UNIT_CONV_NO_UNITS, parts, ids_ngbs_density);
  list[3] =
      io_make_output_field("Ids_ngb_force", LONGLONG, MAX_NUM_OF_NEIGHBOURS,
                           UNIT_CONV_NO_UNITS, parts, ids_ngbs_force);

#endif
}


/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * 
 * In this version, we only want the ids and the offset.
 */
void hydro_write_index(struct part* parts, struct io_props* list,
		       int* num_fields) {

#ifdef WITH_LOGGER
  *num_fields = 2;

  /* List what we want to write */
  list[0] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);

  list[1] = io_make_output_field("Offset", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, last_offset);
#else
  error("Cannot write index without logger");
#endif
  
}


/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void hydro_write_flavour(hid_t h_grpsph) {

  /* Viscosity and thermal conduction */
  io_write_attribute_s(h_grpsph, "Thermal Conductivity Model",
                       "(No treatment) as in Springel (2005)");
  io_write_attribute_s(
      h_grpsph, "Viscosity Model",
      "as in Springel (2005), i.e. Monaghan (1992) with Balsara (1995) switch");
  io_write_attribute_f(h_grpsph, "Viscosity alpha", const_viscosity_alpha);
  io_write_attribute_f(h_grpsph, "Viscosity beta", 3.f);
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
INLINE static int writeEntropyFlag(void) { return 0; }

#endif /* SWIFT_GADGET2_HYDRO_IO_H */
