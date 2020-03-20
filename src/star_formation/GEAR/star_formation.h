/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_STAR_FORMATION_H
#define SWIFT_GEAR_STAR_FORMATION_H

/* Local includes */
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "random.h"
#include "star_formation_struct.h"
#include "units.h"

/**
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * Use the star formation criterion given by eq. 3 in Revaz & Jablonka 2018.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 *
 */
INLINE static int star_formation_is_star_forming(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const struct entropy_floor_properties* restrict entropy_floor) {

  /* Check if collapsing particles */
  if (xp->sf_data.div_v > 0) {
    return 0;
  }

  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, p, xp);

  const float temperature_max = starform->maximal_temperature;

  /* Check the temperature criterion */
  if (temperature > temperature_max) {
    return 0;
  }

  /* Get the required variables */
  //const float sigma2 = p->pressure_floor_data.sigma2 * cosmo->a * cosmo->a;
  const float n_jeans_2_3 = starform->n_jeans_2_3;

  const float h = p->h * kernel_gamma * cosmo->a;
  const float density = hydro_get_physical_density(p, cosmo);

  // TODO use GRACKLE */
  const float mu = hydro_props->mu_neutral;

  /* Compute the density criterion */
  const float coef =
      M_PI_4 / (phys_const->const_newton_G * n_jeans_2_3 * h * h);
  const float density_criterion =
      coef * (hydro_gamma * phys_const->const_boltzmann_k * temperature /
	      (mu * phys_const->const_proton_mass));// +
  //sigma2);

  /* Check the density criterion */
  return density > density_criterion;
}

/**
 * @brief Compute the star-formation rate of a given particle.
 *
 * Nothing to do here. Everything is done in
 * #star_formation_should_convert_to_star.
 *
 * @param p #part.
 * @param xp the #xpart.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* starform, const struct phys_const* phys_const,
    const struct hydro_props* hydro_props, const struct cosmology* cosmo,
    const double dt_star) {}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.

 * Compute the star formation rate from eq. 4 in Revaz & Jablonka 2012.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    struct part* p, struct xpart* xp, const struct star_formation* starform,
    const struct engine* e, const double dt_star) {

  const struct phys_const* phys_const = e->physical_constants;
  const struct cosmology* cosmo = e->cosmology;

  /* Check that we are running a full time step */
  if (dt_star == 0.) {
    return 0;
  }

  /* Get a few variables */
  const float G = phys_const->const_newton_G;
  const float density = hydro_get_physical_density(p, cosmo);

  /* Get the mass of the future possible star */
  const float mass_gas = hydro_get_mass(p);

  /* Compute the probability */
  const float inv_free_fall_time =
      sqrtf(density * 32.f * G * 0.33333333f * M_1_PI);
  float prob = 1.f - expf(-starform->star_formation_efficiency *
                          inv_free_fall_time * dt_star);

  /* Add the mass factor */
  if (starform->n_stars_per_part != 1) {
    const float min_mass =
        starform->mass_stars * starform->min_mass_frac_plus_one;
    const float mass_star =
        mass_gas > min_mass ? starform->mass_stars : mass_gas;

    prob *= mass_gas / mass_star;
  }

  /* Roll the dice... */
  const float random_number =
      random_unit_interval(p->id, e->ti_current, random_number_star_formation);

  /* Can we form a star? */
  return random_number < prob;
}

/**
 * @brief Decides whether a new particle should be created or if the hydro
 * particle needs to be transformed.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param starform The properties of the star formation model.
 *
 * @return 1 if a new spart needs to be created.
 */
INLINE static int star_formation_should_spawn_spart(
    struct part* p, struct xpart* xp, const struct star_formation* starform) {

  /* Check if we are splitting the particles or not */
  if (starform->n_stars_per_part == 1) {
    return 0;
  }

  const float mass_min =
      starform->min_mass_frac_plus_one * starform->mass_stars;
  return hydro_get_mass(p) > mass_min;
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct star_formation* starform, const int with_cosmology) {}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle.
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sp the new created star particle with its properties.
 * @param starform the star formation law properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 * @param convert_part Did we convert a part (or spawned one)?
 */
INLINE static void star_formation_copy_properties(
    struct part* p, const struct xpart* xp, struct spart* sp,
    const struct engine* e, const struct star_formation* starform,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const int convert_part) {

  /* Store the current mass */
  const float mass_gas = hydro_get_mass(p);
  if (!convert_part) {
    const float mass_star = starform->mass_stars;
    const float new_mass_gas = mass_gas - mass_star;

    /* Update the spart */
    sp->mass = mass_star;
    sp->gpart->mass = mass_star;

    /* Update the part */
    hydro_set_mass(p, new_mass_gas);
    p->gpart->mass = new_mass_gas;
  } else {
    sp->mass = mass_gas;
  }
  sp->birth.mass = sp->mass;

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Store the tracers data */
  sp->tracers_data = xp->tracers_data;

  /* Store the birth density in the star particle */
  sp->birth.density = hydro_get_physical_density(p, cosmo);

  /* Store the birth temperature*/
  sp->birth.temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                  cosmo, cooling, p, xp);

  /* Copy the chemistry properties */
  chemistry_copy_star_formation_properties(p, xp, sp);

  /* Copy the progenitor id */
  sp->prog_id = p->id;
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 */
INLINE static void starformation_print_backend(
    const struct star_formation* starform) {
  message("Star formation law is 'GEAR'");
}

/**
 * @brief Finishes the density calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 * @param xp The extra particle to act upon
 * @param cd The global star_formation information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_formation_end_density(
    struct part* restrict p, struct xpart* restrict xp,
    const struct star_formation* cd, const struct cosmology* cosmo) {

#ifdef SPHENIX_SPH
  /* Copy the velocity divergence */
  xp->sf_data.div_v = p->viscosity.div_v;
#else
  /* Copy the velocity divergence */
  xp->sf_data.div_v = p->density.div_v;
#endif
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(struct part* restrict p,
                                      struct xpart* restrict xp,
                                      const struct star_formation* cd,
                                      const struct cosmology* cosmo) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * state to start the density loop.
 *
 * Nothing to do here.
 *
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    struct part* restrict p, const struct star_formation* data) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state at the beginning of the simulation after the ICs have been read.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global star_formation information used for this run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_part(const struct phys_const* restrict phys_const,
                               const struct unit_system* restrict us,
                               const struct cosmology* restrict cosmo,
                               const struct star_formation* data,
                               const struct part* restrict p,
                               struct xpart* restrict xp) {}

/**
 * @brief Split the star formation content of a particle into n pieces
 *
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void star_formation_split_part(
    struct part* p, struct xpart* xp, const double n) {
  error("Loic: to be implemented");
}

/**
 * @brief Deal with the case where no spart are available for star formation.
 *
 * @param e The #engine.
 * @param p The #part.
 * @param xp The #xpart.
 */
__attribute__((always_inline)) INLINE static void
star_formation_no_spart_available(const struct engine* e, const struct part* p,
                                  const struct xpart* xp) {
  error(
      "Failed to get a new particle. Please increase "
      "Scheduler:cell_extra_sparts "
      "or Scheduler:cell_extra_gparts");
}

/**
 * @brief Compute some information for the star formation model based
 * on all the particles that were read in.
 *
 * This is called once on start-up of the code.
 *
 * @param star_form The #star_formation structure.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_stats(struct star_formation* star_form,
                                const struct engine* e) {

  const struct space* s = e->s;
  double avg_mass = 0;

  /* Sum the mass over all the particles. */
  for (size_t i = 0; i < s->nr_parts; i++) {
    avg_mass += hydro_get_mass(&s->parts[i]);
  }

#ifdef WITH_MPI
  /* Compute the contribution from the other nodes. */
  MPI_Allreduce(MPI_IN_PLACE, &avg_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  star_form->mass_stars =
      avg_mass / (e->total_nr_parts * star_form->n_stars_per_part);

  if (e->nodeID == 0) {
    message("Mass new stars: %g", star_form->mass_stars);
  }
}

#endif /* SWIFT_GEAR_STAR_FORMATION_H */
