/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TRACERS_EAGLE_H
#define SWIFT_TRACERS_EAGLE_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "cooling.h"
#include "hydro.h"
#include "part.h"
#include "tracers_struct.h"

/**
 * @brief Update the particle tracers just after it has been initialised at the
 * start of a step.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 */
static INLINE void tracers_after_init(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after it has been drifted.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 */
static INLINE void tracers_after_drift(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after its time-step has been
 * computed.
 *
 * In EAGLE we record the highest temperature reached.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 */
static INLINE void tracers_after_timestep(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {

  /* Current temperature */
  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, p, xp);

  /* New record? */
  if (temperature > xp->tracers_data.maximum_temperature) {

    xp->tracers_data.maximum_temperature = temperature;

    if (with_cosmology) {
      xp->tracers_data.maximum_temperature_scale_factor = cosmo->a;
    } else {
      xp->tracers_data.maximum_temperature_time = time;
    }
  }

  /* Also record maximum (SF) density */
  const float density = hydro_get_physical_density(p, cosmo);
  if (p->SFR > 0 && density > xp->tracers_data.maximum_physical_density)
    xp->tracers_data.maximum_physical_density = density;

}

/**
 * @brief Initialise the tracer data at the start of a calculation.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data (containing the tracers
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 */
static INLINE void tracers_first_init_xpart(
    const struct part *p, struct xpart *xp, const struct unit_system *us,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling) {

  xp->tracers_data.maximum_temperature = -1.f;
  xp->tracers_data.maximum_temperature_time = -1.f;
  xp->tracers_data.maximum_physical_density = -1.f;
  xp->tracers_data.hit_by_SNII_feedback = 0;
  xp->tracers_data.hit_by_AGN_feedback = 0;
  xp->tracers_data.AGN_feedback_energy = 0.f;
  xp->tracers_data.SNII_birth_density = -FLT_MAX;
  xp->tracers_data.SNII_feedback_energy = 0.f;
}

/**
 * @brief Update the particles' tracer data after a stellar feedback
 * event.
 *
 * @param xp The extended particle data.
 */
static INLINE void tracers_after_feedback(
  struct xpart *xp, const double nH_birth, const double delta_energy) {

  xp->tracers_data.hit_by_SNII_feedback++;

  /* Current total birth_density for this particle */
  double gas_nH_birth = xp->tracers_data.SNII_feedback_energy *
      xp->tracers_data.SNII_birth_density;

  /* New total */
  xp->tracers_data.SNII_birth_density =
      (gas_nH_birth + nH_birth * delta_energy) /
      (xp->tracers_data.SNII_feedback_energy + delta_energy);

  xp->tracers_data.SNII_feedback_energy += delta_energy;
}

/**
 * @brief Update the particles' tracer data after an AGN feedback
 * event.
 *
 * @param xp The extended particle data.
 * @param with_cosmology Are we running with cosmology?
 * @param scale_factor The current scale-factor (if running with cosmo)
 * @param time The current time (if running without cosmo)
 * @param delta_energy Amount of energy injected in the feedback event
 * (internal physical units)
 */
static INLINE void tracers_after_black_holes_feedback(
    struct xpart *xp, const int with_cosmology, const float scale_factor,
    const double time, const double delta_energy) {

  xp->tracers_data.hit_by_AGN_feedback++;
  xp->tracers_data.AGN_feedback_energy += delta_energy;
}

/**
 * @brief Split the tracer content of a particle into n pieces
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void tracers_split_part(
    struct part *p, struct xpart *xp, const double n) {

  xp->tracers_data.AGN_feedback_energy /= n;
}

#endif /* SWIFT_TRACERS_EAGLE_H */
