/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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

/* This object's header. */
#include "part.h"
#include "restart.h"
#include "star_formation.h"
#include "star_formation_io.h"
#include "units.h"

/**
 * @brief  Initialises the star formation law properties in the internal
 * unit system.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us the current internal system of units
 * @param hydro_props The propoerties of the hydro scheme.
 * @param starform the properties of the star formation law
 */
void starformation_init(struct swift_params* parameter_file,
                        const struct phys_const* phys_const,
                        const struct unit_system* us,
                        const struct hydro_props* hydro_props,
                        struct star_formation* starform) {

  starformation_init_backend(parameter_file, phys_const, us, hydro_props,
                             starform);
}

/**
 * @brief Print the properties of the star fromation law
 *
 * @param starform the star formation properties.
 */
void starformation_print(const struct star_formation* starform) {

  starformation_print_backend(starform);
}

/**
 * @brief Write an star_formation struct to the given FILE as a stream of
 * bytes.
 *
 * @param starform the star formation struct
 * @param stream the file stream
 */
void starformation_struct_dump(const struct star_formation* starform,
                               FILE* stream) {
  restart_write_blocks((void*)starform, sizeof(struct star_formation), 1,
                       stream, "starformation", "star formation");
}

/**
 * @brief Restore a star_formation struct from the given FILE as a stream of
 * bytes.
 *
 * @param starform the star formation struct
 * @param stream the file stream
 */
void starformation_struct_restore(const struct star_formation* starform,
                                  FILE* stream) {
  restart_read_blocks((void*)starform, sizeof(struct star_formation), 1, stream,
                      NULL, "star formation");
}

/**
 * @brief Move both p and sp in order to avoid a division by zero during the
 * star formation.
 *
 * @param e The #engine.
 * @param c The cell that is currently star forming.
 * @param p The #part generating a star.
 * @param xp The #xpart generating a star.
 * @param sp The new #spart.
 */
void starformation_avoid_divison_by_zero(const struct engine* e, struct cell* c,
                                         struct part* p, struct xpart* xp,
                                         struct spart* sp) {
#ifdef SWIFT_DEBUG_CHECKS
  if (p->x[0] != sp->x[0] || p->x[1] != sp->x[1] || p->x[2] != sp->x[2]) {
    error(
        "Moving particles that are not at the same location."
        " (%g, %g, %g) - (%g, %g, %g)",
        p->x[0], p->x[1], p->x[2], sp->x[0], sp->x[1], sp->x[2]);
  }
#endif

  /* Move a bit the particle in order to avoid
     division by 0.
  */
  const float max_displacement = 0.2;
  const double delta_x =
      2.f * random_unit_interval(p->id, e->ti_current,
                                 (enum random_number_type)0) -
      1.f;
  const double delta_y =
      2.f * random_unit_interval(p->id, e->ti_current,
                                 (enum random_number_type)1) -
      1.f;
  const double delta_z =
      2.f * random_unit_interval(p->id, e->ti_current,
                                 (enum random_number_type)2) -
      1.f;

  sp->x[0] += delta_x * max_displacement * p->h;
  sp->x[1] += delta_y * max_displacement * p->h;
  sp->x[2] += delta_z * max_displacement * p->h;

  /* Copy the position to the gpart */
  sp->gpart->x[0] = sp->x[0];
  sp->gpart->x[1] = sp->x[1];
  sp->gpart->x[2] = sp->x[2];

  /* Do the gas particle. */
  const double mass_ratio = sp->mass / hydro_get_mass(p);
  const double dx[3] = {mass_ratio * delta_x * max_displacement * p->h,
                        mass_ratio * delta_y * max_displacement * p->h,
                        mass_ratio * delta_z * max_displacement * p->h};

  p->x[0] -= dx[0];
  p->x[1] -= dx[1];
  p->x[2] -= dx[2];

  /* Compute offsets since last cell construction */
  xp->x_diff[0] += dx[0];
  xp->x_diff[1] += dx[1];
  xp->x_diff[1] += dx[2];
  xp->x_diff_sort[0] += dx[0];
  xp->x_diff_sort[1] += dx[1];
  xp->x_diff_sort[2] += dx[2];

  /* Copy the position to the gpart */
  p->gpart->x[0] = p->x[0];
  p->gpart->x[1] = p->x[1];
  p->gpart->x[2] = p->x[2];

  const float dx2_part = xp->x_diff[0] * xp->x_diff[0] +
                         xp->x_diff[1] * xp->x_diff[1] +
                         xp->x_diff[2] * xp->x_diff[2];
  const float dx2_sort = xp->x_diff_sort[0] * xp->x_diff_sort[0] +
                         xp->x_diff_sort[1] * xp->x_diff_sort[1] +
                         xp->x_diff_sort[2] * xp->x_diff_sort[2];

  const float dx_part = sqrtf(dx2_part);
  const float dx_sort = sqrtf(dx2_sort);

  for (struct cell* child = c; (child->hydro.dx_max_part < dx_part ||
                                child->hydro.dx_max_sort < dx_sort);
       /* NULL */) {
    child->hydro.dx_max_part = max(child->hydro.dx_max_part, dx_part);
    child->hydro.dx_max_sort = max(child->hydro.dx_max_sort, dx_sort);

    /* Should we go below? */
    if (child == child->hydro.super) break;

    /* Check that we can still go deeper */
    if (!child->split) error("Cannot go deeper");

    /* Get the correct progeny */
    const double pivot[3] = {child->loc[0] + child->width[0] / 2,
                             child->loc[1] + child->width[1] / 2,
                             child->loc[2] + child->width[2] / 2};
    const int bid = (p->x[0] >= pivot[0]) * 4 + (p->x[1] >= pivot[1]) * 2 +
                    (p->x[2] >= pivot[2]);

    c = child->progeny[bid];

    /* Check that we have a child */
    if (c == NULL) {
      error("No child found.");
    }
  }
}
