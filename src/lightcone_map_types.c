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

/* Local includes */
#include "black_holes.h"
#include "engine.h"
#include "gravity.h"
#include "hydro.h"
#include "lightcone_map.h"
#include "part.h"
#include "stars.h"

/* This object's header */
#include "lightcone_map_types.h"


static double angular_smoothing_scale(const double *pos, const double hsml) {
  
  /* Compute distance to particle */
  double dist = 0;
  for(int i=0; i<3; i+=1)
    dist += pos[i]*pos[i];
  dist = sqrt(dist);
  
  /* Avoid trig call for small angles (accurate to about 0.3%) */
  if(dist > 10.0*hsml)
    return hsml/dist;
  else
    return atan(hsml/dist);
}

/**
 * @brief Make a healpix map of projected mass in each pixel
 *
 * @param map the #lightcone_map structure
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
void lightcone_map_total_mass(struct lightcone_map *map, const struct engine *e,
                              const struct gpart *gp, const double a_cross,
                              const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  /* const struct xpart *xparts = s->xparts; */ /* Currently not used */
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const double radius = angular_smoothing_scale(x_cross, p->h);
    lightcone_map_buffer_update(map, x_cross, radius, p->mass);
  } break;
  case swift_type_stars: {
    const struct spart *sp = &sparts[-gp->id_or_neg_offset];
    const double radius = 0.0;
    lightcone_map_buffer_update(map, x_cross, radius, sp->mass);
  } break;
  case swift_type_black_hole: {      
    const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
    const double radius = 0.0;
    lightcone_map_buffer_update(map, x_cross, radius, bp->mass);
  } break;
  case swift_type_dark_matter:
  case swift_type_dark_matter_background:
  case swift_type_neutrino: {
    const double radius = 0.0;
    lightcone_map_buffer_update(map, x_cross, radius, gp->mass);
  } break;
  default:
    /* Unknown type, nothing to do */
    break;
  }
}


/**
 * @brief Make a healpix map of projected gas mass in each pixel
 *
 * @param map the #lightcone_map structure
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
void lightcone_map_gas_mass(struct lightcone_map *map, const struct engine *e,
                            const struct gpart *gp, const double a_cross,
                            const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const double radius = angular_smoothing_scale(x_cross, p->h);
    lightcone_map_buffer_update(map, x_cross, radius, p->mass);
  } break;
  default:
    /* Not gas, nothing to do */
    break;
  }
}

/**
 * @brief Make a healpix map of projected neutrino mass in each pixel
 *
 * @param map the #lightcone_map structure
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
void lightcone_map_neutrino_mass(struct lightcone_map *map, const struct engine *e,
                                 const struct gpart *gp, const double a_cross,
                                 const double x_cross[3]) {
  
  switch (gp->type) {
  case swift_type_neutrino: {
    const double radius = 0.0;
    lightcone_map_buffer_update(map, x_cross, radius, gp->mass);
  } break;
  default:
    /* Not a neutrino, nothing to do */
    break;
  }
}

/**
 * @brief Make a healpix map of the compton y parameter
 *
 * @param map the #lightcone_map structure
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
void lightcone_map_compton_y(struct lightcone_map *map, const struct engine *e,
                             const struct gpart *gp, const double a_cross,
                             const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart* xp = &xparts[-gp->id_or_neg_offset];
      double y_compton = cooling_get_ycompton(phys_const, hydro_props, us, cosmo,
                                       cool_func, p, xp);
      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      double pixel_size_2 = map->pixel_area_steradians;

      double y_for_map = y_compton / (pixel_size_2 * angular_diameter_distance_2);

      const double radius = angular_smoothing_scale(x_cross, p->h);
      lightcone_map_buffer_update(map, x_cross, radius, y_for_map);
    } break;
    default:
      /* Not gas, nothing to do */
      break;
  }
}

/**
 * @brief Make a healpix map of the doppler b parameter
 *
 * @param map the #lightcone_map structure
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
void lightcone_map_doppler_b(struct lightcone_map *map, const struct engine *e,
                             const struct gpart *gp, const double a_cross,
                             const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart* xp = &xparts[-gp->id_or_neg_offset];
      double n_e = cooling_get_electron_pressure(phys_const, hydro_props, us, cosmo,
                                          cool_func, p, xp);

      double rho = hydro_get_physical_density(p, cosmo);

      double m = hydro_get_mass(p);

      const double c = phys_const->const_speed_light_c;

      const double sigma_thompson = phys_const->const_thomson_cross_section;

      double doppler_b_factor = n_e * m * sigma_thompson / (rho * c);

      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      double angular_diameter_distance = pow(angular_diameter_distance_2, 0.5);

      double radial_velocity =
          (p->v[0] * x_cross[0] * a_cross + p->v[1] * x_cross[1] * a_cross +
           p->v[2] * x_cross[2] * a_cross) /
          angular_diameter_distance;

      double pixel_size_2 = map->pixel_area_steradians;

      double b_for_map = doppler_b_factor * radial_velocity /
                  (pixel_size_2 * angular_diameter_distance_2);

      const double radius = angular_smoothing_scale(x_cross, p->h);
      lightcone_map_buffer_update(map, x_cross, radius, b_for_map);
    } break;
    default:
      /* Not gas, nothing to do */
      break;
  }
}
