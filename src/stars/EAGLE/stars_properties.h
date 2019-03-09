/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STAR_PROPERTIES_H
#define SWIFT_EAGLE_STAR_PROPERTIES_H

/* Forward declarations */
struct hydro_props;
struct unit_system;
struct swift_params;
struct phys_const;

/* Local includes */
#include "inline.h"

/* Standard headers */
#include <stdio.h>

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/**
 * @brief Contains all the constants and parameters of the stars scheme
 */
struct stars_props {

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Properties of the EAGLE feedback model */
  struct {

    /*! Temperature change to apply to the gas particles heated by feedback */
    double delta_T;

    /*! Maximal fraction of the available SNe energy to inject */
    double f_E_max;

    /*! Minimal fraction of the available SNe energy to inject */
    double f_E_min;

    /*! Width of the sigmoid used in the feedback energy fraction model. */
    double sigma;

    /*! Anchor point for the metallicity dependance of the feedback energy
     * fraction model. */
    double Z_0;

    /*! Anchor point for the density dependance of the feedback energy fraction
     * model. */
    double rho_0;

    /*! Power-law of the density dependance of the feedback energy fraction
     * model. */
    double n_n;

    /*! Energy of one super-nova in cgs units */
    double E_SNe_cgs;

    /*! Energy of one super-nova in internal units */
    double E_SNe;

    /*! Minimal stellar mass for a SNII in solar masses */
    double SNII_min_mass_Msun;

    /*! Maximal stellar mass for a SNII in solar masses */
    double SNII_max_mass_Msun;

    /*! Time between birth of star and SNII feedback in years */
    double SNII_feedback_delay_years;

    /*! Time between birth of star and SNII feedback in internal units */
    double SNII_feedback_delay;

  } feedback;
};

void stars_props_init(struct stars_props *sp,
                      const struct phys_const *phys_const,
                      const struct unit_system *us, struct swift_params *params,
                      const struct hydro_props *p);

void stars_props_print(const struct stars_props *sp);

#if defined(HAVE_HDF5)
void stars_props_print_snapshot(hid_t h_grpstars, const struct stars_props *sp);
#endif

void stars_props_struct_dump(const struct stars_props *p, FILE *stream);

void stars_props_struct_restore(const struct stars_props *p, FILE *stream);

#endif /* SWIFT_EAGLE_STAR_PROPERTIES_H */
