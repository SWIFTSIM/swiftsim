/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Angus Lepper (angus.lepper@ed.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "velociraptor_interface.h"

/**
 * @brief Initialise VELOCIraptor with input and output file names along with cosmological info needed to run.
 *
 * @param e The #engine.
 *
 */
void init_velociraptor(struct engine *e) {
    struct cosmoinfo cosmo_info;
    struct unitinfo unit_info;
    struct siminfo sim_info;
    
    struct unit_system vel_us;
    struct phys_const vel_const;

    /* Initialize velociraptor unit system and constants */
    units_init(&vel_us, e->parameter_file, "VelociraptorUnitSystem");
    phys_const_init(&vel_us, e->parameter_file, &vel_const);

    /* Set cosmological constants. */
    cosmo_info.atime = e->cosmology->a;
    cosmo_info.littleh = e->cosmology->h;
    cosmo_info.Omega_m = e->cosmology->Omega_m;
    cosmo_info.Omega_b = e->cosmology->Omega_b;
    cosmo_info.Omega_Lambda = e->cosmology->Omega_lambda;
    cosmo_info.Omega_cdm = e->cosmology->Omega_m - e->cosmology->Omega_b;
    cosmo_info.w_de = e->cosmology->w_0;

    message("Scale factor: %e", cosmo_info.atime);
    message("Little h: %e", cosmo_info.littleh);
    message("Omega_m: %e", cosmo_info.Omega_m);
    message("Omega_b: %e", cosmo_info.Omega_b);
    message("Omega_Lambda: %e", cosmo_info.Omega_Lambda);
    message("Omega_cdm: %e", cosmo_info.Omega_cdm);
    message("w_de: %e", cosmo_info.w_de);

    /* Set unit conversions. */
    unit_info.lengthtokpc = units_conversion_factor(e->internal_units, &vel_us, UNIT_CONV_LENGTH); /* 1kpc <=> 3.086e21cm */
    unit_info.velocitytokms = units_conversion_factor(e->internal_units, &vel_us, UNIT_CONV_SPEED); /* 1km/s <=> 1e5cm/s */
    unit_info.masstosolarmass = units_conversion_factor(e->internal_units, &vel_us, UNIT_CONV_MASS); /* 1M_sol <=> 1.99e33g */
    unit_info.gravity = vel_const.const_newton_G; /* TODO: G = 6.67408e-8 (cgs) */
    unit_info.hubbleunit = e->cosmology->H; /* TODO: double check this. */

    message("Length conversion factor: %e", unit_info.lengthtokpc);
    message("Velocity conversion factor: %e", unit_info.velocitytokms);
    message("Mass conversion factor: %e", unit_info.masstosolarmass);
    message("G: %e", unit_info.gravity);
    message("H: %e", unit_info.hubbleunit);

    const int nr_gparts = e->s->nr_gparts;
    
    /* Set simulation information. */
    sim_info.period = unit_info.lengthtokpc * e->s->cdim[0] * e->s->width[0]; /* Physical size of box in VELOCIraptor units (kpc). */
    sim_info.zoomhigresolutionmass = -1.0; /* Placeholder. */
    sim_info.interparticlespacing = sim_info.period / pow(nr_gparts, 1./3.); /* Placeholder. */
    sim_info.icosmologicalsim = (e->policy & engine_policy_cosmology); /* Placeholder. */

    message("Period: %e", sim_info.period);
    message("Zoom high res mass: %e", sim_info.zoomhigresolutionmass);
    message("Inter-particle spacing: %e", sim_info.interparticlespacing);
    message("Cosmological: %d", sim_info.icosmologicalsim);

    InitVelociraptor("stf_input.cfg", "stf_ouput.out", cosmo_info, unit_info, sim_info);
}

/**
 * @brief Run VELOCIraptor with current particle data.
 *
 * @param e The #engine.
 *
 */
void invoke_velociraptor(struct engine *e) {

    struct gpart *gparts = e->s->gparts;
    const int nr_gparts = e->s->nr_gparts;

    InvokeVelociraptor(nr_gparts, gparts);
}
