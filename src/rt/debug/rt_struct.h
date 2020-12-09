/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STRUCT_DEBUG_H
#define SWIFT_RT_STRUCT_DEBUG_H

/**
 * @file src/rt/debug/rt_struct.h
 * @brief Main header file for the debug radiative transfer struct.
 */

/* Additional RT data in hydro particle struct */
struct rt_part_data {

  /* data to store during entire run */
  int calls_tot; /* total number of calls to this particle during entire run */

  /* data to store during one time step */
  int calls_per_step;        /* calls per time step to this particle */
  int iact_stars_inject;     /* how many stars this particle interacted with */
  int calls_self_inject;     /* calls from self inject tasks */
  int calls_pair_inject;     /* calls from pair inject tasks */
  int calls_iact_gradient;   /* calls from gradient interaction loop */
  int calls_iact_gradient_sym;   /* calls from gradient interaction loop */
  int calls_iact_gradient_nonsym;  /* calls from gradient interaction loop */
  int calls_iact_transport;  /* calls from transport interaction loop */
  int calls_iact_transport_sym;  /* calls from transport interaction loop */
  int calls_iact_transport_nonsym; /* calls from transport interaction loop */
  int photon_number_updated; /* calls from ghost1 tasks */

  /* temp for debugging */
  long long neigh_iact_grad[400];
  long long neigh_cell_iact_grad[400];
  int neigh_iact_grad_free;
  long long neigh_iact_transp[400];
  long long neigh_cell_iact_transp[400];
  int neigh_iact_transp_free;
  long long this_cell_grad;
  long long this_cell_transport;


  /* temp for debugging */
  long long hydro_neigh_iact_grad[400];
  long long hydro_neigh_cell_iact_grad[400];
  int hydro_neigh_iact_grad_free;
  long long hydro_neigh_iact_transp[400];
  long long hydro_neigh_cell_iact_transp[400];
  int hydro_neigh_iact_transp_free;
  long long hydro_this_cell_grad;
  long long hydro_this_cell_transport;

  /* temp for debugging */
  int calls_hydro_iact_gradient;
  int calls_hydro_iact_gradient_sym;
  int calls_hydro_iact_gradient_nonsym;
  int calls_hydro_iact_force;
  int calls_hydro_iact_force_sym;
  int calls_hydro_iact_force_nonsym;


  float h_grad;
  float h_transport;
  float h_hydro_grad;
  float h_force;

};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* data to store during entire run */
  int calls_tot; /* total number of calls to this particle during entire run */

  /* data to store during one time step */
  int calls_per_step;    /* calls per time step to this particle */
  int iact_hydro_inject; /* how many hydro particles this particle interacted
                            with */
  int calls_self_inject; /* calls from self inject tasks */
  int calls_pair_inject; /* calls from pair inject tasks */
  int emission_rate_set; /* stellar photon emisison rate has been computed */
};

#endif /* SWIFT_RT_STRUCT_DEBUG_H */
