/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_STRUCT_COLIBRE_H
#define SWIFT_FEEDBACK_STRUCT_COLIBRE_H

#include "chemistry_struct.h"

/**
 * @brief Feedback fields carried by each hydro particles
 */
struct feedback_part_data {};

/**
 * @brief Feedback fields carried by each star particles
 */
struct feedback_spart_data {

  union {

    /**
     * @brief Values collected from the gas neighbours.
     */
    struct {

      /*! Inverse of normalisation factor used for the enrichment */
      float enrichment_weight_inv;

      /*! Total mass (unweighted) of neighbouring gas particles */
      float ngb_mass;

    } to_collect;

    /**
     * @brief Values to be distributed to the gas neighbours.
     *
     * WARNING: The first two elements must be the enrichment_weight and mass!!
     */
    struct {

      /*! Normalisation factor used for the enrichment */
      float enrichment_weight;

      /*! Mass released */
      float mass;

      /*! Total metal mass released */
      float total_metal_mass;

      /*! Total mass released by each element */
      float metal_mass[chemistry_element_count];

      /*! Total mass released due to SNIa */
      float mass_from_SNIa;

      /*! Total metal mass released due to SNIa */
      float metal_mass_from_SNIa;

      /*! Total iron mass released due to SNIa */
      float Fe_mass_from_SNIa;

      /*! Total mass released due to SNII */
      float mass_from_SNII;

      /*! Total metal mass released due to SNII */
      float metal_mass_from_SNII;

      /*! Total mass released due to AGB */
      float mass_from_AGB;

      /*! Total metal mass released due to AGB */
      float metal_mass_from_AGB;

      /*! Energy change due to thermal and kinetic energy of ejecta */
      float energy;

      /*! Probability to heating neighbouring gas particle for SNII feedback */
      float SNII_heating_probability;

      /*! Change in energy from SNII feedback energy injection */
      float SNII_delta_u;

      /*! Probability to heating neighbouring gas particle for SNIa feedback */
      float SNIa_heating_probability;

      /*! Change in energy from SNIa feedback energy injection */
      float SNIa_delta_u;

      /*! Age of star particle in Myr when SNII goes off */
      float SNII_star_age_Myr;

      /*! HII region timer in SU (time since BB for cosmo runs)*/
      float HIIregion_endtime;

      /*! ID of star particle responsible for HII region */
      long long HIIregion_starid;

      /*! HII region probability */
      float HIIregion_probability;

      /*! momentum available at the given timestep */
      float momentum;

      /*! momentum weigth */
      float momentum_weight;

      /*! momentum probability */
      float momentum_probability;

      /*! momentum probability */
      float momentum_delta_v;

    } to_distribute;
  };
};

#endif /* SWIFT_FEEDBACK_STRUCT_COLIBRE_H */
