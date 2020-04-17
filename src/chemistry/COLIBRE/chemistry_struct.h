/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenunuiv.nl)
 *               2020 Camila Correa (c.a.correa@uva.nl)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_COLIBRE_H
#define SWIFT_CHEMISTRY_STRUCT_COLIBRE_H

/**
 * @brief The individual elements traced in the COLIBRE model.
 */

enum chemistry_element {
  chemistry_element_H = 0,
  chemistry_element_He,
  chemistry_element_C,
  chemistry_element_N,
  chemistry_element_O,
  chemistry_element_Ne,
  chemistry_element_Mg,
  chemistry_element_Si,
  chemistry_element_Fe,
  chemistry_element_Eu,
  chemistry_element_count
};

/**
 * @brief Global chemical abundance information in the COLIBRE model.
 */
struct chemistry_global_data {

  /*! Fraction of the particle mass in given elements at the start of the run */
  float initial_metal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals at the start of the run */
  float initial_metal_mass_fraction_total;

  /*! Diffusion constant to calculate Smagorinsky diffusion coefficient */
  /* This constant is a free parameter we read in params file */
  float diffusion_constant;

  /*! Dimensionless pre-factor used in the diffusion time-step condition */
  float chemistry_time_limiter;
};

/**
 * @brief Chemical abundances traced by the #part in the COLIBRE model.
 */
struct chemistry_part_data {

  /*! Fraction of the particle mass in a given element */
  float metal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals */
  float metal_mass_fraction_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Fraction of total gas mass in metals coming from SNIa */
  float metal_mass_fraction_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Fraction of total gas mass in metals coming from AGB */
  float metal_mass_fraction_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Fraction of total gas mass in metals coming from SNII */
  float metal_mass_fraction_from_SNII;

  /*! Fraction of total gas mass in Iron coming from SNIa */
  float iron_mass_fraction_from_SNIa;

  /*! Metal weighted redshift (defined in Wiersma+ 2010, eq. 3) */
  float metal_weighted_redshift;

  /*! Metal weighted redshift */
  float iron_weighted_redshift;

  /*! Tracker of metal mass received by enrichment events times time */
  float metal_mass_tracker;

  /*! Tracker of iron mass received by SNIa events times time */
  float iron_mass_tracker;

  /*! Fraction of the particle mass in a given element accumulated via diffusion
   * diffusion since the last active step */
  float dmetal_mass_fraction[chemistry_element_count];

  /*! Fraction of the particle mass in metals accumulated via diffusion
   * diffusion since the last active step */
  float dmetal_mass_fraction_total;

  /*! Fraction of the particle mass in metals from SNIa accumulated via
   * diffusion diffusion since the last active step */
  float dmetal_mass_fraction_from_SNIa;

  /*! Fraction of the particle mass in metals from AGB accumulated via diffusion
   * diffusion since the last active step */
  float dmetal_mass_fraction_from_AGB;

  /*! Fraction of the particle mass in metals from SNII accumulated via
   * diffusion diffusion since the last active step */
  float dmetal_mass_fraction_from_SNII;

  /*! Fraction of the particle mass in iron from SNIa accumulated via diffusion
   * diffusion since the last active step */
  float diron_mass_fraction_from_SNIa;

  /*! Tensor of the velocity shear */
  /* (Calculated in physical coordinates) */
  float shear_tensor[3][3];

  /*! Diffusion coefficient as defined by the Smagorinsky model */
  /* (Calculated in physical coordinates) */
  float diffusion_coefficient;

  /*! Diffusion rate */
  /* (Calculated in physical coordinates) */
  float diffusion_rate[chemistry_element_count];
};

/**
 * @brief Chemical abundances traced by the #bpart in the EAGLE model.
 */
struct chemistry_bpart_data {

  /*! Mass in a given element */
  float metal_mass[chemistry_element_count];

  /*! Mass in *all* metals */
  float metal_mass_total;

  /*! Mass coming from SNIa */
  float mass_from_SNIa;

  /*! Mass coming from AGB */
  float mass_from_AGB;

  /*! Mass coming from SNII */
  float mass_from_SNII;

  /*! Metal mass coming from SNIa */
  float metal_mass_from_SNIa;

  /*! Metal mass coming from AGB */
  float metal_mass_from_AGB;

  /*! Metal mass coming from SNII */
  float metal_mass_from_SNII;

  /*! Iron mass coming from SNIa */
  float iron_mass_from_SNIa;
};

#endif /* SWIFT_CHEMISTRY_STRUCT_COLIBRE_H */
