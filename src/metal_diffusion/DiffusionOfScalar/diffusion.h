/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Camila Correa (correa@strw.leidenuniv.nl)
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
#ifndef SWIFT_DIFFUSION_H
#define SWIFT_DIFFUSION_H

/**
 * @file src/metal_diffusion/DiffusionOfScalar/diffusion.h
 */


/**
 * @brief Prepares a particle for the smooth metal calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth metallicity tasks
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void diffusion_init_part(
    struct part* restrict p, const struct diffusion_part_data* dp) {
    
    struct diffusion_part_data* cpd = &p->chemistry_data;
    
    for (int i = 0; i < chemistry_element_count; i++) {
        cpd->smoothed_metal_mass_fraction[i] = 0.f;
    }
    
    cpd->smoothed_metal_mass_fraction_total = 0.f;
    cpd->smoothed_iron_mass_fraction_from_SNIa = 0.f;
}
