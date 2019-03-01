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
#ifndef SWIFT_SCALAR_DIFFUSION_STRUCT_SCALAR_H
#define SWIFT_SCALAR_DIFFUSION_STRUCT_SCALAR_H

/**
 * @brief Diffusion particle data traced by the #part.
 */
struct diffusion_part_data {
    
    /*! Scalar S at time step *before* diffusion with neighbouring particles */
    float scalar;
    
    /*! Scalar S at time step *after* diffusion with neighbouring particles  */
    float a_scalar;

    /*! Tensor of the velocity shear */
    float shear_tensor[3][3];
    
    /*! Diffusion coefficient as defined by the Smagorinsky model */
    float diffusion_coefficient;
    
    /*! Diffusion rate */
    float diffusion_rate;
};

#endif /* SWIFT_SCALAR_DIFFUSION_STRUCT_H */

