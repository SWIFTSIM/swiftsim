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
#ifndef SWIFT_DIFFUSION_SCALAR_H
#define SWIFT_DIFFUSION_SCALAR_H

/**
 * @file src/metal_diffusion/DiffusionOfScalar/diffusion.h
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "diffusion_struct.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"

/**
 * @brief Prepares a particle for the diffusion calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various diffusion tasks
 * @param p The particle to act upon
 * @param dp #diffusion_part_data containing information of diffusion arrays.
 */
__attribute__((always_inline)) INLINE static void diffusion_init_part(
    struct part* restrict p) {
    
    struct diffusion_part_data* dp = &p->diffusion_data;

    dp->diffusion_coefficient = 0.f;
    dp->diffusion_rate = 0.f;
    
    for (int k = 0; k < 3; k++) {
        dp->shear_tensor[0][k] = 0.f;
        dp->shear_tensor[1][k] = 0.f;
        dp->shear_tensor[2][k] = 0.f;
    }
}

/**
 * @brief Sets the diffusion properties of the particles to a valid starting point
 * @param data The global diffusion information.
 * @param dp Pointer to the diffusion particle data.
 */
__attribute__((always_inline)) INLINE static void diffusion_first_init_part(
       struct part* restrict p){
    
    // Add initialization of the passive scalar.
        p->diffusion_data.a_scalar = p->diffusion_data.scalar;
    
    diffusion_init_part(p);
}


/**
 * @brief Prints the properties of the diffusion model to stdout.
 * @brief The #diffusion_part_data containing information about the current
 * model. */
static INLINE void diffusion_print_backend(
                                           const struct diffusion_part_data* data) {
    
    message("Diffusion model is 'SCALAR': diffusion of given passive scalar");
}

/**
 * @brief Finishes the calculation of the velocity shear tensor,
 * and calculates the diffusion coefficient at the end of the density loop.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void diffusion_coefficient_end_density(
                                                                        struct part* restrict p, const struct cosmology *cosmo) {
    
    int k;
    /* Some smoothing length multiples. */
    const float h = p->h;
    const float h_inv = 1.0f / h;                       /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
    const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
    
    const float rho_inv = 1.f / p->rho; /* 1 / rho / h^d */
    const float a_inv2 = cosmo->a2_inv;
    
    /* Velocity shear tensor in physical coordinates */
    for (k = 0; k < 3; k++){
        p->diffusion_data.shear_tensor[k][0] *= a_inv2 * h_inv_dim_plus_one * rho_inv;
        p->diffusion_data.shear_tensor[k][1] *= a_inv2 * h_inv_dim_plus_one * rho_inv;
        p->diffusion_data.shear_tensor[k][2] *= a_inv2 * h_inv_dim_plus_one * rho_inv;
    }
    
    float shear_tensor_S[3][3];
    /* Norm of shear_tensor (in physical coordinates) */
    float TShearTensorN = p->diffusion_data.shear_tensor[0][0] + p->diffusion_data.shear_tensor[1][1] + p->diffusion_data.shear_tensor[2][2];
    TShearTensorN *= (1.f/3.f) * rho_inv;
    
    /* I define a new shear_tensor "shear_tensor_S" */
    for (k = 0; k < 3; k++){
        shear_tensor_S[k][0] = 0.5 * rho_inv * (p->diffusion_data.shear_tensor[k][0] + p->diffusion_data.shear_tensor[0][k]);
        if (k == 0) shear_tensor_S[k][0] -= TShearTensorN;
        shear_tensor_S[k][1] = 0.5 * rho_inv * (p->diffusion_data.shear_tensor[k][1] + p->diffusion_data.shear_tensor[1][k]);
        if (k == 1) shear_tensor_S[k][1] -= TShearTensorN;
        shear_tensor_S[k][2] = 0.5 * rho_inv * (p->diffusion_data.shear_tensor[k][2] + p->diffusion_data.shear_tensor[2][k]);
        if (k == 2) shear_tensor_S[k][2] -= TShearTensorN;
    }
    /* Calculate the trace */
    float NormTensor = 0.f;
    for (k = 0; k < 3; k++){
        NormTensor += fabs(shear_tensor_S[k][0] * shear_tensor_S[k][0]) + fabs(shear_tensor_S[k][1] * shear_tensor_S[k][1]) + fabs(shear_tensor_S[k][2] * shear_tensor_S[k][2]);
    }
    NormTensor = sqrt(NormTensor);
    
    // We can now combine to get the diffusion coefficient //
    p->diffusion_data.diffusion_coefficient = p->rho * NormTensor * pow_dimension(h); /* rho / h^d * Norm tensor (physical coordinates) * h^d */

}

/**
 * @brief Updates the passive scalar after diffusion at the end of the graident loop.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void diffusion_end_gradient(struct part* restrict p) {
    
    p->diffusion_data.scalar = p->diffusion_data.a_scalar;
}
#endif
