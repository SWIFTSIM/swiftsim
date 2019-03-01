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
    struct part* restrict p, const struct diffusion_global_data* cd) {
    
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
       const struct diffusion_global_data* data,
       struct part* restrict p){
    
    // Add initialization of the passive scalar.
        p->diffusion_data.scalar = data->initial_scalar;
    
    diffusion_init_part(p, data);
}


/**
 * @brief Initialises the passive scalar array to be diffused.
 *
 * @param parameter_file The parsed parameter file.
 * @param data The properties to initialise.
 */
static INLINE void diffusion_init_backend(struct swift_params* parameter_file,
                                          struct diffusion_global_data* data) {
    
    /* Read the passive scalar from initial conditions */
    data->initial_scalar = parser_get_opt_param_float(parameter_file, "Diffusion:init_passive_scalar", -1);
}

/**
 * @brief Prints the properties of the diffusion model to stdout.
 * @brief The #diffusion_global_data containing information about the current
 * model. */
static INLINE void diffusion_print_backend(
                                           const struct diffusion_global_data* data) {
    
    message("Diffusion model is 'SCALAR': diffusion of given passive scalar");
}

/**
 * @brief Calculates the diffusion coefficient at the end of the density loop.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void diffusion_coefficient_end_density(
                                                                        struct part* restrict p) {
    
    /* Some data from part */
    const float h = p->h;
    const float density = p->rho; /* 1 / h^d * rho */
    struct diffusion_part_data* dp = &p->diffusion_data;
    float shear_tensor_S[3][3];
    int k;
    /* Norm of shear_tensor from density loop */
    float TShearTensorN = (1.f/3.f) * (p->diffusion_data.shear_tensor[0][0] + p->diffusion_data.shear_tensor[1][1] + p->diffusion_data.shear_tensor[2][2]) / density;
    
    /* I define a new shear_tensor "shear_tensor_S" */
    for (k = 0; k < 3; k++){
        shear_tensor_S[k][0] = (0.5 / density) * (p->diffusion_data.shear_tensor[k][0] + p->diffusion_data.shear_tensor[0][k]);
        if (k == 0) shear_tensor_S[k][0] -= TShearTensorN;
        shear_tensor_S[k][1] = (0.5 / density) * (p->diffusion_data.shear_tensor[k][1] + p->diffusion_data.shear_tensor[1][k]);
        if (k == 1) shear_tensor_S[k][1] -= TShearTensorN;
        shear_tensor_S[k][2] = (0.5 / density) * (p->diffusion_data.shear_tensor[k][2] + p->diffusion_data.shear_tensor[2][k]);
        if (k == 2) shear_tensor_S[k][2] -= TShearTensorN;
    }
    /* Calculate the trace */
    float NormTensor = 0.f;
    for (k = 0; k < 3; k++){
        NormTensor += fabs(shear_tensor_S[k][0] * shear_tensor_S[k][0]) + fabs(shear_tensor_S[k][1] * shear_tensor_S[k][1]) + fabs(shear_tensor_S[k][2] * shear_tensor_S[k][2]);
    }
    NormTensor = sqrt(NormTensor);
    
    // ok, combine to get the diffusion coefficient //
    dp->diffusion_coefficient = density * NormTensor * h * h; /*Check code units*/

}
#endif
