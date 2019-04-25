/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_CHEMISTRY_COLIBRE_H
#define SWIFT_CHEMISTRY_COLIBRE_H

/**
 * @file src/chemistry/COLIBRE/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
chemistry_get_element_name(enum chemistry_element elem) {

  static const char* chemistry_element_names[chemistry_element_count] = {
      "Hydrogen", "Helium",    "Carbon",  "Nitrogen", "Oxygen",
      "Neon",     "Magnesium", "Silicon", "Iron"};

  return chemistry_element_names[elem];
}

/**
 * @brief Prepares a particle for the smooth metal calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth metallicity tasks
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part* restrict p, const struct chemistry_global_data* cd) {

    struct diffusion_part_data* dp = &p->diffusion_data;

    dp->diffusion_coefficient = 0.0f;
    
    for (int elem = 0; elem < chemistry_element_count; ++elem){
        dp->diffusion_rate[elem] = 0.0f;
    }
    
    for (int k = 0; k < 3; k++) {
        dp->shear_tensor[0][k] = 0.0f;
        dp->shear_tensor[1][k] = 0.0f;
        dp->shear_tensor[2][k] = 0.0f;
    }
}

/**
 * @brief Finishes the calculation of the velocity shear tensor,
 * and calculates the diffusion coefficient at the end of the density loop.
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 */

__attribute__((always_inline)) INLINE static void chemistry_end_density(
     struct part* restrict p, const struct chemistry_global_data* cd,
     const struct cosmology* cosmo) {

    int k;
    /* Some smoothing length multiples. */
    const float h = p->h;
    const float h_inv = 1.0f / h;                       /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
    const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
    
    const float rho_inv = 1.0f / p->rho; /* 1 / rho */
    /*const float a_inv2 = cosmo->a2_inv;*/
    
    /* Velocity shear tensor (check physical coordinates) */
    for (k = 0; k < 3; k++){
        p->diffusion_data.shear_tensor[k][0] *= h_inv_dim_plus_one * rho_inv;
        p->diffusion_data.shear_tensor[k][1] *= h_inv_dim_plus_one * rho_inv;
        p->diffusion_data.shear_tensor[k][2] *= h_inv_dim_plus_one * rho_inv;
    }
    
    float shear_tensor_S[3][3];
    /* Norm of shear_tensor (in physical coordinates) */
    float TShearTensorN = p->diffusion_data.shear_tensor[0][0] + p->diffusion_data.shear_tensor[1][1] + p->diffusion_data.shear_tensor[2][2];
    TShearTensorN *= (1.0f/3.0f);
    
    /* I define a new shear_tensor "shear_tensor_S" */
    for (k = 0; k < 3; k++){
        shear_tensor_S[k][0] = 0.5 * (p->diffusion_data.shear_tensor[k][0] + p->diffusion_data.shear_tensor[0][k]);
        if (k == 0) shear_tensor_S[k][0] -= TShearTensorN;
        shear_tensor_S[k][1] = 0.5 * (p->diffusion_data.shear_tensor[k][1] + p->diffusion_data.shear_tensor[1][k]);
        if (k == 1) shear_tensor_S[k][1] -= TShearTensorN;
        shear_tensor_S[k][2] = 0.5 * (p->diffusion_data.shear_tensor[k][2] + p->diffusion_data.shear_tensor[2][k]);
        if (k == 2) shear_tensor_S[k][2] -= TShearTensorN;
    }
    /* Calculate the trace */
    float NormTensor = 0.0f;
    for (k = 0; k < 3; k++){
        NormTensor += shear_tensor_S[k][0] * shear_tensor_S[k][0] + shear_tensor_S[k][1] * shear_tensor_S[k][1] + shear_tensor_S[k][2] * shear_tensor_S[k][2];
    }
    NormTensor = sqrt(NormTensor);
    
    // We can now combine to get the diffusion coefficient //
    p->diffusion_data.diffusion_coefficient = p->rho * NormTensor * h * h;  /* rho * Norm tensor (physical coordinates) * h^2 */
    
    /* Velocity shear tensor goes back to zero for next loop */
    for (k = 0; k < 3; k++){
        p->diffusion_data.shear_tensor[k][0] = 0.0f;
        p->diffusion_data.shear_tensor[k][1] = 0.0f;
        p->diffusion_data.shear_tensor[k][2] = 0.0f;
    }
    
    /* Getting ready for diffusion rate calculation in <FORCE LOOP> */
    for (int elem = 0; elem < chemistry_element_count; ++elem){
        p->diffusion_data.diffusion_rate[elem] = 0.0f;
    }
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_has_no_neighbours(struct part* restrict p,
                                 struct xpart* restrict xp,
                                 const struct chemistry_global_data* cd,
                                 const struct cosmology* cosmo) {

  //struct chemistry_part_data* cpd = &p->chemistry_data;
  // CC: Not sure what to do here, better to add message
    error("Needs implementing!");
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param data The global chemistry information.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

  // Add initialization of all other fields in chemistry_part_data struct.
  if (data->initial_metal_mass_fraction_total != -1) {
    p->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;

    for (int elem = 0; elem < chemistry_element_count; ++elem){
        p->chemistry_data.metal_mass_fraction[elem] = data->initial_metal_mass_fraction[elem];
        p->diffusion_data.dmetal_mass_fraction[elem] = data->initial_metal_mass_fraction[elem];
    }
  }
    
  /* I will add a modification to the initial set-up until enrichment is ON */
  /* Let's say metal rich core of solar metallicity, and outskirts of a */
  /* factor of 1000 lower metallicity */
    float a = 1.0f/1000.0f;
    float b = (1.0f - a)/(1.0f-p->chemistry_data.metal_mass_fraction_total)+a;
    
    float r2 = (p->x[0]-2282.34) * (p->x[0]-2282.34);
    r2 += (p->x[1]-2282.34) * (p->x[1]-2282.34);
    if (r2 > 142.57){
        for (int elem = 0; elem < 2; ++elem){
            p->chemistry_data.metal_mass_fraction[elem] *= b;
            p->diffusion_data.dmetal_mass_fraction[elem] *= b;
        }
        p->chemistry_data.metal_mass_fraction_total *= a;
        
        for (int elem = 2; elem < chemistry_element_count; ++elem){
            p->chemistry_data.metal_mass_fraction[elem] *= a;
            p->diffusion_data.dmetal_mass_fraction[elem] *= a;
        }
    }
    
  chemistry_init_part(p, data);
}


/**
 * @brief Initialises the chemistry properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(struct swift_params* parameter_file,
                                          const struct unit_system* us,
                                          const struct phys_const* phys_const,
                                          struct chemistry_global_data* data) {

  /* Read the total metallicity */
  data->initial_metal_mass_fraction_total = parser_get_opt_param_float(
      parameter_file, "COLIBREChemistry:init_abundance_metal", -1);

  if (data->initial_metal_mass_fraction_total != -1) {
    /* Read the individual mass fractions */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      char buffer[50];
      sprintf(buffer, "COLIBREChemistry:init_abundance_%s",
              chemistry_get_element_name((enum chemistry_element)elem));

      data->initial_metal_mass_fraction[elem] =
          parser_get_param_float(parameter_file, buffer);
    }
  }
}

/**
 * @brief Sets the chemistry properties of the sparticles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
                                                                             const struct chemistry_global_data* data, struct spart* restrict sp) {
    
    /* Initialize mass fractions for total metals and each metal individually */
    if (data->initial_metal_mass_fraction_total != -1) {
        sp->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;
        
        for (int elem = 0; elem < chemistry_element_count; ++elem)
        sp->chemistry_data.metal_mass_fraction[elem] =
        data->initial_metal_mass_fraction[elem];
    }
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data* data) {

  message("Chemistry model is 'COLIBRE' tracking %d elements.",
          chemistry_element_count);
}

/**
 * @brief Updates the metal mass fractions after diffusion at the end of the <FORCE LOOP>
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(struct part* restrict p, const struct cosmology* cosmo) {
    
    for (int elem = 0; elem < chemistry_element_count; ++elem){
        p->chemistry_data.metal_mass_fraction[elem] = p->diffusion_data.dmetal_mass_fraction[elem];
    }
    
    struct chemistry_part_data* cpd = &p->chemistry_data;
    
    for (int i = 0; i < chemistry_element_count; i++) {
        /* Final operation on the density (add self-contribution). */
        cpd->smoothed_metal_mass_fraction[i] = cpd->metal_mass_fraction[i];
    }
    
    /* Smooth mass fraction of all metals */
    cpd->smoothed_metal_mass_fraction_total = cpd->metal_mass_fraction_total;
    
    /* Smooth iron mass fraction from SNIa */
    cpd->smoothed_iron_mass_fraction_from_SNIa = cpd->iron_mass_fraction_from_SNIa;
}


#endif /* SWIFT_CHEMISTRY_COLIBRE_H */
