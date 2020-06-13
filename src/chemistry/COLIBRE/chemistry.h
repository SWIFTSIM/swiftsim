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
      "Neon",     "Magnesium", "Silicon", "Iron",     "Europium"};

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

  struct chemistry_part_data* cpd = &p->chemistry_data;

  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    cpd->diffusion_rate[elem] = 0.0f;
  }

  for (int k = 0; k < 3; k++) {
    cpd->shear_tensor[0][k] = 0.0f;
    cpd->shear_tensor[1][k] = 0.0f;
    cpd->shear_tensor[2][k] = 0.0f;
  }

  /* Initializing dmetal arrays for calculation in force loop */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    cpd->dmetal_mass_fraction[elem] = 0.0f;
  }

  /* Also have total metallicity ready for loop */
  cpd->dmetal_mass_fraction_total = 0.0f;

  /* and the metals from the different channels */
  cpd->dmetal_mass_fraction_from_SNIa = 0.0f;
  cpd->dmetal_mass_fraction_from_AGB = 0.0f;
  cpd->dmetal_mass_fraction_from_SNII = 0.0f;
  cpd->diron_mass_fraction_from_SNIa = 0.0f;
  cpd->dEu_mass_fraction_from_NSM = 0.0f;
  cpd->dEu_mass_fraction_from_CEJSN = 0.0f;
  cpd->dEu_mass_fraction_from_collapsar = 0.0f;
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

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
  const float rho = hydro_get_comoving_density(p);
  const float rho_inv = 1.0f / rho; /* 1 / rho */

  /* Convert Velocity shear tensor to physical coordinates */
  const float common_factor = h_inv_dim_plus_one * rho_inv * cosmo->a2_inv;
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      p->chemistry_data.shear_tensor[k][j] *= common_factor;
    }
  }

  /* Now add Hubble flow to diagonal terms */
  p->chemistry_data.shear_tensor[0][0] += cosmo->H;
  p->chemistry_data.shear_tensor[1][1] += cosmo->H;
  p->chemistry_data.shear_tensor[2][2] += cosmo->H;

  /* Calculate the trace of the shear tensor divided by 3 */
  float TShearTensorN = p->chemistry_data.shear_tensor[0][0] +
                        p->chemistry_data.shear_tensor[1][1] +
                        p->chemistry_data.shear_tensor[2][2];
  TShearTensorN *= (1.0f / 3.0f);

  /* Define a new shear_tensor "shear_tensor_S" */
  float shear_tensor_S[3][3];
  for (int k = 0; k < 3; k++) {
    shear_tensor_S[k][0] = 0.5 * (p->chemistry_data.shear_tensor[k][0] +
                                  p->chemistry_data.shear_tensor[0][k]);
    shear_tensor_S[k][1] = 0.5 * (p->chemistry_data.shear_tensor[k][1] +
                                  p->chemistry_data.shear_tensor[1][k]);
    shear_tensor_S[k][2] = 0.5 * (p->chemistry_data.shear_tensor[k][2] +
                                  p->chemistry_data.shear_tensor[2][k]);
  }
  shear_tensor_S[0][0] -= TShearTensorN;
  shear_tensor_S[1][1] -= TShearTensorN;
  shear_tensor_S[2][2] -= TShearTensorN;

  /* Norm of shear tensor S */
  float NormTensor2 = 0.0f;
  for (int k = 0; k < 3; k++) {
    NormTensor2 += shear_tensor_S[k][0] * shear_tensor_S[k][0] +
                   shear_tensor_S[k][1] * shear_tensor_S[k][1] +
                   shear_tensor_S[k][2] * shear_tensor_S[k][2];
  }
  const float NormTensor = sqrtf(NormTensor2);

  /* We can now combine to get the diffusion coefficient (in physical
   * coordinates)
   * This is rho a^-3 * Norm tensor (physical already) * a^2 * h^2 */
  p->chemistry_data.diffusion_coefficient =
      cd->diffusion_constant * rho * NormTensor * h * h * cosmo->a_inv;
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * Set all the fields such that no diffusion occurs.
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

  /* Getting ready for diffusion rate calculation in the force loop */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    p->chemistry_data.diffusion_rate[elem] = 0.0f;
    p->chemistry_data.dmetal_mass_fraction[elem] = 0.0f;
  }

  /* Also have total metallicity ready for loop */
  p->chemistry_data.dmetal_mass_fraction_total = 0.0f;

  /* and the metals from the different channels */
  p->chemistry_data.dmetal_mass_fraction_from_SNIa = 0.0f;
  p->chemistry_data.dmetal_mass_fraction_from_AGB = 0.0f;
  p->chemistry_data.dmetal_mass_fraction_from_SNII = 0.0f;
  p->chemistry_data.diron_mass_fraction_from_SNIa = 0.0f;
  p->chemistry_data.dEu_mass_fraction_from_NSM = 0.0f;
  p->chemistry_data.dEu_mass_fraction_from_CEJSN = 0.0f;
  p->chemistry_data.dEu_mass_fraction_from_collapsar = 0.0f;

  for (int k = 0; k < 3; k++) {
    p->chemistry_data.shear_tensor[0][k] = 0.0f;
    p->chemistry_data.shear_tensor[1][k] = 0.0f;
    p->chemistry_data.shear_tensor[2][k] = 0.0f;
  }
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

    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      p->chemistry_data.metal_mass_fraction[elem] =
          data->initial_metal_mass_fraction[elem];
    }
  }

  chemistry_init_part(p, data);

  /* Setting diffusion coefficient to zero initial value */
  p->chemistry_data.diffusion_coefficient = 0.0f;

  /* Initial value of the trackers */
  p->chemistry_data.metal_mass_tracker = 0.0f;
  p->chemistry_data.iron_mass_tracker = 0.0f;

  /* Dummy initial values to weighted redshits */
  p->chemistry_data.metal_weighted_redshift = -1.f;
  p->chemistry_data.iron_weighted_redshift = -1.f;

  p->chemistry_data.mass_from_NSM = 0.0f;
  p->chemistry_data.mass_from_CEJSN = 0.0f;
  p->chemistry_data.mass_from_collapsar = 0.0f;
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

    /* Let's check that things make sense (broadly) */

    /* H + He + Z should be ~1 */
    float total_frac = data->initial_metal_mass_fraction[chemistry_element_H] +
                       data->initial_metal_mass_fraction[chemistry_element_He] +
                       data->initial_metal_mass_fraction_total;

    if (total_frac < 0.98 || total_frac > 1.02)
      error("The abundances provided seem odd! H + He + Z = %f =/= 1.",
            total_frac);

    /* Sum of metal elements should be <= Z */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      if (elem != chemistry_element_H && elem != chemistry_element_He) {
        total_frac += data->initial_metal_mass_fraction[elem];
      }
    }

    if (total_frac > 1.02 * data->initial_metal_mass_fraction_total)
      error(
          "The abundances provided seem odd! \\sum metal elements (%f) > Z "
          "(%f)",
          total_frac, data->initial_metal_mass_fraction_total);

    /* Sum of all elements should be <= 1 */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      total_frac += data->initial_metal_mass_fraction[elem];
    }

    if (total_frac > 1.02)
      error("The abundances provided seem odd! \\sum elements (%f) > 1",
            total_frac);
  }

  /* Read diffusion constant */
  data->diffusion_constant = parser_get_param_float(
      parameter_file, "COLIBREChemistry:metal_diffusion_constant");

  /* Read time-step limiter */
  data->chemistry_time_limiter = parser_get_param_float(
      parameter_file, "COLIBREChemistry:metal_diffusion_timestep_mult");
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
 * @brief Updates the metal mass fractions after diffusion at the end of the
 * force loop.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part* p, const struct cosmology* cosmo) {

  /* Diffuse each element individually */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    p->chemistry_data.metal_mass_fraction[elem] +=
        p->chemistry_data.dmetal_mass_fraction[elem];
  }

  /* Diffuse the total metal mass fraction */
  p->chemistry_data.metal_mass_fraction_total +=
      p->chemistry_data.dmetal_mass_fraction_total;

  /* And diffuse the various trackers of metal fractions from stellar FB */
  p->chemistry_data.metal_mass_fraction_from_SNIa +=
      p->chemistry_data.dmetal_mass_fraction_from_SNIa;

  p->chemistry_data.metal_mass_fraction_from_AGB +=
      p->chemistry_data.dmetal_mass_fraction_from_AGB;

  p->chemistry_data.metal_mass_fraction_from_SNII +=
      p->chemistry_data.dmetal_mass_fraction_from_SNII;

  p->chemistry_data.iron_mass_fraction_from_SNIa +=
      p->chemistry_data.diron_mass_fraction_from_SNIa;

  const double current_mass = hydro_get_mass(p);

  /* Update metal mass from SNIa for consistency */
  p->chemistry_data.mass_from_SNIa =
      p->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;

  /* Update metal mass from SNII  */
  p->chemistry_data.mass_from_SNII =
      p->chemistry_data.metal_mass_fraction_from_SNII * current_mass;

  /* Update metal mass from AGB  */
  p->chemistry_data.mass_from_AGB =
      p->chemistry_data.metal_mass_fraction_from_AGB * current_mass;

  /* Update europium mass for channels NSM, CEJSN and collapsars  */
  const double current_Eu_mass =
      current_mass *
      p->chemistry_data.metal_mass_fraction[chemistry_element_Eu];

  p->chemistry_data.mass_from_NSM +=
      p->chemistry_data.dEu_mass_fraction_from_NSM * current_Eu_mass;

  p->chemistry_data.mass_from_CEJSN +=
      p->chemistry_data.dEu_mass_fraction_from_CEJSN * current_Eu_mass;

  p->chemistry_data.mass_from_collapsar +=
      p->chemistry_data.dEu_mass_fraction_from_collapsar * current_Eu_mass;

  /* Make sure the total metallicity is >= 0 */
  p->chemistry_data.metal_mass_fraction_total =
      max(p->chemistry_data.metal_mass_fraction_total, 0.f);
}

/**
 * @brief Computes the metal diffution time-step.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float chemistry_timestep(
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* cd, const struct part* restrict p) {

  /* h and rho in physical units */
  const float h_phys = cosmo->a * p->h;
  const float rho_phys = p->rho * cosmo->a3_inv;

  /*Diff. coeff. in physical units */
  const float coeff = p->chemistry_data.diffusion_coefficient;

  const float dt_diff = cd->chemistry_time_limiter * h_phys * h_phys *
                        rho_phys / (coeff + FLT_MIN);

  /* Convert back to co-moving coordinates */
  return dt_diff * cosmo->a2_inv;
}

/**
 * @brief Initialise the chemistry properties of a black hole with
 * the chemistry properties of the gas it is born from.
 *
 * Black holes don't store fractions so we need to use element masses.
 *
 * @param bp_data The black hole data to initialise.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_bpart_from_part(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {

  bp_data->metal_mass_total = p_data->metal_mass_fraction_total * gas_mass;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] = p_data->metal_mass_fraction[i] * gas_mass;
  }
  bp_data->mass_from_SNIa = p_data->mass_from_SNIa;
  bp_data->mass_from_SNII = p_data->mass_from_SNII;
  bp_data->mass_from_AGB = p_data->mass_from_AGB;
  bp_data->metal_mass_from_SNIa =
      p_data->metal_mass_fraction_from_SNIa * gas_mass;
  bp_data->metal_mass_from_SNII =
      p_data->metal_mass_fraction_from_SNII * gas_mass;
  bp_data->metal_mass_from_AGB =
      p_data->metal_mass_fraction_from_AGB * gas_mass;
  bp_data->iron_mass_from_SNIa =
      p_data->iron_mass_fraction_from_SNIa * gas_mass;
  bp_data->mass_from_NSM = p_data->mass_from_NSM;
  bp_data->mass_from_CEJSN = p_data->mass_from_CEJSN;
  bp_data->mass_from_collapsar = p_data->mass_from_collapsar;
}

/**
 * @brief Add the chemistry data of a gas particle to a black hole.
 *
 * Black holes don't store fractions so we need to add element masses.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_part_to_bpart(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {

  bp_data->metal_mass_total += p_data->metal_mass_fraction_total * gas_mass;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] += p_data->metal_mass_fraction[i] * gas_mass;
  }

  bp_data->mass_from_SNIa += p_data->mass_from_SNIa;
  bp_data->mass_from_SNII += p_data->mass_from_SNII;
  bp_data->mass_from_AGB += p_data->mass_from_AGB;

  bp_data->metal_mass_from_SNIa +=
      p_data->metal_mass_fraction_from_SNIa * gas_mass;
  bp_data->metal_mass_from_SNII +=
      p_data->metal_mass_fraction_from_SNII * gas_mass;
  bp_data->metal_mass_from_AGB +=
      p_data->metal_mass_fraction_from_AGB * gas_mass;
  bp_data->iron_mass_from_SNIa +=
      p_data->iron_mass_fraction_from_SNIa * gas_mass;

  bp_data->mass_from_NSM += p_data->mass_from_NSM;
  bp_data->mass_from_CEJSN += p_data->mass_from_CEJSN;
  bp_data->mass_from_collapsar += p_data->mass_from_collapsar;
}

/**
 * @brief Transfer chemistry data of a gas particle to a black hole.
 *
 * In contrast to `chemistry_add_part_to_bpart`, only a fraction of the
 * masses stored in the gas particle are transferred here. Absolute masses
 * of the gas particle are adjusted as well.
 * Black holes don't store fractions so we need to add element masses.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param nibble_mass The mass to be removed from the gas particle.
 * @param nibble_fraction The fraction of the (original) mass of the gas
 *        particle that is removed.
 */
__attribute__((always_inline)) INLINE static void
chemistry_transfer_part_to_bpart(
    struct chemistry_bpart_data* bp_data,
    struct chemistry_part_data* p_data, const double nibble_mass,
    const double nibble_fraction) {

  bp_data->metal_mass_total += p_data->metal_mass_fraction_total * nibble_mass;
  for (int i = 0; i < chemistry_element_count; ++i)
    bp_data->metal_mass[i] += p_data->metal_mass_fraction[i] * nibble_mass;

  bp_data->mass_from_SNIa += p_data->mass_from_SNIa * nibble_fraction;
  bp_data->mass_from_SNII += p_data->mass_from_SNII * nibble_fraction;
  bp_data->mass_from_AGB += p_data->mass_from_AGB * nibble_fraction;

  /* Absolute masses, so need to reduce the gas particle */
  p_data->mass_from_SNIa -= p_data->mass_from_SNIa * nibble_fraction;
  p_data->mass_from_SNII -= p_data->mass_from_SNII * nibble_fraction;
  p_data->mass_from_AGB -= p_data->mass_from_AGB * nibble_fraction;

  bp_data->metal_mass_from_SNIa +=
      p_data->metal_mass_fraction_from_SNIa * nibble_mass;
  bp_data->metal_mass_from_SNII +=
      p_data->metal_mass_fraction_from_SNII * nibble_mass;
  bp_data->metal_mass_from_AGB +=
      p_data->metal_mass_fraction_from_AGB * nibble_mass;
  bp_data->iron_mass_from_SNIa +=
      p_data->iron_mass_fraction_from_SNIa * nibble_mass;

  bp_data->mass_from_NSM += p_data->mass_from_NSM * nibble_fraction;
  bp_data->mass_from_CEJSN += p_data->mass_from_CEJSN * nibble_fraction;
  bp_data->mass_from_collapsar += p_data->mass_from_collapsar * nibble_fraction;

  /* Again, these are absolute masses so we need to adjust the gas as well */
  p_data->mass_from_NSM -= p_data->mass_from_NSM * nibble_fraction;
  p_data->mass_from_CEJSN -= p_data->mass_from_CEJSN * nibble_fraction;
  p_data->mass_from_collapsar -= p_data->mass_from_collapsar * nibble_fraction;

  /* TODO: Need to think about, and possibly implement, other transfers
   * involving e.g. diffusion or mass-weighted redshift of Fe enrichment */
}

/**
 * @brief Add the chemistry data of a black hole to another one.
 *
 * @param bp_data The black hole data to add to.
 * @param swallowed_data The black hole data to use.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_bpart_to_bpart(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_bpart_data* swallowed_data) {

  bp_data->metal_mass_total += swallowed_data->metal_mass_total;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] += swallowed_data->metal_mass[i];
  }
  bp_data->mass_from_SNIa += swallowed_data->mass_from_SNIa;
  bp_data->mass_from_SNII += swallowed_data->mass_from_SNII;
  bp_data->mass_from_AGB += swallowed_data->mass_from_AGB;
  bp_data->metal_mass_from_SNIa += swallowed_data->metal_mass_from_SNIa;
  bp_data->metal_mass_from_SNII += swallowed_data->metal_mass_from_SNII;
  bp_data->metal_mass_from_AGB += swallowed_data->metal_mass_from_AGB;
  bp_data->iron_mass_from_SNIa += swallowed_data->iron_mass_from_SNIa;
  bp_data->mass_from_NSM += swallowed_data->mass_from_NSM;
  bp_data->mass_from_CEJSN += swallowed_data->mass_from_CEJSN;
  bp_data->mass_from_collapsar += swallowed_data->mass_from_collapsar;
}

/**
 * @brief Split the metal content of a particle into n pieces
 *
 * We only need to split the fields that are not fractions.
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void chemistry_split_part(
    struct part* p, const double n) {
  p->chemistry_data.mass_from_SNIa /= n;
  p->chemistry_data.mass_from_SNII /= n;
  p->chemistry_data.mass_from_AGB /= n;
  p->chemistry_data.mass_from_NSM /= n;
  p->chemistry_data.mass_from_CEJSN /= n;
  p->chemistry_data.mass_from_collapsar /= n;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in cooling related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_cooling(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in cooling related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_cooling(const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in star formation related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in star formation related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

#endif /* SWIFT_CHEMISTRY_COLIBRE_H */
