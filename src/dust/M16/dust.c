/* Local includes. */

#include "chemistry.h"
#include "dust_properties.h"
#include "dust_yield_tables.h"
#include "dust.h"
#include <string.h>
/* /\** */
/*  * @brief redistribute any dust mass back to element abundances */
/*  * on star particle formation according to dust composition, to  */
/*  * represent astration  */
/*  * */
/*  * @param p The gas particles.   */
/*  * @param dp Global dust parameters for initialisation. */
/*  *\/ */
/* void redistribute_dust_masses(const struct part* p,  */
/* 			      struct dustevo_props *dp) { */
/*   /\**  */
/*    * iterate through grain species and element types and */
/*    * redistribute dust abundances to element abundances, */
/*    * according to composition array */
/*    *\/ */
/*   ; */
/* } */

/**
 * @brief Prints the dust evolution model to stdout.
 *
 * @param dust #dustevo_props struct.
 */
void dustevo_print_backend(const struct dustevo_props *dp) {
    message("Running with an approximated McKinnon et al. (2016) dust evolution model.");
}

/**
 * @brief initialise structure housing global dust parametrisation.
 * In particular, flags and values set in the parameter file, 
 * and any hard-coded properties
 *
 * @param dp Global dust parameters for initialisation
 * @param params The parsed parameter file.
 * @param phys_const The physical constants in internal units.
 * @param us The current internal system of units.
 */
void dustevo_props_init_backend(struct dustevo_props* dp,
				struct swift_params* params,
				struct feedback_props* fp,
				struct cooling_function_data* cooling,
				const struct phys_const* phys_const,
				const struct unit_system* us) {

/**
 * initialise structure housing global dust parametrisation.
 * First, set global parameters: 
 *
 * initialise dustevo_props struct object and store:
 * - Parameter file flags (e.g. with_accretion, with_sputtering)
 * - Parameter file values (e.g diffusion boost, clumping factor)
 * - Hard coded values (e.g. grain radii, density, composition)
 *
 *
 * Then, using the yield tables from feedback_props, compute 
 * and set dust yields for each of the  dust species 
 *
 * Finally apply corections to the cooling and heating rates
 * to remove implicit dust, via the scale_out_table_depletion
 * function.
 **/  
  
  /* read some parameters */

  /* <!! SET IN PARAMETER FILE> */
  dp->initial_grain_mass_fraction[grain_species_C] = 0.;
  dp->initial_grain_mass_fraction[grain_species_O] = 0.;
  dp->initial_grain_mass_fraction[grain_species_Mg] = 0.;
  dp->initial_grain_mass_fraction[grain_species_Si] = 0.;
  dp->initial_grain_mass_fraction[grain_species_Fe] = 0.;

  /* set assumed abundance patterns to Wiersma et al (2009a) */
  memset(dp->abundance_pattern, 0, sizeof dp->abundance_pattern);
  dp->abundance_pattern[chemistry_element_H] = 7.0649785e-01 / 0.0127;
  dp->abundance_pattern[chemistry_element_He] = 2.8055534e-01  / 0.0127;
  dp->abundance_pattern[chemistry_element_C] = 2.0665436e-03  / 0.0127;
  dp->abundance_pattern[chemistry_element_N] = 8.3562563e-04  / 0.0127;
  dp->abundance_pattern[chemistry_element_O] = 5.4926244e-03  / 0.0127;
  dp->abundance_pattern[chemistry_element_Ne] = 1.4144605e-03  / 0.0127;
  dp->abundance_pattern[chemistry_element_Mg] = 5.9070642e-04 / 0.0127;
  dp->abundance_pattern[chemistry_element_Si] = 6.8258739e-04 / 0.0127;
  dp->abundance_pattern[chemistry_element_Fe] = 1.1032152e-03 / 0.0127;

  /* set element atomic weights */
  memset(dp->atomic_weight, 0, sizeof dp->atomic_weight);
  dp->atomic_weight[chemistry_element_H] = 1.0079;
  dp->atomic_weight[chemistry_element_He] = 4.0026;
  dp->atomic_weight[chemistry_element_C] = 12.0107;
  dp->atomic_weight[chemistry_element_N] = 14.0067;
  dp->atomic_weight[chemistry_element_O] = 15.9994;
  dp->atomic_weight[chemistry_element_Ne] = 20.1797;
  dp->atomic_weight[chemistry_element_Mg] = 24.305;
  dp->atomic_weight[chemistry_element_Si] = 28.0855;
  dp->atomic_weight[chemistry_element_Fe] = 55.845;

  /* set condensation fractions */
  memset(dp->condensation_frac, 0, sizeof dp->condensation_frac);
  dp->condensation_frac[grain_species_C] = 0.5;
  dp->condensation_frac[grain_species_O] = 0.;
  dp->condensation_frac[grain_species_Mg] = 0.8;
  dp->condensation_frac[grain_species_Si] = 0.8;
  dp->condensation_frac[grain_species_Fe] = 0.8;


  /* set grain composition */
  float gcomp1[1] = {1.};
  float gcomp2[1] = {1.};
  float gcomp3[1] = {1.};
  float gcomp4[1] = {1.};
  float gcomp5[1] = {1.};

  dp->grain_element_mfrac[grain_species_C] = gcomp1;
  dp->grain_element_mfrac[grain_species_O] = gcomp2;
  dp->grain_element_mfrac[grain_species_Mg] = gcomp3;
  dp->grain_element_mfrac[grain_species_Si] = gcomp4;
  dp->grain_element_mfrac[grain_species_Fe] = gcomp5;

  /* set chemistry indices composition */
  int gidx1[1] = {chemistry_element_C};
  int gidx2[1] = {chemistry_element_O};
  int gidx3[1] = {chemistry_element_Mg};
  int gidx4[1] = {chemistry_element_Si};
  int gidx5[1] = {chemistry_element_Fe};

  dp->grain_element_indices[grain_species_C] = gidx1;
  dp->grain_element_indices[grain_species_O] = gidx2;
  dp->grain_element_indices[grain_species_Mg] = gidx3;
  dp->grain_element_indices[grain_species_Si] = gidx4;
  dp->grain_element_indices[grain_species_Fe] = gidx5;
  
  /* set element count contributing to each grain */
  dp->grain_element_count[grain_species_C] = 1;
  dp->grain_element_count[grain_species_O] = 1;
  dp->grain_element_count[grain_species_Mg] = 1;
  dp->grain_element_count[grain_species_Si] = 1;
  dp->grain_element_count[grain_species_Fe] = 1;

  /** NOTE 1: total metallicity yields untouched here, so Z represents the conserved dust + gas phase metals **/

  /** NOTE 2: only the IMF resampled tables are modified in fp, while plain .yield arrays are unchanged (the 
   *  original yield tables are only used to compute these and are already modified via the SNII yield factors) **/


  initialise_dyield_tables(fp, dp);
  message("%f", dp->dyield_SNII.yield_IMF_resampled[4999]);
  compute_SNII_dyield(fp, dp);
  compute_AGB_dyield(fp, dp);
  print_dyield_tables(fp, dp);

}
