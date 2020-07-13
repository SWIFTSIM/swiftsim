#ifndef SWIFT_DUST_M16_PROPERTIES_H
#define SWIFT_DUST_M16_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Standard includes */
#include <hdf5.h>

/* Local includes. */
#include "chemistry_struct.h" 
#include "feedback_properties.h"
#include "cooling_struct.h"
#include "units.h"
#include "physical_constants.h"

/**
 * @brief Stores AGB and SNII dust yields
 */

struct dust_yield_table {

  /* Array to read dust yield tables into */
  double *yield;

  /* Array to store IMF-resampled dust yield tables */
  double *yield_IMF_resampled;

};

/**
 * @brief Properties of the dust evolution model.
 */
struct dustevo_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we doing grain destruction by sputtering? */
  int with_sputtering;

  /*! Are we doing grain destruction by SNII? */
  int with_SNII_destruction;

  /*! Are we doing grain growth by accretion? */
  int with_accretion;

  /*! Are we using the subgrid T and rho? */
  int with_subgrid_props;

  /*! Are we pairing the dust fractions to cooling? */
  int pair_to_cooling;

  /* ------------ Global parameters ------------- */

  /*! Clumping factor assumed for accretion at a given density (default 1.)*/
  float clumping_factor;

  /*! Boost (> 1.) or reduction (< 1.) factor applied to dust diffusion rates (default 1.) */
  float diffusion_rate_boost;

  /* ----------- Correcting cooling tables ---------- */

  /* array of element fractions assumed to be in the dust-phase */
  float *logfD;

  /* ------------ Dust yield tables ----------------- */

  /* Yield tables for AGB and SNII  */
  struct dust_yield_table dyield_AGB;
  struct dust_yield_table dyield_SNII;
};

/**
 * @brief Rescales COLIBRE cooling / heating rates using depletions.
 *
 * This is because values in abundances arrays become entirely
 * gas-phase when running with dust (ie. no implicit dust as in 
 * COLIBRE and COLIBRE-CHIMES cooling)
 *
 * @param dustevo_properties  
 */
static INLINE void generate_dust_yield_tables(struct dustevo_props *dp,
					      struct feedback_props *fp) {

  /**
   * First read in tables for Dell'Aglia+19 AGB dust yields. These are 
   * in the same format as the Wiersma+09 chemical yields.
   * yields are stored in dp->dyield_AGB 
   *
   * Then, remove the yielded dust mass from the corresponding chemical 
   * elements, according to each grain composition
   *
   * Finally, SNII dust yields are computed directly from chemical yields
   * following Zhukovska, Gail & Trieloff 2008. Again, corresponding mass
   * is removed from chemical yields, accounting for modification factors
   * 
   * Currently, no SNIa dust creation is assumed
   **/  
  ;
}

static INLINE void scale_out_table_depletion(struct cooling_function_data* cooling){
  /**
   * Here, iterate through the 5 axes of the cooling table, scaling out
   * depletion factors, as: 
   *
   * Theating[idx] -= log10(1-pow(10, table->log10fD[idx]))
   *
   * and
   *
   * Tcooling[idx] -= log10(1-pow(10, table->log10fD[idx]))
   *
   * Where table->log10fD are depletion factors that need to be read in 
   * from the COLIBRE tables.
   *
   * Better to house this function in the cooling/COLIBRE and cooling/CHIMES
   * code? Need to modify cooling code anyway to read table depletions.
   **/
  message("Scale out COLIBRE depletion from cooling/heating rates...");

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
static INLINE void dustevo_props_init_backend(struct dustevo_props* dp,
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
  
  message("Initialising backend...");
  scale_out_table_depletion(cooling);
  /* read some parameters */

}



#endif /* SWIFT_DUST_M16_PROPERTIES_H */
