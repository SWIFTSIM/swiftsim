#ifndef SWIFT_DUST_M16_PROPERTIES_H
#define SWIFT_DUST_M16_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"
#include "inline.h"

/**
 * @brief Stores AGB and SNII dust yields
 */

struct dust_yield_tables {

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

  /* ------------ Dust yield tables ----------------- */

  /* Yield tables for AGB and SNII  */
  struct dust_yield_table dyield_AGB;
  struct dust_yield_table dyield_SNII;


};

/**
 * @brief initialise structure housing global dust parametrisation.
 * In particular, flags and values set in the parameter file, 
 * and any hard-coded propertues
 *
 * @param dustevo_properties Global dust parameters for initialisation
 * @param params The parsed parameter file.
 * @param phys_const The physical constants in internal units.
 * @param us The current internal system of units.
 */
void dustevo_props_init(struct dustevo_props *dustevo_properties,
			struct swift_params *params,
			const struct phys_const *phys_const,
			const struct unit_system *us) {}

#endif /* SWIFT_DUST_M16_PROPERTIES_H */
