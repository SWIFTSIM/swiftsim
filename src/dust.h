/* Currently a modified version of cooling.h */

#ifndef SWIFT_DUSTEVO_H
#define SWIFT_DUSTEVO_H

/**
 * @file src/dust.h
 * @brief dust evolution function declarations
 */

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "parser.h"
#include "physical_constants.h"
#include "restart.h"
#include "space.h"
#include "units.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct feedback_props;
struct space;

/* void dustevo_sputter_part(void); */
struct dustevo_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we actually evolving dust? This flag can be turned off to emulate the master branch */
  int with_dust_evolution;

  /*! What dust model implementation are we using? */
  int model_type;

  /*! Are we doing grain destruction by sputtering? */
  int with_sputtering;

  /*! Are we doing grain destruction by SNII? */
  int with_SNII_destruction;

  /*! Are we doing grain growth by accretion? */
  int with_accretion;

  /*! Are we using the subgrid T and rho? */
  int with_subgrid_props;

  /*! Are we actually cooling? */
  int with_cooling_on;

  /*! Are we pairing the dust fractions to cooling? */
  int pair_to_cooling;

  /* ------------ Fixed parameters ------------ */
  
  /* SNII events per mass of star formation in cgs (assumes Chabrier 2000 IMF) */
  float specific_numSNII;

  /* ------------ grain species element mass ratios ------------- */

  float comp_Gra[chemistry_element_count];
  float comp_Ide[chemistry_element_count];
  float comp_Sil[chemistry_element_count];
  float comp_Mgd[chemistry_element_count]; 
  float comp_Fed[chemistry_element_count];

  /* ------------ mapping refractory elements to grains ------------ */
  int refractory_idx[5];
  int dust_idx[5];

  /* ------------ constants that are used a lot ---------- */
  /* const float number_density_to_cgs; */ // setting here not working - see dust.c

};

void redistribute_dust_to_element_abundances(struct spart* sp,
					     const struct dustevo_props *dp);

void evolve_dust_part(const struct phys_const *phys_const,
		      const struct unit_system *us,
		      const struct cosmology *cosmo,
		      const struct hydro_props *hydro_properties,
		      const struct entropy_floor_properties *floor_props,
		      const struct cooling_function_data *cooling,
		      const struct dustevo_props *dp,
		      struct part *restrict p, struct xpart *restrict xp,
		      const float dt, const float dt_therm, const double time);

void evolve_dust_part_m16(const struct phys_const *phys_const,
		      const struct unit_system *us,
		      const struct cosmology *cosmo,
		      const struct hydro_props *hydro_properties,
		      const struct entropy_floor_properties *floor_props,
		      const struct cooling_function_data *cooling,
		      const struct dustevo_props *dp,
		      struct part *restrict p, struct xpart *restrict xp,
		      const float dt, const float dt_therm, const double time);


void dustevo_sputter_part(const struct phys_const *phys_const,
			  const struct unit_system *us,
			  const struct cosmology *cosmo,
			  const struct hydro_props *hydro_properties,
			  const struct entropy_floor_properties *floor_props,
			  const struct cooling_function_data *cooling,
			  const struct dustevo_props *dp,
			  struct part *restrict p, struct xpart *restrict xp,
			  const float dt, const float dt_therm, const double time);

void dustevo_accretion_part(const struct phys_const *phys_const,
			    const struct unit_system *us,
			    const struct cosmology *cosmo,
			    const struct hydro_props *hydro_properties,
			    const struct entropy_floor_properties *floor_props,
			    const struct cooling_function_data *cooling,
			    const struct dustevo_props *dp,
			    struct part *restrict p, struct xpart *restrict xp,
			    const float dt, const float dt_therm, const double time);

void dustevo_SNIIdestruction_part(const struct phys_const *phys_const,
				  const struct unit_system *us,
				  const struct cosmology *cosmo,
				  const struct hydro_props *hydro_properties,
				  const struct entropy_floor_properties *floor_props,
				  const struct cooling_function_data *cooling,
				  const struct dustevo_props *dp,
				  struct part *restrict p, struct xpart *restrict xp,
				  const float dt, const float dt_therm, const double time);

void dustevo_props_init(struct dustevo_props *dustevo_properties,
			struct swift_params *params,
			const struct phys_const *phys_const,
			const struct unit_system *us);

void dustevo_struct_restore(const struct dustevo_props* dustevo, FILE* stream);

void dustevo_struct_dump(const struct dustevo_props* dustevo, FILE* stream);

#endif /* SWIFT_DUSTEVO_H */
