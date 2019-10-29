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

  /*! Are we doing grain destruction by sputtering? */
  int with_sputtering;

  /*! Are we doing grain destrubtion by SNII? */
  int with_SNII_destruction;

  /*! Are we doing grain growth by accretion? */
  int with_accretion;

  /*! Are we actually cooling? */
  int with_cooling_on;

  /* ------------ Fixed parameters ------------ */
  
  /* Number of SNII per gram of star formation (assumes Chabrier 2000 IMF) */
  float specific_numSNII_cgs;
};

void dustevo_sputter_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm, const double time);

void dustevo_accretion_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm, const double time);

void dustevo_props_init(struct dustevo_props *dustevo_properties,
			struct swift_params *params,
			const struct phys_const *phys_const);

void dustevo_struct_restore(const struct dustevo_props* dustevo, FILE* stream);

void dustevo_struct_dump(const struct dustevo_props* dustevo, FILE* stream);

#endif /* SWIFT_DUSTEVO_H */
