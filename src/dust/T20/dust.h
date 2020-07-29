#ifndef SWIFT_DUST_T20_H
#define SWIFT_DUST_T20_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "parser.h"
#include "physical_constants.h"
#include "restart.h"
#include "space.h"
#include "units.h"
#include "dust_properties.h"
#include "chemistry_struct.h"

void redistribute_dust_masses(struct spart* sp,
			      const struct part* p, 
			      const struct dustevo_props *dp);

void dustevo_print_backend(const struct dustevo_props *dp);

void evolve_dust_part(const struct phys_const *phys_const,
		      const struct unit_system *us,
		      const struct cosmology *cosmo,
		      const struct hydro_props *hydro_properties,
		      const struct entropy_floor_properties *floor_props,
		      const struct cooling_function_data *cooling,
		      const struct dustevo_props *dp,
		      struct part *restrict p, struct xpart *restrict xp,
		      const float dt, const float dt_therm, const double time);


/**
 * @brief Return a string containing the name of a given #grain_species.
 */
__attribute__((always_inline)) INLINE static const char*
dust_get_grain_name(enum grain_species grain) {

  static const char* grain_species_names[grain_species_count] = {
      "Graphite", "Silicates"};

  return grain_species_names[grain];
}


/**
 * @brief Prepares a particle for the smooth dust calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth dust tasks
 *
 * @param p The particle to act upon
 * @param dp #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void dust_init_part(
    struct part* restrict p, const struct dustevo_props* dp) {

  /* DO NOTHING FOR NOW - AWAITING DIFFUSION */
  ;
}

__attribute__((always_inline)) INLINE static void dust_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct dustevo_props* dp,
    struct part* restrict p,
    struct xpart* restrict xp) {

  // Add initialization of fields in dust_part_data struct.
    for (int grain = 0; grain < grain_species_count; ++grain) {
      p->dust_data.grain_mass_fraction[grain] =
          dp->initial_grain_mass_fraction[grain];
    }
  

    dust_init_part(p, dp);
}

#endif /* SWIFT_DUST_T20_PROPERTIES_H */
