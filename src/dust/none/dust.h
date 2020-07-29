#ifndef SWIFT_DUST_NONE_H
#define SWIFT_DUST_NONE_H

/* #include "dust_struct.h" */

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "parser.h"
#include "physical_constants.h"
#include "restart.h"
#include "space.h"
#include "units.h"
#include "dust_properties.h"

/**
 * @brief redistribute any dust mass back to element abundances
 * on star particle formation according to dust composition, to 
 * represent astration 
 *
 * Nothing here.
 *
 * @param p The gas particles.  
 * @param dp Global dust parameters for initialisation.
 */
static INLINE void redistribute_dust_masses(struct spart* sp,
					    const struct part* p, 
					    const struct dustevo_props *dp) {}

static INLINE void dustevo_print_backend(const struct dustevo_props *dp) {
  message("Running without dust evolution modelling.");
}

static INLINE void evolve_dust_part(const struct phys_const *phys_const,
				    const struct unit_system *us,
				    const struct cosmology *cosmo,
				    const struct hydro_props *hydro_properties,
				    const struct entropy_floor_properties *floor_props,
				    const struct cooling_function_data *cooling,
				    const struct dustevo_props *dp,
				    struct part *restrict p, struct xpart *restrict xp,
				    const float dt, const float dt_therm, const double time) {}


/**
 * @brief Prepares a particle for the smooth dust calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth dust tasks
 *
 * Nothing here.
 *
 * @param p The particle to act upon
 * @param dp #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void dust_init_part(
    struct part* restrict p, const struct dustevo_props* dp) {
}

__attribute__((always_inline)) INLINE static void dust_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct dustevo_props* dp,
    struct part* restrict p,
    struct xpart* restrict xp) {
}



#endif /* SWIFT_DUST_NONE_H */
