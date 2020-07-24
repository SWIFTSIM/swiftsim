#ifndef SWIFT_DUST_T20_H
#define SWIFT_DUST_T20_H

#include "dust_properties.h"

void redistribute_dust_masses(const struct part* p, 
			      struct dustevo_props *dp);

void dustevo_print_backend(const struct dustevo_props *dp);

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
