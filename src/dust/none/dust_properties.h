#ifndef SWIFT_DUST_NONE_PROPERTIES_H
#define SWIFT_DUST_NONE_PROPERTIES_H

#include "chemistry_struct.h"
#include "feedback_properties.h"
#include "cooling_struct.h"
#include "units.h"
#include "physical_constants.h"

/**
 * @brief Properties of the dust evolution model.
 *
 * Nothing here.
 */
struct dustevo_props {
  float *logfD;
};

/**
 * @brief initialise structure housing global dust parametrisation.
 * In particular, flags and values set in the parameter file, 
 * and any hard-coded properties.
 *
 * Nothing here.
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
}

#endif /* SWIFT_DUST_NONE_PROPERTIES_H */
