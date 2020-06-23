#include "cooling.h"

/**
 * @brief initialise structure housing global dust parametrisation.
 * In particular, flags and values set in the parameter file, 
 * any hard-coded propertues and dust yields
 *
 * @param dp Global dust parameters for initialisation
 * @param params The parsed parameter file.
 * @param phys_const The physical constants in internal units.
 * @param us The current internal system of units.
 */
void dustevo_props_init(struct dustevo_props *dp,
			struct swift_params *params,
			struct feedback_properties *fp,
			struct cooling_function_data *cooling,
			const struct phys_const *phys_const,
			const struct unit_system *us) {

  dustevo_props_init_backend(dp, params, fp, cooling, phys_const, us);
}

