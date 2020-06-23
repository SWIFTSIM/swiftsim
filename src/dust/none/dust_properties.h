#ifndef SWIFT_DUST_NONE_PROPERTIES_H
#define SWIFT_DUST_NONE_PROPERTIES_H

/**
 * @brief Properties of the dust evolution model.
 *
 * Nothing here.
 */
struct dustevo_props {};

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
void dustevo_props_init_backend(struct dustevo_props *dp,
				struct swift_params *params,
				const struct phys_const *phys_const,
				const struct unit_system *us) {};

#endif /* SWIFT_DUST_NONE_PROPERTIES_H */
