#include "dust.h"

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
void dustevo_props_init(struct dustevo_props* dp,
			struct swift_params* params,
			struct feedback_props* fp,
			struct cooling_function_data* cooling,
			const struct phys_const* phys_const,
			const struct unit_system* us) {

  dustevo_props_init_backend(dp, params, fp, cooling, phys_const, us);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * Calls cooling_print_backend for the chosen cooling function.
 *
 * @param cooling The properties of the cooling function.
 */
void dustevo_print(const struct dustevo_props* dp) {

  dustevo_print_backend(dp);
}

void dustevo_struct_restore(const struct dustevo_props* dustevo, FILE* stream) {
  restart_read_blocks((void*)dustevo, sizeof(struct dustevo_props), 1, stream,
                      NULL, "dustevo function");
}

void dustevo_struct_dump(const struct dustevo_props* dustevo,
                           FILE* stream) {
  restart_write_blocks((void*)dustevo, sizeof(struct dustevo_props),
                       1, stream, "dustevo", "dustevo function");
}
