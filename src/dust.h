#ifndef SWIFT_DUST_H
#define SWIFT_DUST_H

/**
 * @file src/dust.h
 * @brief Branches between the different dust evolution models.
 */

/* Config parameters. */
#include "../config.h"
#include "chemistry_struct.h"

/* Import the correct dust definition */
#if defined(DUST_NONE)
#include "./dust/none/dust.h"
#include "./dust/none/dust_struct.h"
#include "./dust/none/dust_properties.h"
#elif defined(DUST_T20)
#include "./dust/T20/dust.h"
#include "./dust/T20/dust_struct.h"
#include "./dust/T20/dust_properties.h"
#elif defined(DUST_M16)
#include "./dust/M16/dust.h"
#include "./dust/M16/dust_struct.h"
#include "./dust/M16/dust_properties.h"
#else
#error "Invalid choice of dust model."
#endif

/* Common functions */

void scale_out_table_depletion(struct cooling_function_data* cooling);

void dustevo_props_init(struct dustevo_props *dp,
			struct swift_params *params,
			struct feedback_properties *fp,
			struct cooling_function_data *cooling,
			const struct phys_const *phys_const,
			const struct unit_system *us);


#endif /* SWIFT_DUST_H */
