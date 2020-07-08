#ifndef SWIFT_DUST_H
#define SWIFT_DUST_H

/**
 * @file src/dust.h
 * @brief Branches between the different dust evolution models.
 */

/* Config parameters. */
#include "../config.h"

/* Import the correct dust definition */
#if defined(DUST_NONE)
#include "./dust/none/dust.h"
#include "./dust/none/dust_struct.h"
#include "./dust/none/dust_properties.h"
#include "./dust/none/dust_yield_tables.h"
#elif defined(DUST_T20)
#include "./dust/T20/dust.h"
#include "./dust/T20/dust_struct.h"
#include "./dust/T20/dust_properties.h"
#include "./dust/T20/dust_yield_tables.h"
#elif defined(DUST_M16)
#include "./dust/M16/dust.h"
#include "./dust/M16/dust_struct.h"
#include "./dust/M16/dust_properties.h"
#else
#error "Invalid choice of dust model."
#endif

/* additional imports */
#include "chemistry_struct.h"
#include "feedback_properties.h"
#include "cooling_struct.h"
#include "units.h"
#include "physical_constants.h"
#include "parser.h"
#include "restart.h"
#include "space.h"



/* Common functions */

void scale_out_table_depletion(struct cooling_function_data* cooling);

void dustevo_props_init(struct dustevo_props *dp,
			struct swift_params *params,
			struct feedback_props *fp,
			struct cooling_function_data *cooling,
			const struct phys_const *phys_const,
			const struct unit_system *us);

void dustevo_struct_dump(const struct dustevo_props* dustevo,
			 FILE* stream);

void dustevo_struct_restore(const struct dustevo_props* dustevo, FILE* stream);

#endif /* SWIFT_DUST_H */
