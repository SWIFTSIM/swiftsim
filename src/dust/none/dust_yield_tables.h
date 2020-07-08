#ifndef SWIFT_DUST_NONE_TABLES_H
#define SWIFT_DUST_NONE_TABLES_H

#include <hdf5.h>
#include "inline.h"
#include "dust.h"

static INLINE void depletion_correct_rates(struct cooling_function_data *cooling,
					   struct dustevo_props *dp) {
}

static INLINE void read_colibre_depletion(hid_t id, struct dustevo_props *dp) {
}

#endif /* SWIFT_DUST_NONE_TABLES_H */
