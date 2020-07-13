#ifndef SWIFT_DUST_T20_H
#define SWIFT_DUST_T20_H

#include "dust_properties.h"

void redistribute_dust_masses(const struct part* p, 
			      struct dustevo_props *dp);

void dustevo_print_backend(const struct dustevo_props *dp);

#endif /* SWIFT_DUST_T20_PROPERTIES_H */
