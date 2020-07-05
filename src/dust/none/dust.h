#ifndef SWIFT_DUST_NONE_H
#define SWIFT_DUST_NONE_H

/* #include "dust_struct.h" */
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
static INLINE void redistribute_dust_masses(const struct part* p, 
					    struct dustevo_props *dp) {
}

#endif /* SWIFT_DUST_NONE_H */
