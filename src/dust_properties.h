#ifndef SWIFT_DUST_PROPERTIES_H
#define SWIFT_DUST_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

#if defined(DUST_NONE)
#include "./dust/none/dust_properties.h"
#elif defined(DUST_T20)
#include "./dust/T20/dust_properties.h"
#elif defined(DUST_M16)
#include "./dust/M16/dust_properties.h"
#else
#error "Invalid choice of dust evolution model"
#endif


#endif /* SWIFT_DUST_PROPERTIES_H */
