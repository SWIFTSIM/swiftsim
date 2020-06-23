#ifndef SWIFT_DUST_STRUCT_NONE_H
#define SWIFT_DUST_STRUCT_NONE_H

/**
 * @file src/dust/none/dust_struct.h
 * @brief Empty infrastructure for running without dust evolution
 */


/**
 * @brief Global dust modelling information.
 *
 * Nothing here.
 */
struct dustevo_props {};

/**
 * @brief The individual grain species traced in the model.
 */
enum grain_species { grain_species_count = 0 };

#endif /* SWIFT_DUST_STRUCT_NONE_H */
