#ifndef SWIFT_DUST_IO_NONE_H
#define SWIFT_DUST_IO_NONE_H

#include "dust.h"
#include "io_properties.h"

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * No output to add
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int dust_write_particles(const struct part* parts,
                                            struct io_props* list,
                                            const int with_cosmology) {
  return 0;
}

/**
 * @brief Writes the current model of SPH to the file
 *
 * Nothing here.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The #engine.
 */
INLINE static void dust_write_flavour(hid_t h_grp, hid_t h_grp_columns,
				      const struct engine* e) {}

/**
 * @brief Create and write array for mapping elements to dust grains
 *
 * Nothing here.
 *
 * @param h_grp The HDF5 group in which to write
 * @param e The #engine.
 */
INLINE static void dust_write_composition(hid_t h_grp,
					  const struct engine* e) {}

#endif /* SWIFT_DUST_NONE_PROPERTIES_H */
