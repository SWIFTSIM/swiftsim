#ifndef SWIFT_DUST_NONE_STRUCT_H
#define SWIFT_DUST_NONE_STRUCT_H

/**
 * @brief The individual grain species traced in the model.
 */
enum grain_species { grain_species_count = 0 };

/**
 * @brief dust properties traced by the #part. Empty.
 */
struct dust_part_data {
  /* <!! Can avoid setting somehow, as in chemistry/none?> */
  float grain_mass_fraction[grain_species_count];
};

/**
 * @brief dust properties traced by the #bpart. Empty.
 */
struct dust_bpart_data {
  float grain_mass_fraction[grain_species_count];
};

#endif /* SWIFT_DUST_NONE_STRUCT_H */
