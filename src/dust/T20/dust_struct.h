#ifndef SWIFT_DUST_STRUCT_T20_H
#define SWIFT_DUST_STRUCT_T20_H

/**
 * @brief The individual grain species traced in the model.
 */
enum grain_species { 
  grain_species_graphite = 0,
  grain_species_silicate,
  grain_species_count
 };

/**
 * @brief dust properties traced by the #part in the T20 model.
 */
struct dust_part_data {

  /*! Fraction of the particle mass in a given species of grain */
  float grain_mass_fraction[grain_species_count];

  /*! Fraction of the particle mass in *all* metals  <!! UNNEEDED? > */
  //float grain_mass_fraction_total;
};

/**
 * @brief dust properties traced by the #bpart in the T20 model.
 */
struct dust_bpart_data {

  /*! Mass in a given species of grain */
  float grain_mass[grain_species_count];

  /*! Mass in *all* grains <!! UNNEEDED?> */
  //float grain_mass_total;
};

#endif /* SWIFT_DUST_STRUCT_T20_H */
