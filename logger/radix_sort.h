/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef LOGGER_LOGGER_RADIX_SORT_H
#define LOGGER_LOGGER_RADIX_SORT_H

#include "logger_index.h"

/* Number of bits to consider at a time.
 * If you wish to change this, you will need to modify some other piece of codes
 */
#define RADIX_NUMBER_BITS 4

/**
 * @brief Get the index of the bucket for the radix-counting sort.
 *
 * @param data The #index_data to consider,
 * @param i The index of the radix loop.
 */
__attribute__((always_inline)) INLINE static int radix_sort_get_bucket(
    const struct index_data *data, int i) {
  /* Consider 4 bits at a time */
  const int n = RADIX_NUMBER_BITS * i;
  const long long mask =
      (1LL << n) + (1LL << (n + 1)) + (1LL << (n + 2)) + (1LL << (n + 3));

  /* Now keep only the requested bits and
   * translate the bits to the least significant part. */
  return (data->id & mask) >> n;
}

void counting_sort(struct index_data *data, struct index_data *output, size_t N,
                   int i);
void radix_sort(struct index_data *data, size_t N);

#endif  // LOGGER_LOGGER_RADIX_SORT_H
