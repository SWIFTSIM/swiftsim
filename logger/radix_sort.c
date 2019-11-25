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

/* Include corresponding header */
#include "radix_sort.h"

/* Include local headers */
#include "logger_index.h"

/**
 * @brief Sort the data with the radix + counting sort algorithm.
 */
void radix_sort(struct index_data *data, size_t N) {

  /* Allocate temporary memory */
  struct index_data *sorted =
      (struct index_data *)malloc(N * sizeof(struct index_data));

  if (sorted == NULL) {
    error("Failed to allocate temporary array for radix sort");
  }

  /* Create an additional pointer to the data */
  struct index_data *data_counting = data;

  /* Loop over all the bits */
  const int n_loop = 8 * sizeof(long long) / RADIX_NUMBER_BITS;
  for (int i = 0; i < n_loop; i++) {
    counting_sort(data_counting, sorted, N, i);

    /* Swap the pointers */
    struct index_data *tmp = sorted;

    sorted = data_counting;
    data_counting = tmp;
  }

  /* Copy the data back to the correct array */
  if (data_counting != data) {
    memcpy(data, data_counting, N * sizeof(struct index_data));
    sorted = data_counting;
  }

  /* Free the allocated memory */
  free(sorted);
}

/**
 * @brief Do a counting sort on only a few bits (subsort routine of radix).
 *
 * @param data The #index_data to sort.
 * @param output The sorted #index_data (already allocated).
 * @param N The number of element in data.
 * @param i The index of the radix loop.
 */
void counting_sort(struct index_data *data, struct index_data *output, size_t N,
                   int i) {
  const int n_buckets = 1 << RADIX_NUMBER_BITS;

  /* Initialize the bucket counter */
  long long counter[n_buckets];
  bzero(counter, n_buckets * sizeof(long long));

  /* Count the number of element in each bucket */
  for (size_t j = 0; j < N; j++) {
    const int k = radix_sort_get_bucket(&data[j], i);
    counter[k]++;
  }

  /* Accumulate the counter */
  for (int j = 1; j < n_buckets; j++) {
    counter[j] += counter[j - 1];
  }

  /* Now sort */
  for (size_t j1 = N; j1 > 0; j1--) {
    const size_t j = j1 - 1;

    /* Get position in the sorted array */
    const int bucket = radix_sort_get_bucket(&data[j], i);
    const size_t new_pos = counter[bucket] - 1;

    /* Set the element */
    output[new_pos] = data[j];

    /* Update the counters */
    counter[bucket]--;
  }
}
