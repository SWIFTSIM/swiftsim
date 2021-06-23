/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MPI_PART_H
#define SWIFT_MPI_PART_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "error.h"
#include "helper_macros.h"

#define populate_offsets_(x) (((char *)&temp.x - (char *)&temp) / sizeof(char)),
#define populate_lengths_(x) (sizeof(temp.x)),

/**
 * @brief Register an MPI_Datatype as an indexed subset of a given struct.
 *
 * For a structure X{ a, b, c, d}; we can create a type that picks out only
 * b and c (say) by calling the macro with (X, &MPI_type, b, c);. This will
 * work for any types for the fields a,b,c,d.
 *
 * A maximum of 10 fields in X can be picked out in this way.
 *
 * @param type The structure type from which we construct the MPI type.
 * @param MPI_type (pointer) Pointer to the MPI_Datatype to create.
 * @param ... The list of fields in the structure to register as part of the
 * subset.
 */
#define create_indexed_mpi_type(type, MPI_type, ...)                         \
  ({                                                                         \
    type temp;                                                               \
                                                                             \
    /* List of size of each field in the struct. */                          \
    const int lengths[MACRO_NUM_ARGUMENTS(__VA_ARGS__)] = {                  \
        MACRO_FOR_EACH(populate_lengths_, __VA_ARGS__)};                     \
                                                                             \
    /* List of offsets of each field from the start of the struct. */        \
    const int offsets[MACRO_NUM_ARGUMENTS(__VA_ARGS__)] = {                  \
        MACRO_FOR_EACH(populate_offsets_, __VA_ARGS__)};                     \
                                                                             \
    /* Create, register and resize the new type */                           \
    if (MPI_Type_indexed(MACRO_NUM_ARGUMENTS(__VA_ARGS__), lengths, offsets, \
                         MPI_BYTE, MPI_type) != MPI_SUCCESS) {               \
      error("Failed to create indexed MPI type for " #type);                 \
    }                                                                        \
    if (MPI_Type_create_resized(*MPI_type, 0, sizeof(type), MPI_type) !=     \
        MPI_SUCCESS) {                                                       \
      error("Failed to resize MPI type for " #type);                         \
    }                                                                        \
    if (MPI_Type_commit(MPI_type) != MPI_SUCCESS) {                          \
      error("Failed to commit indexed MPI type for " #type);                 \
    }                                                                        \
  })

#endif /* SWIFT_MPI_PART_H */
