/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

#include "logger_chemistry.h"

/* Define the size of all the fields. */
#define member_size(type, member) sizeof(((type *)0)->member)

/* Hydro part */
const int chemistry_logger_field_size_part[chemistry_logger_field_part_count] = {
    member_size(struct part, chemistry_data.smoothed_metal_mass_fraction) +
    member_size(struct part, chemistry_data.metal_mass),
};

int chemistry_logger_local_to_global_part[chemistry_logger_field_part_count];

/* Stellar part */
const int chemistry_logger_field_size_spart[chemistry_logger_field_spart_count] = {
    member_size(struct spart, chemistry_data.metal_mass_fraction),
};

int chemistry_logger_local_to_global_spart[chemistry_logger_field_spart_count];
