/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Camila Correa (correa@strw.leidenuniv.nl)
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
#ifndef SWIFT_DIFFUSION_IO_H
#define SWIFT_DIFFUSION_IO_H

#include "diffusion.h"
#include "io_properties.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void diffusion_read_particles(struct part* parts,
                                        struct io_props* list) {
    
    /* List what we want to read */
    list[0] = io_make_input_field("PassiveScalar", FLOAT, 1, COMPULSORY, UNIT_CONV_NO_UNITS,
                                  parts, diffusion_data.scalar);
    return 1;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void diffusion_write_particles(const struct part* parts,
                                         struct io_props* list) {
    
    /* List what we want to write */
    list[0] = io_make_output_field("PassiveScalar", FLOAT, 1, COMPULSORY, UNIT_CONV_NO_UNITS,
                                   parts, diffusion_data.scalar);
    
    list[1] = io_make_output_field("VelocityShearTensor", FLOAT, 3, 3, UNIT_CONV_NO_UNITS, parts, diffusion_data.shear_tensor);
    
    list[2] = io_make_output_field("DiffusionCoefficient", FLOAT, 1, UNIT_CONV_NO_UNITS, parts, diffusion_data.diffusion_coefficient);
    
    list[3] = io_make_output_field("DiffusionRate", FLOAT, 1, UNIT_CONV_NO_UNITS, parts, diffusion_data.diffusion_rate);
    
    return 4;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void diffusion_write_flavour(hid_t h_grp) {
}
#endif

#endif /* SWIFT_DIFFUSION_IO_H */
