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
#ifndef SWIFT_HDF5_FUNCTIONS_GEAR_H
#define SWIFT_HDF5_FUNCTIONS_GEAR_H

/* Local includes. */
#include "chemistry.h"
#include "inline.h"

__attribute__((always_inline)) INLINE static void h5_open_group(
    struct swift_params* params, char *group_name, hid_t *file_id, hid_t *group_id) {

  /* Get filename. */
  char filename[256];
  parser_get_param_string(params, "GEARFeedback:YieldsTable",
			  filename);

  /* Open file. */
  *file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (*file_id < 0) error("unable to open file %s.\n", filename);

  /* Open group. */
  *group_id = H5Gopen(*file_id, group_name, H5P_DEFAULT);
  if (*group_id < 0) error("unable to open group %s.\n", group_name);
  
}


__attribute__((always_inline)) INLINE static void h5_close_group(
    hid_t file_id, hid_t group_id) {

  /* Close group */
  hid_t status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing file.");
}
#endif // SWIFT_HDF5_FUNCTIONS_GEAR_H
