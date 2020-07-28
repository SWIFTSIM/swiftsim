/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020  Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "io_compression.h"

/* Local includes. */
#include "error.h"

/* Some standard headers. */
#include <stdlib.h>
#include <string.h>

/* Compression level names. */
const char* lossy_compression_schemes_names[compression_level_count] = {
    "off", "on", "D-scale-1", "D-scale-3", "D-scale-6", "f-mantissa-10"};

enum lossy_compression_schemes compression_scheme_from_name(const char* name) {

  for (int i = 0; i < compression_level_count; ++i) {
    if (strcmp(name, lossy_compression_schemes_names[i]) == 0)
      return (enum lossy_compression_schemes)i;
  }

  error("Invalid lossy compression scheme name: '%s'", name);
  return (enum lossy_compression_schemes)0;
}

void set_hdf5_lossy_compression(hid_t* h_prop, hid_t* h_type,
                                const enum lossy_compression_schemes comp,
                                const char* field_name) {

  if (comp == compression_do_not_write) {
    error(
        "Applying a lossy compression filter to a field labelled as 'do not "
        "write'");
  }

  else if (comp == compression_write_d_scale_6) {

    /* Scale filter with a scaling by 10^6 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 6);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_3) {

    /* Scale filter with a scaling by 10^3 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 3);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_d_scale_1) {

    /* Scale filter with a scaling by 10^1 */

    hid_t h_err = H5Pset_scaleoffset(*h_prop, H5Z_SO_FLOAT_DSCALE, 1);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_f_mantissa_10) {

    /* Float numbers with 10-bits mantissa and 8-bits exponent */

    /* Note a regular IEEE-754 float has:
     * - size = 4
     * - m_size = 23
     * - e_size = 8
     * i.e. 23 + 8 + 1 (the sign bit) == 32 bits (== 4 bytes) */

    const int size = 4;
    const int m_size = 10;
    const int e_size = 8;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    H5Tclose(*h_type);
    *h_type = H5Tcopy(H5T_IEEE_F32BE);
    hid_t h_err = H5Tset_fields(*h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(*h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(*h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(*h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(*h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(*h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }
}
