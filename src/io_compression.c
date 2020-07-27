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
    "off", "on", "D-scale-1", "D-scale-3", "D-scale-6"};

enum lossy_compression_schemes compression_scheme_from_name(const char* name) {

  for (int i = 0; i < compression_level_count; ++i) {
    if (strcmp(name, lossy_compression_schemes_names[i]) == 0)
      return (enum lossy_compression_schemes)i;
  }

  error("Invalid lossy compression scheme name: '%s'", name);
  return 0;
}

void set_hdf5_lossy_compression(hid_t h_prop, hid_t h_type,
                                const enum lossy_compression_schemes comp,
                                const char* field_name) {

  if (comp == compression_write_dscale_6) {

    /* Scale filter with a scaling by 10^6 */

    hid_t h_err = H5Pset_scaleoffset(h_prop, H5Z_SO_FLOAT_DSCALE, 6);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_dscale_3) {

    /* Scale filter with a scaling by 10^3 */

    hid_t h_err = H5Pset_scaleoffset(h_prop, H5Z_SO_FLOAT_DSCALE, 3);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else if (comp == compression_write_dscale_1) {

    /* Scale filter with a scaling by 10^1 */

    hid_t h_err = H5Pset_scaleoffset(h_prop, H5Z_SO_FLOAT_DSCALE, 1);
    if (h_err < 0)
      error("Error while setting scale-offset filter for field '%s'.",
            field_name);
  }

  else {
  }

#if 0
  
  if (strcmp(field_name, "Temperatures") == 0) {

    /* Numbers > 1 with 10-bits mantissa and 6-bits exponent */

    const int size = 4;
    const int m_size = 10;
    const int e_size = 6;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = 0;  // (1 << (e_size - 1)) - 1;

    h_type = H5Tcopy(H5T_IEEE_F32BE);
    hid_t h_err = H5Tset_fields(h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  if (strcmp(field_name, "Densities") == 0) {

    /* Numbers with 10-bits mantissa and 7-bits exponent */

    const int size = 4;
    const int m_size = 10;
    const int e_size = 7;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size - 1)) - 1;

    h_type = H5Tcopy(H5T_IEEE_F32BE);
    hid_t h_err = H5Tset_fields(h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }

  if (strcmp(field_name, "MetalMassFractions") == 0) {

    /* Numbers with 10-bits mantissa and 6-bits exponent biased such
       that the maximal valid number is 2. */

    const int size = 4;
    const int m_size = 10;
    const int e_size = 6;
    const int offset = 0;
    const int precision = m_size + e_size + 1;
    const int e_pos = offset + m_size;
    const int s_pos = e_pos + e_size;
    const int m_pos = offset;
    const int bias = (1 << (e_size)) - 2;

    h_type = H5Tcopy(H5T_IEEE_F32BE);
    hid_t h_err = H5Tset_fields(h_type, s_pos, e_pos, e_size, m_pos, m_size);
    if (h_err < 0)
      error("Error while setting type properties for field '%s'.", field_name);

    h_err = H5Tset_offset(h_type, offset);
    if (h_err < 0)
      error("Error while setting type offset properties for field '%s'.",
            field_name);

    h_err = H5Tset_precision(h_type, precision);
    if (h_err < 0)
      error("Error while setting type precision properties for field '%s'.",
            field_name);

    h_err = H5Tset_size(h_type, size);
    if (h_err < 0)
      error("Error while setting type size properties for field '%s'.",
            field_name);

    h_err = H5Tset_ebias(h_type, bias);
    if (h_err < 0)
      error("Error while setting type bias properties for field '%s'.",
            field_name);

    h_err = H5Pset_nbit(h_prop);
    if (h_err < 0)
      error("Error while setting n-bit filter for field '%s'.", field_name);
  }
#endif
}
