/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_GEAR_H
#define SWIFT_RT_IO_GEAR_H

#define RT_LABELS_SIZE 20

/* #include "io_properties.h" */

/**
 * @file src/rt/GEAR/rt_io.h
 * @brief Main header file for GEAR M1 Closure radiative transfer
 * scheme IO routines.
 */

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_particles(const struct part* parts,
                                    struct io_props* list) {

  /* List what we want to read */

  char fieldname[30];
  int count = 0;
  for (int phg = 0; phg < RT_NGROUPS; phg++) {
    sprintf(fieldname, "RTPhotonGroup%dEnergies", phg + 1);
    list[count++] =
        io_make_input_field(fieldname, FLOAT, 1, OPTIONAL, UNIT_CONV_NO_UNITS,
                            parts, rt_data.conserved[phg].E);
    sprintf(fieldname, "RTPhotonGroup%dFluxes", phg + 1);
    list[count++] =
        io_make_input_field(fieldname, FLOAT, 3, OPTIONAL, UNIT_CONV_NO_UNITS,
                            parts, rt_data.conserved[phg].F);
  }

  return count;
}

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_stars(const struct spart* sparts,
                                struct io_props* list) {
  return 0;
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 */
INLINE static void rt_convert_conserved_photon_energies(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[g] = part->rt_data.conserved[g].E;
  }
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 */
INLINE static void rt_convert_conserved_photon_fluxes(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  int i = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[i++] = part->rt_data.conserved[g].F[0];
    ret[i++] = part->rt_data.conserved[g].F[1];
    ret[i++] = part->rt_data.conserved[g].F[2];
  }
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "RTPhotonEnergies", FLOAT, RT_NGROUPS, UNIT_CONV_NO_UNITS, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_energies,
      "Photon Energies (all groups)");
  list[1] = io_make_output_field_convert_part(
      "RTPhotonFluxes", FLOAT, 3 * RT_NGROUPS, UNIT_CONV_NO_UNITS, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_fluxes,
      "Photon Fluxes (all groups; x, y, "
      "and z coordinates)");

  return 2;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {
  return 0;
}

/**
 * @brief Write the RT model properties to the snapshot.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The engine
 * @param rtp The #rt_props
 */
INLINE static void rt_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                    const struct engine* e,
                                    const struct rt_props* rtp) {

#if defined(HAVE_HDF5)

  /* Write scheme name */
  if (rtp->hydro_controlled_injection) {
    io_write_attribute_s(h_grp, "RT Scheme",
                         RT_IMPLEMENTATION ", hydro controlled injection");
  } else {
    io_write_attribute_s(h_grp, "RT Scheme", RT_IMPLEMENTATION);
  }

  /* Write photon goup counts */
  io_write_attribute_i(h_grp, "RT photon group number", RT_NGROUPS);

  /* Write photon goup bin edges */
  io_write_attribute(h_grp, "RT photon group edges", FLOAT, rtp->photon_groups,
                     RT_NGROUPS);

  /* If without RT, we have nothing more to do. */
  const int with_rt = e->policy & engine_policy_rt;
  if (!with_rt) return;

  /* Write photon group names now */

  /* Generate Energy Group names */
  char names_energy[RT_NGROUPS * RT_LABELS_SIZE];
  for (int g = 0; g < RT_NGROUPS; g++) {
    char newEname[RT_LABELS_SIZE];
    sprintf(newEname, "Group%dEnergy", g + 1);
    strcpy(names_energy + g * RT_LABELS_SIZE, newEname);
  }

  /* Now write them down */
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, RT_LABELS_SIZE);

  hsize_t dimsE[1] = {RT_NGROUPS};
  hid_t spaceE = H5Screate_simple(1, dimsE, NULL);
  hid_t dsetE = H5Dcreate(h_grp_columns, "RTPhotonGroupEnergies", type, spaceE,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dsetE, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names_energy);
  H5Dclose(dsetE);

  /* Generate Fluxes Group Names */
  char names_fluxes[3 * RT_NGROUPS * RT_LABELS_SIZE];
  int i = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    char newFnameX[RT_LABELS_SIZE];
    sprintf(newFnameX, "Group%dFluxx", g + 1);
    strcpy(names_fluxes + i * RT_LABELS_SIZE, newFnameX);
    i++;
    char newFnameY[RT_LABELS_SIZE];
    sprintf(newFnameY, "Group%dFluxy", g + 1);
    strcpy(names_fluxes + i * RT_LABELS_SIZE, newFnameY);
    i++;
    char newFnameZ[RT_LABELS_SIZE];
    sprintf(newFnameZ, "Group%dFluxz", g + 1);
    strcpy(names_fluxes + i * RT_LABELS_SIZE, newFnameZ);
    i++;
  }

  /* Now write them down */
  hsize_t dimsF[1] = {3 * RT_NGROUPS};
  hid_t spaceF = H5Screate_simple(1, dimsF, NULL);
  hid_t dsetF = H5Dcreate(h_grp_columns, "RTPhotonGroupFluxes", type, spaceF,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dsetF, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names_fluxes);
  H5Dclose(dsetF);

  H5Tclose(type);

#endif
}

#endif /* SWIFT_RT_IO_GEAR_H */
