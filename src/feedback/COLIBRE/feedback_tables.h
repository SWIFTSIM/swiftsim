/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COLIBRE_FEEDBACK_TABLES_H
#define SWIFT_COLIBRE_FEEDBACK_TABLES_H

/**
 * @brief Allocates and reads early feedback tables, flattens and stores them in
 * the feedback properties.
 *
 * @param feedback_props the #feedback_props data struct to read the tables
 * into.
 */
INLINE static void read_feedback_tables(struct feedback_props *fp) {

#ifdef HAVE_HDF5

  hid_t dataset;
  herr_t status;

  hid_t tempfile_id =
      H5Fopen(fp->early_feedback_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0)
    error("unable to open file %s\n", fp->early_feedback_table_path);

  /* read sizes of array dimensions */
  dataset = H5Dopen(tempfile_id, "Header/NMETALLICITIES", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &fp->HII_nr_metbins);
  if (status < 0) error("error reading number of metallicities \n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset: number of metallicities \n");

  dataset = H5Dopen(tempfile_id, "Header/NAGES", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   &fp->HII_nr_agebins);
  if (status < 0) error("error reading number of ages \n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset: number of ages\n");

  /* allocate arrays */
  if (posix_memalign((void **)&fp->HII_log10_Zbins, SWIFT_STRUCT_ALIGNMENT,
                     fp->HII_nr_metbins * sizeof(float)) != 0)
    error("Failed to allocate metallicity bins\n");
  if (posix_memalign((void **)&fp->HII_agebins, SWIFT_STRUCT_ALIGNMENT,
                     fp->HII_nr_agebins * sizeof(float)) != 0)
    error("Failed to allocate age bins\n");
  if (posix_memalign((void **)&fp->HII_log10_Qcum, SWIFT_STRUCT_ALIGNMENT,
                     fp->HII_nr_metbins * fp->HII_nr_agebins * sizeof(float)) !=
      0)
    error("Failed to allocate Q array\n");
  if (posix_memalign((void **)&fp->SW_log10_Pcum, SWIFT_STRUCT_ALIGNMENT,
                     fp->HII_nr_metbins * fp->HII_nr_agebins * sizeof(float)) !=
      0)
    error("Failed to allocate P array\n");

  /* read in the metallicity bins */
  dataset = H5Dopen(tempfile_id, "MetallicityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   fp->HII_log10_Zbins);
  if (status < 0) error("error reading metallicity bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset: metallicity bins");

  /* read in the stellar ages bins */
  dataset = H5Dopen(tempfile_id, "AgeBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   fp->HII_agebins);
  if (status < 0) error("error reading age bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset: age bins");

  /* read in cumulative ionizing photons */
  dataset = H5Dopen(tempfile_id, "logQcumulative", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   fp->HII_log10_Qcum);
  if (status < 0) error("error reading Q\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset: logQcumulative");

  /* read in cumulative momentum input */
  dataset = H5Dopen(tempfile_id, "logPcumulative", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   fp->SW_log10_Pcum);
  if (status < 0) error("error reading P\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset: logPcumulative");

  /* Close the file */
  status = H5Fclose(tempfile_id);
  if (status < 0) error("error closing file");

#else
  error("Need HDF5 to read early feedback tables");
#endif
}

#endif /* SWIFT_COLIBRE_FEEDBACK_TABLES_H */
