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
#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_particle.h"
#include "logger_python_tools.h"
#include "logger_reader.h"
#include "logger_time.h"

#ifdef HAVE_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Read the minimal and maximal time.
 *
 * <b>basename</b> Base name of the logger files.
 *
 * <b>verbose</b> Verbose level.
 *
 * <b>returns</b> tuple containing min and max time.
 */
static PyObject *getTimeLimits(__attribute__((unused)) PyObject *self,
                               PyObject *args) {

  /* declare variables. */
  char *basename = NULL;

  int verbose = 0;

  /* parse arguments. */
  if (!PyArg_ParseTuple(args, "s|i", &basename, &verbose)) return NULL;

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, basename, verbose);

  if (verbose > 1) message("Reading time limits.");

  /* Get the time limits */
  double time_min = logger_reader_get_time_begin(&reader);
  double time_max = logger_reader_get_time_end(&reader);

  /* Free the memory. */
  logger_reader_free(&reader);

  /* Create the output */
  PyObject *out = PyTuple_New(2);
  PyTuple_SetItem(out, 0, PyFloat_FromDouble(time_min));
  PyTuple_SetItem(out, 1, PyFloat_FromDouble(time_max));

  return (PyObject *)out;
}

/**
 * @brief Reverse offset in log file
 *
 * <b>filename</b> string filename of the log file
 * <b>verbose</b> Verbose level
 */
static PyObject *pyReverseOffset(__attribute__((unused)) PyObject *self,
                                 PyObject *args) {
  /* input variables. */
  char *filename = NULL;

  int verbose = 0;

  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "s|i", &filename, &verbose)) return NULL;

  /* initialize the reader which reverse the offset if necessary. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);

  /* Free the reader. */
  logger_reader_free(&reader);

  return Py_BuildValue("");
}


__attribute__((always_inline)) INLINE static PyObject
*logger_loader_create_output(void **output, const int *field_indices, const int n_fields,
                             uint64_t n_tot) {

  struct logger_python_field python_fields[100];
  gravity_logger_generate_python(python_fields);

  PyObject *list = PyList_New(n_fields);
  struct logger_python_field *tmp;

  for(int i = 0; i < n_fields; i++) {
    tmp = NULL;
    for(int j = 0; j < gravity_logger_field_count; j++) {
      if (field_indices[i] == gravity_logger_mask_id[j]) {
        tmp = &python_fields[j];
        break;
      }
    }
    if (tmp == NULL) {
      error("Failed to find the required field");
    }

    PyObject *tmp_array = NULL;
    if (tmp->dimension > 1) {
      npy_intp dims[2] = {n_tot, tmp->dimension};
      tmp_array = PyArray_SimpleNewFromData(2, dims, tmp->typenum, output[i]);
    }
    else {
      npy_intp dims = n_tot;
      tmp_array = PyArray_SimpleNewFromData(1, &dims, NPY_DOUBLE, output[i]);
    }

    PyList_SetItem(list, i, tmp_array);
  }

  return list;
}


/**
 */
static PyObject *pyGetParticleData(__attribute__((unused)) PyObject *self,
                                   PyObject *args) {
  message("WARNING NOT DEALING WITH SPECIAL FLAGS AND TIMESTEP CORRECTLY");
  /* input variables. */
  char *filename = NULL;

  int verbose = 1;
  PyObject *fields = NULL;
  double time = 0;

  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "sOd|i", &filename, &fields, &time, &verbose))
    return NULL;

  if (!PyList_Check(fields)) {
    error("Expecting a list of fields");
  }

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);
  const struct header *h = &reader.log.header;

  /* Get the fields in the logger format. */
  const int n_fields = PyList_Size(fields);
  int *field_indices = (int *) malloc(n_fields * sizeof(int));
  for(int i = 0; i < n_fields; i++) {
    field_indices[i] = -1;

    PyObject *field = PyList_GetItem(fields, i);
    if (!PyUnicode_Check(field)) {
      error("Expecting list of string for the fields");
    }

    Py_ssize_t size = 0;
    const char *field_name = PyUnicode_AsUTF8AndSize(field, &size);

    for(int j = 0; j < h->masks_count; j++) {
      if (strcmp(field_name, h->masks[j].name) == 0) {
        field_indices[i] = j;
        break;
      }
    }
    if (field_indices[i] == -1) {
      error("Failed to find the field %s", field_name);
    }
  }

  /* Set the time. */
  logger_reader_set_time(&reader, time);

  /* Get the number of particles. */
  int n_type = 0;
  const uint64_t *n_part = logger_reader_get_number_particles(&reader, &n_type);
  uint64_t n_tot = 0;
  for(int i = 0; i < n_type; i++) {
    n_tot += n_part[i];
  }

  /* Allocate the output memory */
  void **output = malloc(n_fields * sizeof(void));
  for(int i = 0; i < n_fields; i++) {
    output[i] = malloc(n_tot * h->masks[field_indices[i]].size);
  }

  /* Read the particles. */
  logger_reader_read_all_particles(&reader, time, logger_reader_const,
                                   field_indices, n_fields, output,
                                   n_part);

  PyObject *array = logger_loader_create_output(output, field_indices, n_fields, n_tot);

  /* Free the reader. */
  logger_reader_free(&reader);
  free(field_indices);

  return array;
}

/* definition of the method table. */

static PyMethodDef libloggerMethods[] = {
    {"reverseOffset", pyReverseOffset, METH_VARARGS,
     "Reverse the offset (from pointing backward to forward).\n\n"
     "Parameters\n"
     "----------\n\n"
     "filename: str\n"
     "  The filename of the log file.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n"},
    {"getTimeLimits", getTimeLimits, METH_VARARGS,
     "Read the time limits of the simulation.\n\n"
     "Parameters\n"
     "----------\n\n"
     "basename: str\n"
     "  The basename of the index files.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n\n"
     "Returns\n"
     "-------\n\n"
     "times: tuple\n"
     "  time min, time max\n"},
    {"get_particle_data", pyGetParticleData, METH_VARARGS,
     "TODO"},

    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef libloggermodule = {
    PyModuleDef_HEAD_INIT,
    "liblogger",
    "Module reading a SWIFTsim logger snapshot",
    -1,
    libloggerMethods,
    NULL, /* m_slots */
    NULL, /* m_traverse */
    NULL, /* m_clear */
    NULL  /* m_free */
};

PyMODINIT_FUNC PyInit_liblogger(void) {

  /* Create the module. */
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  /* Deal with SWIFT clock */
  clocks_set_cpufreq(0);

  /* Import numpy. */
  import_array();

  return m;
}

#endif  // HAVE_PYTHON
