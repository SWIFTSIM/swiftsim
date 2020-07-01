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

typedef struct {
  PyObject_HEAD struct logger_particle part;
} PyLoggerParticle;

typedef struct {
  PyObject_HEAD struct logger_gparticle gpart;
} PyLoggerGParticle;

typedef struct {
  PyObject_HEAD struct logger_sparticle spart;
} PyLoggerSParticle;

static PyTypeObject PyLoggerParticle_Type = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "logger.Particle",
    .tp_basicsize = sizeof(PyLoggerParticle),
    .tp_doc =
        "This is the type for the hydro particles."
        "It does nothing but cleanup the interface.",
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyType_GenericNew,
};

static PyTypeObject PyLoggerGParticle_Type = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "logger.GParticle",
    .tp_basicsize = sizeof(PyLoggerGParticle),
    .tp_doc =
        "This is the type for the gravity particles."
        "It does nothing but cleanup the interface.",
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyType_GenericNew,
};

static PyTypeObject PyLoggerSParticle_Type = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "logger.SParticle",
    .tp_basicsize = sizeof(PyLoggerSParticle),
    .tp_doc =
        "This is the type for the stars particles."
        "It does nothing but cleanup the interface.",
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = PyType_GenericNew,
};

PyArray_Descr *logger_particle_descr;
PyArray_Descr *logger_gparticle_descr;
PyArray_Descr *logger_sparticle_descr;

/**
 * @brief Convert the dictionary received in input into a logger array.
 *
 * @param dict The dictionary received
 *
 * @return The array.
 */
__attribute__((always_inline)) INLINE static struct logger_particle_array
convert_and_check_input(PyObject *dict) {
  /* Check that we received a dictionnary */
  if (!PyDict_Check(dict)) {
    error("Expecting a dictionnary");
  }

  PyArrayObject *parts = (PyArrayObject *)PyDict_GetItemString(dict, "gas");
  PyArrayObject *gparts =
      (PyArrayObject *)PyDict_GetItemString(dict, "dark_matter");
  PyArrayObject *sparts = (PyArrayObject *)PyDict_GetItemString(dict, "stars");

  struct logger_particle_array out;
  logger_particle_array_init(&out);

  /* Deal with the hydro particles. */
  if (parts != NULL) {
    /* Check parts */
    if (!PyArray_Check(parts)) {
      error("Expecting a numpy array for particles.");
    }

    if (PyArray_NDIM(parts) != 1) {
      error("Expecting a 1D array of particles.");
    }

    if (PyArray_TYPE(parts) != logger_particle_descr->type_num) {
      error("Expecting an array of particles.");
    }

    /* Set the hydro */
    out.hydro.parts = PyArray_DATA(parts);
    out.hydro.n = PyArray_DIM(parts, 0);
    out.hydro.allocated_size = out.hydro.n;
  }

  /* Deal with the gravity particles. */
  if (gparts != NULL) {
    /* Check gparts */
    if (!PyArray_Check(gparts)) {
      error("Expecting a numpy array for gparticles.");
    }

    if (PyArray_NDIM(gparts) != 1) {
      error("Expecting a 1D array of gparticles.");
    }

    if (PyArray_TYPE(gparts) != logger_gparticle_descr->type_num) {
      error("Expecting an array of gparticles.");
    }

    /* Set the grav */
    out.grav.parts = PyArray_DATA(gparts);
    out.grav.n = PyArray_DIM(gparts, 0);
    out.grav.allocated_size = out.grav.n;
  }

  /* Deal with the stellar particles. */
  if (sparts != NULL) {
    /* Check parts */
    if (!PyArray_Check(sparts)) {
      error("Expecting a numpy array for sparticles.");
    }

    if (PyArray_NDIM(sparts) != 1) {
      error("Expecting a 1D array of sparticles.");
    }

    if (PyArray_TYPE(sparts) != logger_sparticle_descr->type_num) {
      error("Expecting an array of sparticles.");
    }

    /* Set the stars */
    out.stars.parts = PyArray_DATA(sparts);
    out.stars.n = PyArray_DIM(sparts, 0);
    out.stars.allocated_size = out.stars.n;
  }

  return out;
}

/** @brief Convert a #logger_particle_array into a dictionary. */
__attribute__((always_inline)) INLINE static PyObject *create_output(
    struct logger_particle_array *array) {
  /* Create the output */
  PyObject *output = PyDict_New();

  /* Check if the array has the correct size */
  if (array->hydro.n != array->hydro.allocated_size ||
      array->grav.n != array->grav.allocated_size ||
      array->stars.n != array->stars.allocated_size) {
    error("Cannot return an array that contains more memory than particles.");
  }

  /* Create the entries in the dictionary. */
  if (array->hydro.n != 0) {
    /* Create the python object */
    npy_intp n = array->hydro.n;
    PyObject *parts = (PyObject *)PyArray_SimpleNewFromData(
        /* nd */ 1, &n, logger_particle_descr->type_num, array->hydro.parts);

    /* Add it to the dict */
    PyDict_SetItem(output, PyUnicode_FromString("gas"), parts);
  }
  if (array->grav.n != 0) {
    /* Create the python object */
    npy_intp n = array->grav.n;
    PyObject *gparts = (PyObject *)PyArray_SimpleNewFromData(
        /* nd */ 1, &n, logger_gparticle_descr->type_num, array->grav.parts);

    /* Add it to the dict */
    PyDict_SetItem(output, PyUnicode_FromString("dark_matter"), gparts);
  }
  if (array->stars.n != 0) {
    /* Create the python object */
    npy_intp n = array->stars.n;
    PyObject *sparts = (PyObject *)PyArray_SimpleNewFromData(
        /* nd */ 1, &n, logger_sparticle_descr->type_num, array->stars.parts);

    /* Add it to the dict */
    PyDict_SetItem(output, PyUnicode_FromString("stars"), sparts);
  }

  return output;
}

/**
 * @brief load data from the index files.
 *
 * <b>basename</b> Base name of the logger files.
 *
 * <b>time</b> The time requested.
 *
 * <b>verbose</b> Verbose level.
 *
 * <b>returns</b> dictionnary containing the data read.
 */
static PyObject *loadSnapshotAtTime(__attribute__((unused)) PyObject *self,
                                    PyObject *args) {

  /* declare variables. */
  char *basename = NULL;

  double time = 0;
  int verbose = 0;

  /* parse arguments. */
  if (!PyArg_ParseTuple(args, "sd|i", &basename, &time, &verbose)) return NULL;

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, basename, verbose);

  if (verbose > 1) message("Reading particles.");

  /* Set the reading time */
  logger_reader_set_time(&reader, time);

  /* Get the number of particles */
  int n_type = 0;
  const uint64_t *n_parts =
      logger_reader_get_number_particles(&reader, &n_type);

  /* Get the number of particles */
  if (verbose > 0) {
    message("Found [%lu, %lu, %lu, %lu, %lu %lu] particles", n_parts[0],
            n_parts[1], n_parts[2], n_parts[3], n_parts[4], n_parts[5]);
  }

  /* Allocate the output memory */
  struct logger_particle_array array;
  logger_particle_array_allocate(&array, n_parts[swift_type_gas],
                                 n_parts[swift_type_dark_matter],
                                 n_parts[swift_type_stars], /* empty */ 0);

  /* Allows to use threads */
  Py_BEGIN_ALLOW_THREADS;

  /* Read the particle. */
  logger_reader_read_all_particles(&reader, time, logger_reader_const, &array);

  /* No need of threads anymore */
  Py_END_ALLOW_THREADS;

  /* Free the memory. */
  logger_reader_free(&reader);

  /* Create the output */
  PyObject *output = create_output(&array);

  /* No need to free the array as we are returning it */
  return output;
}

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

/**
 * @brief Move forward in time an array of particles.
 *
 * <b>filename</b> string filename of the log file.
 * <b>parts</b> The dictionary containing the particles to evolve.
 * <b>time</b> Time requested for the particles.
 * <b>verbose</b> Verbose level
 *
 * <b>returns</b> The evolved array of particles.
 */
static PyObject *pyMoveForwardInTime(__attribute__((unused)) PyObject *self,
                                     PyObject *args) {
  /* input variables. */
  char *filename = NULL;

  int verbose = 0;
  PyObject *dict = NULL;
  double time = 0;

  /* parse the arguments. */
  if (!PyArg_ParseTuple(args, "sOd|i", &filename, &dict, &time, &verbose))
    return NULL;

  struct logger_particle_array array = convert_and_check_input(dict);

  /* Create the next array. */
  struct logger_particle_array next;
  logger_particle_array_allocate(&next, array.hydro.n, array.grav.n,
                                 array.stars.n, /* empty */ 0);

  /* initialize the reader. */
  struct logger_reader reader;
  logger_reader_init(&reader, filename, verbose);

  /* Move the particles around the requested offset */
  logger_reader_move_forward(&reader, &array, &next, time,
                             /* should_interpolate */ 1);

  /* Free the reader. */
  logger_reader_free(&reader);

  /* Create the output */
  PyObject *output = create_output(&next);

  return output;
}

/* definition of the method table. */

static PyMethodDef libloggerMethods[] = {
    {"loadSnapshotAtTime", loadSnapshotAtTime, METH_VARARGS,
     "Load a snapshot directly from the logger using the index files.\n\n"
     "Parameters\n"
     "----------\n\n"
     "basename: str\n"
     "  The basename of the index files.\n\n"
     "time: double\n"
     "  The (double) time of the snapshot.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n\n"
     "Returns\n"
     "-------\n\n"
     "snapshot: dict\n"
     "  The full output generated for the whole file.\n"},
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
    {"moveForwardInTime", pyMoveForwardInTime, METH_VARARGS,
     "Move the particles forward in time.\n\n"
     "Parameters\n"
     "----------\n\n"
     "basename: str\n"
     "  The basename of the index files.\n\n"
     "parts: dict\n"
     "  The dictionary containing all the particles.\n\n"
     "time: double\n"
     "  The requested time for the particles.\n\n"
     "verbose: int, optional\n"
     "  The verbose level of the loader.\n\n"
     "Returns\n"
     "-------\n\n"
     "parts: dict\n"
     "  The dictionary containing the interpolated particles.\n"},

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

/**
 * @brief Creates the particle types for all type of particles.
 *
 * @param m The logger module.
 *
 * @return Success code (< 0 means failed).
 */
int pylogger_particle_create_typeobject(PyObject *m) {

  /* Do the hydro. */
  Py_INCREF(&PyLoggerParticle_Type);
  if (PyModule_AddObject(m, "PyLoggerParticle",
                         (PyObject *)&PyLoggerParticle_Type) < 0) {
    Py_DECREF(&PyLoggerParticle_Type);
    return -1;
  }

  /* Do the gravity. */
  Py_INCREF(&PyLoggerGParticle_Type);
  if (PyModule_AddObject(m, "PyGLoggerParticle",
                         (PyObject *)&PyLoggerGParticle_Type) < 0) {
    Py_DECREF(&PyLoggerGParticle_Type);
    return -2;
  }

  /* Do the stars. */
  Py_INCREF(&PyLoggerSParticle_Type);
  if (PyModule_AddObject(m, "PyLoggerSParticle",
                         (PyObject *)&PyLoggerSParticle_Type) < 0) {
    Py_DECREF(&PyLoggerSParticle_Type);
    return -3;
  }

  return 0;
}

/**
 * @brief Defines the numpy descriptor for the (hydro) particle.
 */
void pylogger_particle_define_descr(void) {

  /* Get the fields */
  struct logger_python_field list[100];

  int num_fields = logger_particles_generate_python(list);

  /* Generate the different fields for the dtype. */
  PyObject *names = PyTuple_New(num_fields);
  PyObject *formats = PyTuple_New(num_fields);
  PyObject *offsets = PyTuple_New(num_fields);
  for (int i = 0; i < num_fields; i++) {
    PyTuple_SetItem(names, i, PyUnicode_FromString(list[i].name));
    PyTuple_SetItem(formats, i, PyUnicode_FromString(list[i].type));
    PyTuple_SetItem(offsets, i, PyLong_FromLong(list[i].offset));
  }

  /* Generate the dtype. */
  PyObject *fields = PyDict_New();
  PyDict_SetItemString(fields, "names", names);
  PyDict_SetItemString(fields, "formats", formats);
  PyDict_SetItemString(fields, "offsets", offsets);

  /* Convert it to a descriptor. */
  PyArray_DescrConverter(fields, &logger_particle_descr);

  /* Update the size of the object (take into account non accessible fields). */
  logger_particle_descr->elsize = sizeof(struct logger_particle);

  /* Give a type to the object. */
  logger_particle_descr->typeobj = &PyLoggerParticle_Type;

  /* Register the new data type. */
  PyArray_RegisterDataType(logger_particle_descr);
}

/**
 * @brief Defines the numpy descriptor for the (gravity) particle.
 */
void pylogger_gparticle_define_descr(void) {

  /* Get the fields */
  struct logger_python_field list[100];

  int num_fields = logger_gparticles_generate_python(list);

  /* Generate the different fields for the dtype. */
  PyObject *names = PyTuple_New(num_fields);
  PyObject *formats = PyTuple_New(num_fields);
  PyObject *offsets = PyTuple_New(num_fields);
  for (int i = 0; i < num_fields; i++) {
    PyTuple_SetItem(names, i, PyUnicode_FromString(list[i].name));
    PyTuple_SetItem(formats, i, PyUnicode_FromString(list[i].type));
    PyTuple_SetItem(offsets, i, PyLong_FromLong(list[i].offset));
  }

  /* Generate the dtype. */
  PyObject *fields = PyDict_New();
  PyDict_SetItemString(fields, "names", names);
  PyDict_SetItemString(fields, "formats", formats);
  PyDict_SetItemString(fields, "offsets", offsets);

  /* Convert it to a descriptor. */
  PyArray_DescrConverter(fields, &logger_gparticle_descr);

  /* Update the size of the object (take into account non accessible fields). */
  logger_gparticle_descr->elsize = sizeof(struct logger_gparticle);

  /* Give a type to the object. */
  logger_gparticle_descr->typeobj = &PyLoggerGParticle_Type;

  /* Register the new data type */
  PyArray_RegisterDataType(logger_gparticle_descr);
}

/**
 * @brief Defines the numpy descriptor for the (stars) particle.
 */
void pylogger_sparticle_define_descr(void) {

  /* Get the fields */
  struct logger_python_field list[100];

  int num_fields = logger_sparticles_generate_python(list);

  /* Generate the different fields for the dtype. */
  PyObject *names = PyTuple_New(num_fields);
  PyObject *formats = PyTuple_New(num_fields);
  PyObject *offsets = PyTuple_New(num_fields);
  for (int i = 0; i < num_fields; i++) {
    PyTuple_SetItem(names, i, PyUnicode_FromString(list[i].name));
    PyTuple_SetItem(formats, i, PyUnicode_FromString(list[i].type));
    PyTuple_SetItem(offsets, i, PyLong_FromLong(list[i].offset));
  }

  /* Generate the dtype. */
  PyObject *fields = PyDict_New();
  PyDict_SetItemString(fields, "names", names);
  PyDict_SetItemString(fields, "formats", formats);
  PyDict_SetItemString(fields, "offsets", offsets);

  /* Convert it to a descriptor. */
  PyArray_DescrConverter(fields, &logger_sparticle_descr);

  /* Update the size of the object (take into account non accessible fields). */
  logger_sparticle_descr->elsize = sizeof(struct logger_sparticle);

  /* Give a type to the object. */
  logger_sparticle_descr->typeobj = &PyLoggerSParticle_Type;

  /* Register the new data type */
  PyArray_RegisterDataType(logger_sparticle_descr);
}

/**
 * @brief Defines the numpy descriptor for the all particles.
 */
void pylogger_all_particle_define_descr(void) {
  pylogger_particle_define_descr();
  pylogger_gparticle_define_descr();
  pylogger_sparticle_define_descr();
}

PyMODINIT_FUNC PyInit_liblogger(void) {

  /* Finalize the initialization of the types. */
  if (PyType_Ready(&PyLoggerParticle_Type) < 0) return NULL;
  if (PyType_Ready(&PyLoggerGParticle_Type) < 0) return NULL;
  if (PyType_Ready(&PyLoggerSParticle_Type) < 0) return NULL;

  /* Create the module. */
  PyObject *m;
  m = PyModule_Create(&libloggermodule);
  if (m == NULL) return NULL;

  /* Deal with SWIFT clock */
  clocks_set_cpufreq(0);

  /* Import numpy. */
  import_array();

  /* Create the types. */
  int test = pylogger_particle_create_typeobject(m);
  if (test < 0) {
    Py_DECREF(m);
    return NULL;
  }

  /* Define the descr of the logger_particle */
  pylogger_all_particle_define_descr();

  return m;
}

#endif  // HAVE_PYTHON
