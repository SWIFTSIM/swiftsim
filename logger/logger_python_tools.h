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
#ifndef SWIFT_LOGGER_PYTHON_TOOLS_H
#define SWIFT_LOGGER_PYTHON_TOOLS_H

#include "../config.h"

/* Include the tools */
#include "logger_tools.h"

#if HAVE_PYTHON
/* Include numpy */
#include <numpy/arrayobject.h>

/* Structure that allows the user to easily define the
   fields in a numpy array. */
struct logger_python_field {

  /* The name of the field */
  char name[STRING_SIZE];

  /* The offset of the field in the structure */
  size_t offset;

  /* The type of data (following numpy definition) */
  int type;

  /* The number of dimension for a field */
  int dimension;
};

#define logger_loader_python_field(name, parts, field, dim, type) \
  logger_loader_python_field_function(                            \
      name, ((char *)(&parts[0].field) - (char *)parts), dim, type)

/**
 * @brief Generate a #logger_python_field structure.
 *
 * @param name The name of the field.
 * @param offset The offset of the field in the corresponding structure.
 * @param dimension The number of dimension for the field.
 * @param type The numpy data type (e.g. NPY_FLOAT32, NPY_DOUBLE, NPY_LONGLONG,
 * ...)
 *
 * @return The initialized structure.
 */
__attribute__((always_inline)) INLINE static struct logger_python_field
logger_loader_python_field_function(char *name, size_t offset, int dimension,
                                    int type) {
  struct logger_python_field ret;

  strcpy(ret.name, name);
  ret.offset = offset;
  ret.type = type;
  ret.dimension = dimension;

  return ret;
}

/**
 * @brief Add a field in a numpy descriptor.
 *
 * @param fields The dictionary containing the fields.
 * @param name The name of the field to add.
 * @param offset The offset of the field in the data structure.
 * @param type The data type of the field.
 */
__attribute__((always_inline)) INLINE static void create_field(
    PyObject *fields, char *name, size_t offset, int type) {

  PyObject *tuple = PyTuple_New(2);
  PyTuple_SetItem(tuple, 0, (PyObject *)PyArray_DescrFromType(type));
  PyTuple_SetItem(tuple, 1, PyLong_FromSize_t(offset));
  PyDict_SetItem(fields, PyUnicode_FromString(name), tuple);

}

/**
 * @brief Add a N-dimensional field in a numpy descriptor.
 *
 * @param fields The dictionary containing the fields.
 * @param name The name of the field to add.
 * @param offset The offset of the field in the data structure.
 * @param type The data type of the field.
 * @param dim The number of dimension for the field.
 */
__attribute__((always_inline)) INLINE static void create_field_ndim(
    PyObject *fields, char *name, size_t offset, int type, int dim) {

  /* Create the N-D descriptor */
    PyArray_Descr *vec = PyArray_DescrNewFromType(type);
    vec->subarray = malloc(sizeof(PyArray_ArrayDescr));
    vec->subarray->base = PyArray_DescrFromType(type);
    vec->subarray->shape = PyTuple_New(1);
    PyTuple_SetItem(vec->subarray->shape, 0, PyLong_FromSize_t(dim));

    /* Create the field */
    PyObject *tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple, 0, (PyObject *)vec);
    PyTuple_SetItem(tuple, 1, PyLong_FromSize_t(offset));
    PyDict_SetItem(fields, PyUnicode_FromString(name), tuple);
}

#endif  // HAVE_PYTHON

#endif  // SWIFT_LOGGER_PYTHON_TOOLS_H
