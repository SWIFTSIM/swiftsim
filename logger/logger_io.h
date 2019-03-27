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
/**
 * @file logger_io.h
 * @brief This file contains basic IO function.
 */
#ifndef __LOGGER_LOGGER_IO_H__
#define __LOGGER_LOGGER_IO_H__

#include "logger_header.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>

size_t io_get_file_size(int fd);
void *io_mmap_file(char *filename, size_t *file_size);
void io_munmap_file(void *map, size_t file_size);

/**
 * @brief read a mask with its offset.
 *
 * @param h #header file structure.
 * @param data Pointer to the data to read.
 * @param offset position in the file.
 * @param mask mask read from the data.
 * @param diff_offset offset difference to previous/next corresponding record read.
 * from the data.
 *
 * @return offset after the record header.
 */
__attribute__((always_inline)) INLINE static size_t io_read_mask(
    const struct header *h, void *data, size_t offset, size_t *mask,
    size_t *diff_offset) {
  /* read mask */
  if (mask) {
    *mask = 0;
    memcpy(mask, data + offset, LOGGER_MASK_SIZE);
  }
  offset += LOGGER_MASK_SIZE;

  /* read offset */
  if (diff_offset) {
    *diff_offset = 0;
    memcpy(diff_offset, data + offset, LOGGER_OFFSET_SIZE);
  }
  offset += LOGGER_OFFSET_SIZE;

  return offset;
}

/**
 * @brief read a single value in a file.
 *
 * @param data Pointer to the data to read.
 * @param size size of the data to read.
 * @param p pointer where to store the data.
 * @param offset position to read.

 * @return offset after the data written.
 */
__attribute__((always_inline)) INLINE static size_t io_read_data(
    void *data, const size_t size, void *p, size_t offset) {
  memcpy(p, data + offset, size);
  return offset + size;
};

/**
 * @brief write a single value in a file.
 *
 * @param data Pointer to the data to read.
 * @param size size of the data to write.
 * @param p pointer to the data.
 * @param offset position to write.
 *
 * @return offset after the data written.
 */
__attribute__((always_inline)) INLINE static size_t io_write_data(
    void *data, const size_t size, const void *p, size_t offset) {
  memcpy(data + offset, p, size);

  return offset + size;
};

#endif  // __LOGGER_LOGGER_IO_H__
