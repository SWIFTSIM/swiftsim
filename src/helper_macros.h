/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_HELPER_MACROS_H
#define SWIFT_HELPER_MACROS_H

/* Config parameters. */
#include "../config.h"

#define SEQUENCE_10_N_(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define SEQUENCE_REV_10_N_() 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
#define NUM_ARGUMENTS_(...) SEQUENCE_10_N_(__VA_ARGS__)

/**
 * @brief Returns the number of arguments passed to the macro.
 *
 * Maximal number returned is 10.
 */
#define MACRO_NUM_ARGUMENTS(...) \
  NUM_ARGUMENTS_(__VA_ARGS__, SEQUENCE_REV_10_N_())

#define FE_0_(f)
#define FE_1_(f, X) f(X)
#define FE_2_(f, X, ...) f(X) FE_1_(f, __VA_ARGS__)
#define FE_3_(f, X, ...) f(X) FE_2_(f, __VA_ARGS__)
#define FE_4_(f, X, ...) f(X) FE_3_(f, __VA_ARGS__)
#define FE_5_(f, X, ...) f(X) FE_4_(f, __VA_ARGS__)
#define FE_6_(f, X, ...) f(X) FE_5_(f, __VA_ARGS__)
#define FE_7_(f, X, ...) f(X) FE_6_(f, __VA_ARGS__)
#define FE_8_(f, X, ...) f(X) FE_7_(f, __VA_ARGS__)
#define FE_9_(f, X, ...) f(X) FE_8_(f, __VA_ARGS__)
#define FE_10_(f, X, ...) f(X) FE_9_(f, __VA_ARGS__)

#define GET_MACRO_(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, NAME_, ...) \
  NAME_

/**
 * @brief Apply an action to all the other arguments of the macro.
 *
 * Maximal number of other arguments is 10.
 */
#define MACRO_FOR_EACH(action, ...)                                      \
  GET_MACRO_(_0, __VA_ARGS__, FE_10_, FE_9_, FE_8_, FE_7_, FE_6_, FE_5_, \
             FE_4_, FE_3_, FE_2_, FE_1_, FE_0_)                          \
  (action, __VA_ARGS__)

#endif /* SWIFT_HELPER_MACROS_H */
