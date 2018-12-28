/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_ATOMIC_H
#define SWIFT_ATOMIC_H

#if STDC_VERSION >= 201112L
  #if !defined(__GNUC__) ||( __GNUC__ > 4 || \
    (__GNUC__ == 4 && (__GNUC_MINOR >= 9 )))
    #define SWIFT_MODERN_ATOMICS
  #endif
#endif

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "inline.h"
#include "minmax.h"

//#if STDC_VERSION >= 201112L

#ifdef NOT_DEFINED
#include <stdatomic.h>

#define atomic_add(v, i) atomic_fetch_add(v, i)
#define atomic_sub(v, i) atomic_fetch_sub(v, i)
#define atomic_or(v, i) atomic_fetch_or(v, i)
#define atomic_inc(v) atomic_fetch_add(v, 1)
#define atomic_dec(v) atomic_fetch_sub(v, 1)
#define
#ifdef SWIFT_MODERN_ATOMICS
typedef _Atomic(float) atomic_float;
typedef _Atomic(double) atomic_double;
#define atomic_cas(obj,expected,desired) _Generic((desired), \
  int: bool_atomic_compare_and_swap_i, \
  long long: bool_atomic_compare_and_swap_ll, \
  float: bool_atomic_compare_and_swap_f, \
  default: bool_atomic_compare_and_swap_ll \
  )(obj,expected,desired)
#else
typedef float atomic_float;
typedef double atomic_double;
#define atomic_cas(obj,expected,desired) _Generic((desired), \
  int: bool_atomic_compare_and_swap_i, \
  long long: bool_atomic_compare_and_swap_ll, \
  default: val_atomic_compare_and_swap_ll \
  )(obj,expected,desired)
#endif
#define atomic_swap(v, n) atomic_exchange(v,n)

_Bool bool_atomic_compare_and_swap_i( atomic_int *obj, int *expected, int desired){
  int preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong(obj, expected, desired);
  *expected = preexpected;
  return retval;
}

_Bool bool_atomic_compare_and_swap_ll( atomic_llong *obj, long long *expected, long long desired){
  long long preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong(obj, expected, desired);
  *expected = preexpected;
  return retval;
}

#ifdef SWIFT_MODERN_ATOMICS
_Bool bool_atomic_compare_and_swap_f( atomic_float *obj, float *expected, float desired){
  float preexpected = *expected;
  _Bool retval = atomic_compare_exchange_strong( obj, expected, desired);
  *expected = preexpected;
  return retval;
}

#endif
/**
 * @brief Atomic min operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point min that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
#ifdef SWIFT_MODERN_ATOMICS
static void atomic_min_f( atomic_float *obj, float const y){

  float test_val, old_val, new_val;
  old_val = *obj;

  do{
    test_val = old_val;
    new_val = fmin(old_val, y);
  }while(!atomic_cas( obj, &test_val, new_val));

}
#else
static void atomic_min_f( volatile float *const address, float const y) {

  atomic_int *const int_ptr = (atomic_int *) address;

  typedef union{
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = *address;

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = fmin(old_val.as_float, y);
  } while (!atomic_cas(int_ptr, &test_val.as_int, new_val.as_int));

}

#endif

/**
 * @brief Atomic max operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point max that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
#ifdef SWIFT_MODERN_ATOMICS
static void atomic_max_f( atomic_float *obj, float const y){

  float test_val, old_val, new_val;
  old_val = *obj;

  do{
    test_val = old_val;
    new_val = fmax(old_val, y);
  }while(!atomic_cas(obj, &test_val, new_val));

}
#else
static void atomic_max_f( volatile float *const address, float const y) {

  atomic_int *const int_ptr = (atomic_int *) address;

  typedef union{
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = *address;

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = fmax(old_val.as_float, y);
  } while (!atomic_cas(int_ptr, &test_val.as_int, new_val.as_int));

}
#endif

#ifdef SWIFT_MODERN_ATOMICS
#define atomic_add_f(v, i) atomic_fetch_add(v, i)
#else
/**
 * @brief Atomic add operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_add_f(
    volatile float *const address, const float y) {

  atomic_int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = *address;

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = old_val.as_float + y;
  } while (!atomic_cas(int_ptr, &test_val.as_int, new_val.as_int));
}
#endif

/**
 * @brief Atomic add operation on doubles.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
#ifdef SWIFT_MODERN_ATOMICS
#define atomic_add_d(v, i) atomic_fetch_add(v, i)
#else
__attribute__((always_inline)) INLINE static void atomic_add_d(
    volatile double *const address, const double y) {

  atomic_llong *const long_long_ptr = (long long *)address;

  typedef union {
    double as_double;
    long long as_long_long;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_double = *address;

  do {
    test_val.as_long_long = old_val.as_long_long;
    new_val.as_double = old_val.as_double + y;
  } while (!atomic_cas(long_long_ptr, &test_val.as_long_long, new_val.as_long_long));
}
#endif

#else
//Old GNU99 implementation
typedef int volatile atomic_int;
typedef unsigned int volatile atomic_uint;
typedef size_t volatile atomic_size_t;
typedef short atomic_short;
typedef double atomic_double;
typedef float atomic_float;
#define atomic_add(v, i) __sync_fetch_and_add(v, i)
#define atomic_sub(v, i) __sync_fetch_and_sub(v, i)
#define atomic_or(v, i) __sync_fetch_and_or(v, i)
#define atomic_inc(v) atomic_add(v, 1)
#define atomic_dec(v) atomic_sub(v, 1)
#define atomic_cas(v, o, n) __sync_bool_compare_and_swap(v, o, n)
#define atomic_load(v) __sync_val_compare_and_swap(v,0,0)
#define atomic_swap(v, n) __sync_lock_test_and_set(v, n)

/**
 * @brief Atomic min operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point min that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_min_f(
    volatile float *const address, const float y) {

  int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = *address;

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = min(old_val.as_float, y);
  } while (!atomic_cas(int_ptr, test_val.as_int, new_val.as_int));
}

/**
 * @brief Atomic max operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point max that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_max_f(
    volatile float *const address, const float y) {

  int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = *address;

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = max(old_val.as_float, y);
  } while (!atomic_cas(int_ptr, test_val.as_int, new_val.as_int));
}

/**
 * @brief Atomic add operation on floats.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_add_f(
    volatile float *const address, const float y) {

  int *const int_ptr = (int *)address;

  typedef union {
    float as_float;
    int as_int;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_float = *address;

  do {
    test_val.as_int = old_val.as_int;
    new_val.as_float = old_val.as_float + y;
  } while ( !atomic_cas(int_ptr, test_val.as_int, new_val.as_int));
}

/**
 * @brief Atomic add operation on doubles.
 *
 * This is a text-book implementation based on an atomic CAS.
 *
 * We create a temporary union to cope with the int-only atomic CAS
 * and the floating-point add that we want.
 *
 * @param address The address to update.
 * @param y The value to update the address with.
 */
__attribute__((always_inline)) INLINE static void atomic_add_d(
    volatile double *const address, const double y) {

  long long *const long_long_ptr = (long long *)address;

  typedef union {
    double as_double;
    long long as_long_long;
  } cast_type;

  cast_type test_val, old_val, new_val;
  old_val.as_double = *address;

  do {
    test_val.as_long_long = old_val.as_long_long;
    new_val.as_double = old_val.as_double + y;
  } while ( !atomic_cas(long_long_ptr, test_val.as_long_long, new_val.as_long_long));
}
#endif

#endif /* SWIFT_ATOMIC_H */
