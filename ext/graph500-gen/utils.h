/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef UTILS_H
#define UTILS_H

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#include <stddef.h>
#include <mpi.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "splittable_mrg.h"
#include "utils.h"
#include "graph_generator.h"
#include "permutation_gen.h"




inline void* xmalloc(size_t n) {
  void* p = malloc(n);
  if (!p) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", n);
    abort();
  }
  return p;
}

inline void* xcalloc(size_t n, size_t k) {
  void* p = calloc(n, k);
  if (!p) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", n);
    abort();
  }
  return p;
}

/* Get a number in [0, n) in an unbiased way. */
#ifdef __MTA__
#pragma mta inline
#endif
inline uint_fast64_t random_up_to(mrg_state* st, uint_fast64_t n) {
  /* PRNG returns values in [0, 0x7FFFFFFF) */
  /* Two iters returns values in [0, 0x3FFFFFFF00000001) */
  assert (n > 0 && n <= UINT64_C(0x3FFFFFFF00000001));
  if (n == 1) {
    return 0;
  } else if (n <= UINT64_C(0x7FFFFFFF)) {
    uint_fast64_t acc_value_limit = (UINT64_C(0x7FFFFFFF) / n) * n; /* Round down to multiple of n */
    while (1) {
      uint_fast64_t acc = mrg_get_uint_orig(st);
      if (acc >= acc_value_limit) continue;
      return acc % n;
    }
  } else if (n <= UINT64_C(0x3FFFFFFF00000001)) {
    uint_fast64_t acc_value_limit = (UINT64_C(0x3FFFFFFF00000001) / n) * n; /* Round down to multiple of n */
    while (1) {
      uint_fast64_t acc = mrg_get_uint_orig(st) * UINT64_C(0x7FFFFFFF);
      acc += mrg_get_uint_orig(st); /* Do this separately to get fixed ordering. */
      if (acc >= acc_value_limit) continue;
      return acc % n;
    }
  } else {
    /* Should have been caught before */
    return 0;
  }
}

/* Spread the two 64-bit numbers into five nonzero values in the correct
 * range. */
inline void make_mrg_seed(uint64_t userseed1, uint64_t userseed2, uint_fast32_t* seed) {
  seed[0] = (userseed1 & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed1 >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed2 & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed2 >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed2 >> 60) << 4) + (userseed1 >> 60) + 1;
}



/* Compare-and-swap; return 1 if successful or 0 otherwise. */
#ifdef __MTA__
#pragma mta inline
static inline int int64_t_cas(volatile int64_t* p, int64_t oldval, int64_t newval) {
  int64_t val = readfe(p);
  if (val == oldval) {
    writeef(p, newval);
    return 1;
  } else {
    writeef(p, val);
    return 0;
  }
}
#elif defined(GRAPH_GENERATOR_MPI) || defined(GRAPH_GENERATOR_SEQ)
/* Sequential */
#ifdef _MSC_VER
static _inline int int64_t_cas(int64_t* p, int64_t oldval, int64_t newval){
#else
static inline int int64_t_cas(int64_t* p, int64_t oldval, int64_t newval) {
#endif
  if (*p == oldval) {
    *p = newval;
    return 1;
  } else {
    return 0;
  }
}
#elif defined(GRAPH_GENERATOR_OMP)
#ifndef _MSC_VER
/* GCC intrinsic */
static inline int int64_t_cas(volatile int64_t* p, int64_t oldval, int64_t newval) {
  return __sync_bool_compare_and_swap(p, oldval, newval);
}
#else
static _inline int int64_t_cas(volatile int64_t* p, int64_t oldval, int64_t newval) {
  if (*p == oldval)
    {
        *p = newval;
        return 1; 
    } else
    { 
        return 0; 
    } 
}
#endif
#else
#ifndef _MSC_VER
#error "Need to define int64_t_cas() for your system"
#else
static _inline int int64_t_cas(volatile int64_t* p, int64_t oldval, int64_t newval) {
  if (*p == oldval)
    {
        *p = newval;
        return 1; 
    } else
    { 
        return 0; 
    } 
}
#endif
#endif


#endif /* UTILS_H */
