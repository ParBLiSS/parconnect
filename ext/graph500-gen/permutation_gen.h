/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef PERMUTATION_GEN_H
#define PERMUTATION_GEN_H

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include <stdint.h>
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "splittable_mrg.h"
#include "graph_generator.h"
#include "permutation_gen.h"
#include "utils.h"

typedef struct slot_data {
  int64_t index, value;
} slot_data;

/* This code defines a simple closed-indexing hash table.  It is used to speed
 * up the rand_sort algorithm given below.  Elements with -1 as index are
 * unused; others are used. */

static inline void hashtable_insert(slot_data* ht, int64_t ht_size, int64_t index, int64_t value, int64_t hashval) {
  int64_t i;
  for (i = hashval; i < ht_size; ++i) {
    if (int64_t_cas(&ht[i].index, (int64_t)(-1), index)) {
      ht[i].value = value;
      return;
    }
  }
  for (i = 0; i < hashval; ++i) {
    if (int64_t_cas(&ht[i].index, (int64_t)(-1), index)) {
      ht[i].value = value;
      return;
    }
  }
  assert (!"Should not happen: overflow in hash table");
}

static inline int hashtable_count_key(const slot_data* ht, int64_t ht_size, int64_t index, int64_t hashval) {
  int c = 0;
  int64_t i;
  for (i = hashval; i < ht_size && ht[i].index != (int64_t)(-1); ++i) {
    if (ht[i].index == index) ++c;
  }
  if (i == ht_size) {
    for (i = 0; i < hashval && ht[i].index != (int64_t)(-1); ++i) {
      if (ht[i].index == index) ++c;
    }
  }
  return c;
}

static inline int hashtable_get_values(const slot_data* ht, int64_t ht_size, int64_t index, int64_t hashval, int64_t* result) {
  int x = 0;
  int64_t i;
  for (i = hashval; i < ht_size && ht[i].index != (int64_t)(-1); ++i) {
    if (ht[i].index == index) {
      result[x++] = ht[i].value;
    }
  }
  if (i == ht_size) {
    for (i = 0; i < hashval && ht[i].index != (int64_t)(-1); ++i) {
      if (ht[i].index == index) {
        result[x++] = ht[i].value;
      }
    }
  }
  return x;
}

static inline void selection_sort(int64_t* a, int64_t n) {
  int64_t i, j;
  if (n <= 1) return;
  for (i = 0; i + 1 < n; ++i) {
    int64_t minpos = i;
    for (j = i + 1; j < n; ++j) {
      if (a[j] < a[minpos]) minpos = j;
    }
    if (minpos != i) {
      int64_t t = a[minpos];
      a[minpos] = a[i];
      a[i] = t;
    }
  }
}

static inline void randomly_permute(int64_t* a, int64_t n, mrg_state* st) {
  int64_t i, j;
  if (n <= 1) return;
  for (i = n - 1; i > 0; --i) {
    j = random_up_to(st, i + 1);
    if (i != j) {
      int64_t t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
  }
}

/* Exclusive prefix sum on ints; returns sum of overall input array */
static inline int int_prefix_sum(int* out, const int* in, size_t n) {
  size_t i;
  if (n == 0) return 0;
  out[0] = 0;
  for (i = 1; i < n; ++i) out[i] = out[i - 1] + in[i - 1];
  return out[n - 1] + in[n - 1];
}

/* For MPI distributed memory. */
inline void rand_sort_mpi(MPI_Comm comm, mrg_state* st, int64_t n,
                   int64_t* result_size_ptr,
                   int64_t** result_ptr /* Allocated using xmalloc() by
                   rand_sort_mpi */) {
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  /* Make MPI data type for slot_data. */
  MPI_Datatype slot_data_type;
  {
    int blocklens[] = {1, 1};
    MPI_Aint temp_base, indices[2];
    slot_data temp;
    MPI_Get_address(&temp, &temp_base);
    MPI_Get_address(&temp.index, &indices[0]);
    MPI_Get_address(&temp.value, &indices[1]);
    indices[0] -= temp_base;
    indices[1] -= temp_base;
    MPI_Datatype old_types[] = {INT64_T_MPI_TYPE, INT64_T_MPI_TYPE};
    MPI_Type_create_struct(2, blocklens, indices, old_types, &slot_data_type);
    MPI_Type_commit(&slot_data_type);
  }

  int64_t total_hash_table_size = 2 * n + 128; /* Must be >n, preferably larger for performance */

  /* Hash table is distributed by blocks: first (total_hash_table_size % size)
   * are of size (total_hash_table_size / size + 1), rest are of size
   * (total_hash_table_size / size).  This distribution is necessary so that
   * the permutation can easily be assembled at the end of the function. */
  int64_t ht_base_block_size = total_hash_table_size / size;
  int ht_block_size_cutoff_rank = total_hash_table_size % size;
  int64_t ht_block_size_cutoff_index = ht_block_size_cutoff_rank * (ht_base_block_size + 1);
  int64_t ht_my_size = ht_base_block_size + (rank < ht_block_size_cutoff_rank);
  int64_t ht_my_start = (rank < ht_block_size_cutoff_rank) ?
                           rank * (ht_base_block_size + 1) :
                           ht_block_size_cutoff_index + (rank - ht_block_size_cutoff_rank) * ht_base_block_size;
  int64_t ht_my_end = ht_my_start + ht_my_size;
#define HT_OWNER(e) \
    (((e) < ht_block_size_cutoff_index) ? \
     (e) / (ht_base_block_size + 1) : \
     ht_block_size_cutoff_rank + ((e) - ht_block_size_cutoff_index) / ht_base_block_size)
#define HT_LOCAL(e) ((e) - ht_my_start)

  /* Input elements to scramble are distributed cyclically for simplicity;
   * their distribution does not matter. */
  int64_t elt_my_size = (n / size) + (rank < n % size);

  int64_t i;

  /* Cache the key-value pairs to avoid PRNG skip operations.  Count the number
   * of pairs going to each destination processor. */
  slot_data* kv_pairs = (slot_data*)xmalloc(elt_my_size * sizeof(slot_data));
  int* outcounts = (int*)xcalloc(size, sizeof(int)); /* Relies on zero-init */
  for (i = 0; i < elt_my_size; ++i) {
    mrg_state new_st = *st;
    mrg_skip(&new_st, 1, i * size + rank, 0);
    int64_t index = (int64_t)random_up_to(&new_st, total_hash_table_size);
    int64_t owner = HT_OWNER(index);
    assert (owner >= 0 && owner < size);
    ++outcounts[owner];
    kv_pairs[i].index = index;
    kv_pairs[i].value = i * size + rank;
  }

  int* outdispls = (int*)xmalloc(size * sizeof(int));
  int total_outcount = int_prefix_sum(outdispls, outcounts, size);

  slot_data* outdata = (slot_data*)xmalloc(total_outcount * sizeof(slot_data));
  int* outoffsets = (int*)xmalloc(size * sizeof(int));
  memcpy(outoffsets, outdispls, size * sizeof(int));

  /* Put the key-value pairs into the output buffer, sorted by destination, to
   * get ready for MPI_Alltoallv. */
  for (i = 0; i < elt_my_size; ++i) {
    int64_t index = kv_pairs[i].index;
    int64_t owner = HT_OWNER(index);
    outdata[outoffsets[owner]] = kv_pairs[i];
    ++outoffsets[owner];
  }
  free(kv_pairs); kv_pairs = NULL;
  for (i = 0; i < size; ++i) {
    assert (outoffsets[i] == outdispls[i] + outcounts[i]);
  }
  free(outoffsets); outoffsets = NULL;

  int* incounts = (int*)xmalloc(size * sizeof(int));

  /* Send data counts. */
  MPI_Alltoall(outcounts, 1, MPI_INT, incounts, 1, MPI_INT, comm);

  int* indispls = (int*)xmalloc(size * sizeof(int));
  int total_incount = int_prefix_sum(indispls, incounts, size);

  slot_data* indata = (slot_data*)xmalloc(total_incount * sizeof(slot_data));

  /* Send data to put into hash table. */
  MPI_Alltoallv(outdata, outcounts, outdispls, slot_data_type,
                indata, incounts, indispls, slot_data_type,
                comm);

  free(outdata); outdata = NULL;
  free(outcounts); outcounts = NULL;
  free(outdispls); outdispls = NULL;
  free(incounts); incounts = NULL;
  free(indispls); indispls = NULL;
  MPI_Type_free(&slot_data_type);

  /* Create the local part of the hash table. */
  slot_data* ht = (slot_data*)xmalloc(ht_my_size * sizeof(slot_data));
  for (i = ht_my_start; i < ht_my_end; ++i) {
    ht[HT_LOCAL(i)].index = (int64_t)(-1); /* Unused */
  }
  for (i = 0; i < total_incount; ++i) {
    int64_t index = indata[i].index, value = indata[i].value;
    assert (HT_OWNER(index) == rank);
    hashtable_insert(ht, ht_my_size, index, value, HT_LOCAL(index));
  }

  free(indata); indata = NULL;

  /* Make the local part of the result.  Most of the rest of this code is
   * similar to the shared-memory/XMT version above. */
  int64_t* result = (int64_t*)xmalloc(total_incount * sizeof(int64_t));
  *result_ptr = result;
  *result_size_ptr = total_incount;

  int64_t* bucket_counts = (int64_t*)xmalloc(ht_my_size * sizeof(int64_t));
  for (i = ht_my_start; i < ht_my_end; ++i) {
    /* Count all elements with same index. */
    bucket_counts[HT_LOCAL(i)] = hashtable_count_key(ht, ht_my_size, i, HT_LOCAL(i));
  }
  /* bucket_counts replaced by its prefix sum (start of each bucket in output array) */
  int64_t* bucket_starts_in_result = bucket_counts;
  int64_t running_sum = 0;
  for (i = 0; i < ht_my_size; ++i) {
    int64_t old_running_sum = running_sum;
    running_sum += bucket_counts[i];
    bucket_counts[i] = old_running_sum;
  }
  assert (running_sum == total_incount);
  bucket_counts = NULL;
  for (i = ht_my_start; i < ht_my_end; ++i) {
    int64_t result_start_idx = bucket_starts_in_result[HT_LOCAL(i)];
    int64_t* temp = result + result_start_idx;
    /* Gather up all elements with same key. */
    int64_t bi = (int64_t)hashtable_get_values(ht, ht_my_size, i, HT_LOCAL(i), temp);
    if (bi > 1) {
      /* Selection sort them (for consistency in parallel implementations). */
      selection_sort(temp, bi);
      /* Randomly permute them. */
      mrg_state new_st = *st;
      mrg_skip(&new_st, 1, i, 100);
      randomly_permute(temp, bi, &new_st);
    }
  }
  free(ht); ht = NULL;
  free(bucket_starts_in_result); bucket_starts_in_result = NULL;
}
#undef HT_OWNER
#undef HT_LOCAL


#endif /* PERMUTATION_GEN_H */
