/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef MAKE_GRAPH_H
#define MAKE_GRAPH_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

//Own includes
#include "graph_generator.h"
#include "permutation_gen.h"
#include "apply_permutation_mpi.h"
#include "scramble_edges.h"
#include "utils.h"

inline void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, const double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
  int64_t N, M;
  int rank, size;

  N = (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, log_numverts);
  M = desired_nedges;

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int64_t nedges = compute_edge_array_size(rank, size, M);
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* local_edges = (generated_edge*)xmalloc(nedges * sizeof(generated_edge));
#else
  int64_t* local_edges = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#endif

  //double start = MPI_Wtime();
  generate_kronecker(rank, size, seed, log_numverts, M, initiator, local_edges);
  //double gen_time = MPI_Wtime() - start;

  int64_t* local_vertex_perm = NULL;

  mrg_state state;
  mrg_seed(&state, seed);
  //start = MPI_Wtime();
  int64_t perm_local_size;
  rand_sort_mpi(MPI_COMM_WORLD, &state, N, &perm_local_size, &local_vertex_perm);
  //double perm_gen_time = MPI_Wtime() - start;

  /* Copy the edge endpoints into the result array if necessary. */
  int64_t* result;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
  for (i = 0; i < nedges; ++i) {
    if (local_edges[i].multiplicity != 0) {
      result[i * 2] = local_edges[i].src;
      result[i * 2 + 1] = local_edges[i].tgt;
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(local_edges); local_edges = NULL;
#else
  result = local_edges;
  *result_ptr = result;
  local_edges = NULL; /* Freed by caller */
#endif

  /* Apply vertex permutation to graph. */
  //start = MPI_Wtime();
  apply_permutation_mpi(MPI_COMM_WORLD, perm_local_size, local_vertex_perm, N, nedges, result);
  //double perm_apply_time = MPI_Wtime() - start;

  free(local_vertex_perm); local_vertex_perm = NULL;

  /* Randomly mix up the order of the edges. */
  //start = MPI_Wtime();
  int64_t* new_result;
  int64_t nedges_out;
  scramble_edges_mpi(MPI_COMM_WORLD, userseed1, userseed2, nedges, result, &nedges_out, &new_result);
  //double edge_scramble_time = MPI_Wtime() - start;

  free(result); result = NULL;

  *result_ptr = new_result;
  *nedges_ptr = nedges_out;

  /*if (rank == 0) {*/
    //fprintf(stdout, "unpermuted_graph_generation:    %f s\n", gen_time);
    //fprintf(stdout, "vertex_permutation_generation:  %f s\n", perm_gen_time);
    //fprintf(stdout, "vertex_permutation_application: %f s\n", perm_apply_time);
    //fprintf(stdout, "edge_scrambling:                %f s\n", edge_scramble_time);
  /*}*/
}

/* PRNG interface for implementations; takes seed in same format as given by
 * users, and creates a vector of doubles in a reproducible (and
 * random-access) way. */
inline void make_random_numbers(
       /* in */ int64_t nvalues    /* Number of values to generate */,
       /* in */ uint64_t userseed1 /* Arbitrary 64-bit seed value */,
       /* in */ uint64_t userseed2 /* Arbitrary 64-bit seed value */,
       /* in */ int64_t position   /* Start index in random number stream */,
       /* out */ double* result    /* Returned array of values */
) {
  int64_t i;
  uint_fast32_t seed[5];
  mrg_state st;
  make_mrg_seed(userseed1, userseed2, seed);


  mrg_seed(&st, seed);

  mrg_skip(&st, 2, 0, 2 * position); /* Each double takes two PRNG outputs */

  for (i = 0; i < nvalues; ++i) {
    result[i] = mrg_get_double_orig(&st);
  }
}


#endif /* MAKE_GRAPH_H */
