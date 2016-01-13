/*
 * Copyright 2016 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file    check_ccl_coloring_undirected.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   GTest Unit Tests for connected component labeling using coloring
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#include <mpi.h>

//Own includes
#include "coloring/labelProp.hpp"

//External includes
#include "mxx/comm.hpp"
#include "gtest.h"

INITIALIZE_EASYLOGGINGPP

/**
 * @brief       coloring of undirected graph with a small chain
 * @details     builds a small undirected chain, 
 *              test if program returns 1 as the component count
 */
TEST(connColoring, smallUndirectedChain) {

  mxx::comm c = mxx::comm();

  //Declare a edgeList vector to save edges
  std::vector< std::pair<int64_t, int64_t> > edgeList;

  using nodeIdType = int64_t;
  using pIdType = uint32_t;

  //Start adding the edges
  if (c.rank() == 0) {

    //Chain (chain 1-2-...1000)
    for(int i = 1; i < 1000 ; i++)
    {
      edgeList.emplace_back(i, i+1);
      edgeList.emplace_back(i+1, i);
    }
  }

  std::random_shuffle(edgeList.begin(), edgeList.end());
  conn::coloring::ccl<nodeIdType> cclInstance(edgeList, c);
  cclInstance.compute();
  auto component_count = cclInstance.computeComponentCount();
  ASSERT_EQ(1, component_count);
}

/**
 * @brief       coloring of undirected graph with 3 components
 * @details     builds a small test graph with three components, 
 *              test if program returns 3 as the component count
 */
TEST(connColoring, smallUndirected) {

  mxx::comm c = mxx::comm();

  //Declare a edgeList vector to save edges
  std::vector< std::pair<int64_t, int64_t> > edgeList;

  using nodeIdType = int64_t;

  //Start adding the edges
  if (c.rank() == 0) {

    //First component (2,3,4,11)
    //2--11
    edgeList.emplace_back(2,11);
    edgeList.emplace_back(11,2);

    //2--3
    edgeList.emplace_back(2,3);
    edgeList.emplace_back(3,2);

    //2--4
    edgeList.emplace_back(2,4);
    edgeList.emplace_back(4,2);

    //3--4
    edgeList.emplace_back(3,4);
    edgeList.emplace_back(4,3);

    //Second component (5,6,8,10)
    //5--6
    edgeList.emplace_back(5,6);
    edgeList.emplace_back(6,5);

    //5--8
    edgeList.emplace_back(5,8);
    edgeList.emplace_back(8,5);

    //6--10
    edgeList.emplace_back(6,10);
    edgeList.emplace_back(10,6);

    //6--8
    edgeList.emplace_back(6,8);
    edgeList.emplace_back(8,6);

    //Third component (50,51,52)
    for(int i = 50; i < 52 ; i++)
    {
      edgeList.emplace_back(i, i+1);
      edgeList.emplace_back(i+1, i);
    }
  }

  std::random_shuffle(edgeList.begin(), edgeList.end());

  //Since we have a graph of small size, use less processes
  c.with_subset(c.rank() < 4, [&](const mxx::comm& comm){
      conn::coloring::ccl<nodeIdType> cclInstance(edgeList, comm);
      cclInstance.compute();
      auto component_count = cclInstance.computeComponentCount();
      ASSERT_EQ(3, component_count);
      });
}

/**
 * @brief       coloring of undirected graph with 3 components
 * @details     builds a medium test graph with three components, 
 *              test if program returns 3 as the component count
 */
TEST(connColoring, mediumUndirected) {

  mxx::comm c = mxx::comm();

  //Declare a edgeList vector to save edges
  std::vector< std::pair<uint64_t, uint64_t> > edgeList;

  //Start adding the edges
  if (c.rank() == 0) {

    //First component (2,3,4,11)
    //2--11
    edgeList.emplace_back(2,11);
    edgeList.emplace_back(11,2);

    //2--3
    edgeList.emplace_back(2,3);
    edgeList.emplace_back(3,2);

    //2--4
    edgeList.emplace_back(2,4);
    edgeList.emplace_back(4,2);

    //3--4
    edgeList.emplace_back(3,4);
    edgeList.emplace_back(4,3);

    //Second component (5,6,8,10)
    //5--6
    edgeList.emplace_back(5,6);
    edgeList.emplace_back(6,5);

    //5--8
    edgeList.emplace_back(5,8);
    edgeList.emplace_back(8,5);

    //6--10
    edgeList.emplace_back(6,10);
    edgeList.emplace_back(10,6);

    //6--8
    edgeList.emplace_back(6,8);
    edgeList.emplace_back(8,6);

    //Third component (chain 50-51-...1000)
    for(int i = 50; i < 1000 ; i++)
    {
      edgeList.emplace_back(i, i+1);
      edgeList.emplace_back(i+1, i);
    }
  }

  std::random_shuffle(edgeList.begin(), edgeList.end());
  conn::coloring::ccl<> cclInstance(edgeList, c);
  cclInstance.compute();
  auto component_count = cclInstance.computeComponentCount();
  ASSERT_EQ(3, component_count);
}

