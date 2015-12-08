/**
 * @file    test_bfsRunner.cpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks the correctness of BFS runs 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>

//Own includes
#include "gtest.h"
#include "bfs/bfsRunner.hpp"

#include "mxx/utils.hpp"

INITIALIZE_EASYLOGGINGPP

/**
 * @brief     Each rank initializes a chain graph of length 50, 
 *            and we run BFS over that p times
 */
TEST(bfsRunCheck, multipleUndirectedChains) {

  mxx::comm comm = mxx::comm();

  //Type to use for vertices
  using vertexIdType = int64_t;

  //Distributed edge list
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  std::size_t offset = 50*comm.rank();

  //Each rank builds undirected chain of length 50
  //[0---49], [50---99] and so on
  for(int i = 0; i < 49; i ++)
  {
    edgeList.emplace_back(i    +offset, i+1  +offset);
    edgeList.emplace_back(i+1  +offset, i    +offset);
  }

  //Count of vertices
  std::size_t nVertices = 50*comm.size();

  {
    conn::bfs::bfsSupport<vertexIdType> bfsInstance(edgeList, nVertices, comm);
    std::vector<std::size_t> componentCountsResult;
    bfsInstance.runBFSIterations(comm.size(), componentCountsResult); 

    std::vector<std::size_t> componentCountsExpected(comm.size(), 50);

    //Check the size of the vector returned
    bool check1 = (componentCountsResult.size() == comm.size());

    //Check the values of the vector returned
    bool check2 = std::equal(componentCountsResult.begin(), componentCountsResult.end(), componentCountsExpected.begin());

    ASSERT_EQ(check1, true); 
    ASSERT_EQ(check2, true);
  }

}

