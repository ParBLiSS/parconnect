/**
 * @file    test_graphgen.cpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   GTest Unit Tests related to the graph generation
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */


#include <mpi.h>
#include <algorithm>

#include "utils/logging.hpp"
#include "utils/prettyprint.hpp"
#include "mxx/comm.hpp"
#include "gtest.h"

#include "graphGen/common/reduceIds.hpp"

INITIALIZE_EASYLOGGINGPP

/**
 * @brief       checking edge id reduction for a small example graph
 */
TEST(graphGen, reduceIdSmallGraph) {

  mxx::comm c = mxx::comm();

  using nodeIdType = int64_t;

  std::vector< std::pair<nodeIdType, nodeIdType> > edgeList;

  //Start adding the edges
  //Example : For 4 processes, offsets will be 400, 300, 200, 100
  nodeIdType offset = 100 * (c.size() - c.rank());

  //Each rank adds 5 bidirectional edges
  //8 unique vertices produced by each rank
  for(int i = 1; i < 6; i ++)
  {
    edgeList.emplace_back(3*i + offset, 6*i + offset);
    edgeList.emplace_back(6*i + offset, 3*i + offset);
  }

  std::size_t totalUniqueVerticesGlobal = 8 * c.size();

  std::random_shuffle(edgeList.begin(), edgeList.end());

  //Initialize the distributed vector for saving unique vertices
  std::vector<nodeIdType> uniqueVertexList;

  //Call the function
  conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, c);

  //Total size of uniqueVertexList
  auto globalSizeUnique = conn::graphGen::globalSizeOfVector(uniqueVertexList, c);

  //Each rank should contain following edge
  auto edge = std::make_pair((nodeIdType)1 + 8*c.rank(), (nodeIdType)2 + 8*c.rank());

  //Edge should exist on this rank
  bool edgeContained = (std::find(edgeList.begin(), edgeList.end(), edge) != edgeList.end());

  ASSERT_EQ(totalUniqueVerticesGlobal, globalSizeUnique);
  ASSERT_TRUE(edgeContained);
}


