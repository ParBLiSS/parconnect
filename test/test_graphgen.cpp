/**
 * @file    test_graphgen.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   GTest Unit Tests related to the graph generation
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */


#include <mpi.h>
#include <algorithm>

//Own includes
#include "utils/logging.hpp"
#include "utils/prettyprint.hpp"
#include "graphGen/common/reduceIds.hpp"
#include "graphGen/graph500/graph500Gen.hpp"

//External includes
#include "gtest.h"
#include "mxx/comm.hpp"

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

  //We have added 8 unique vertices on each rank
  std::size_t totalUniqueVerticesGlobal = 8 * c.size();

  std::random_shuffle(edgeList.begin(), edgeList.end());

  //Initialize the distributed vector for saving unique vertices
  std::vector<nodeIdType> uniqueVertexList;

  //Call the function
  conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, c);

  //Total size of uniqueVertexList or the count of unique vertices in the whole graph
  auto globalSizeUnique = conn::graphGen::globalSizeOfVector(uniqueVertexList, c);

  //Each rank should contain following edge like (0,1) on proc 0 and (8, 9) on proc 1
  //This is assuming the edges were sorted during the reduction
  auto edge = std::make_pair((nodeIdType)0 + 8*c.rank(), (nodeIdType)1 + 8*c.rank());

  //Edge should exist on this rank
  bool edgeContained = (std::find(edgeList.begin(), edgeList.end(), edge) != edgeList.end());

  ASSERT_EQ(totalUniqueVerticesGlobal, globalSizeUnique);
  ASSERT_TRUE(edgeContained);
}

/**
 * @brief       Run graph500 generator and graph reducer
 *              Fetch the value of vertex count in the graph
 *              Check all the vertex ids in the graph are less than above value
 */
TEST(graphGen, reduceIdgraph500) {

  mxx::comm c = mxx::comm();

  using nodeIdType = int64_t;

  std::vector< std::pair<nodeIdType, nodeIdType> > edgeList;

  //bool value, true if all the vertex ids are less than vertex count in the graph
  bool vertexIdCheck = true;

  {
    //Object of the graph500 generator class
    conn::graphGen::graph500Gen g;

    int scale = 11, edgefactor = 16;

    //Populate the edgeList (scale = 10, edgefactor = 16, Undirected for adding edges both ways)
    g.populateEdgeList(edgeList, scale, edgefactor, c); 

    //Initialize the distributed vector for saving unique vertices
    std::vector<nodeIdType> uniqueVertexList;

    //Call the function
    conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, c);

    //Total size of uniqueVertexList
    auto globalSizeUnique = conn::graphGen::globalSizeOfVector(uniqueVertexList, c);

    LOG_IF(!c.rank(), INFO ) << "Graph500 scale = " << scale; 
    LOG_IF(!c.rank(), INFO ) << "Unique vertex count = " << globalSizeUnique; 

    for(auto &e : edgeList)
    {
      vertexIdCheck = vertexIdCheck && (std::max(e.first, e.second) < globalSizeUnique);
    }
  }

  ASSERT_TRUE(vertexIdCheck);
}

