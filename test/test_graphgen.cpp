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
#include "utils/commonfuncs.hpp"
#include "graphGen/common/reduceIds.hpp"
#include "graphGen/graph500/graph500Gen.hpp"
#include "graphGen/fileIO/graphReader.hpp"

//External includes
#include "extutils/logging.hpp"
#include "gtest.h"
#include "mxx/comm.hpp"

//To get the path to src/test/data folder 
#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)
#define PROJECT_TEST_DATA_FOLDER STR(PROJECT_TEST_DATA_DIR)

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
 *              Fetch the value of vertex count ie |V| in the graph
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
    conn::graphGen::Graph500Gen g;

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

/*
 * @brief   Test the correctness of parallel file io graph reader
 *          Supply the file as input containing a graph of directed
 *          chain {1-2-3...1201}
 *          FILE : src/test/data/graphDirChain.txt
 *          Read the graph and test the presence of all the edges we expect
 */
TEST(graphGen, graphFileIO) {

  mxx::comm comm = mxx::comm();

  //File name
  //Make sure this path to file is correct 
  std::string fileName = PROJECT_TEST_DATA_FOLDER;
  fileName = fileName + "/graphDirChain.txt";

  //Above graph is a directed chain, we need to include reverse own our own
  bool addreverseEdgeWhileParsing = true;

  using vertexIdType = int64_t;

  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  {
    conn::graphGen::GraphFileParser<char *, vertexIdType> g(edgeList, addreverseEdgeWhileParsing, fileName, comm);

    g.populateEdgeList();
  }

  //Gather complete edgeList on rank 0
  auto fullEdgeList = mxx::gatherv(edgeList, 0, comm);

  if(!comm.rank())
  {
    const int SRC = 0, DEST = 1;
    std::sort(fullEdgeList.begin(), fullEdgeList.end(), conn::utils::TpleComp2Layers<SRC, DEST>());

    ASSERT_EQ(fullEdgeList.size(), 2400);

    //Check for the presence of 1-2. 2-3, 3-4...1200-1201 at index 0, 2, ...
    int j = 1;
    for(int i = 0; i < fullEdgeList.size(); i += 2)
    {
      ASSERT_EQ(fullEdgeList[i].first, j);
      ASSERT_EQ(fullEdgeList[i].second, j + 1);

      j++;
    }


    //and then check for the presence of 2-1 3-2 ... 1201-1200 at index 1, 3 ...
    j = 1;
    for(int i = 1; i < fullEdgeList.size(); i += 2)
    {
      ASSERT_EQ(fullEdgeList[i].first, j + 1);
      ASSERT_EQ(fullEdgeList[i].second, j);

      j++;
    }

  }
}
