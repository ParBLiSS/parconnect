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
