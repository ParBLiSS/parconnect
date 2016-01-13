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
 * @file    time_ccl_coloring_chains.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks the performance of coloring approach on undirected chain graphs
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */
//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/undirectedChain/undirectedChainGen.hpp"
#include "coloring/labelProp.hpp"

//External includes
#include "extutils/logging.hpp"
#include "extutils/argvparser.hpp"
#include "mxx/reduction.hpp"
#include "mxx/utils.hpp"

INITIALIZE_EASYLOGGINGPP
using namespace std;
using namespace CommandLineProcessing;

int main(int argc, char** argv)
{
  // Initialize the MPI library:
  MPI_Init(&argc, &argv);

  //Initialize the communicator
  mxx::comm comm;

  //Print mpi rank distribution
  mxx::print_node_distribution();


  LOG_IF(!comm.rank(), INFO) << "Code computes connected components using coloring on the chain of given length";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Code computes connected components using coloring on the chain of given length");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("length", "length of the chain (# nodes)", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);


  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  std::size_t length = std::stoi(cmd.optionValue("length"));

  //Object of the graph500 generator class
  conn::graphGen::UndirectedChainGen g;

  using nodeIdType = uint64_t;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<nodeIdType, nodeIdType> > edgeList;

  //Create the edgeList
  g.populateEdgeList(edgeList, length, conn::graphGen::UndirectedChainGen::LOWTOHIGH_IDS, comm); 

  LOG_IF(!comm.rank(), INFO) << "Chain size " << length;

  {
    //Perform the coloring
    conn::coloring::ccl<nodeIdType> cclInstance(edgeList, comm);
    cclInstance.compute();
  }


  MPI_Finalize();
  return(0);
}
