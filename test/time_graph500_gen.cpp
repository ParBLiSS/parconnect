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
 * @file    time_graph500_gen.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks if graph500 graph generater works.
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/graph500/graph500Gen.hpp"
#include "graphGen/common/utils.hpp"

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
  mxx::comm comm = mxx::comm();

  //Print mpi rank distribution
  mxx::print_node_distribution();

  LOG_IF(!comm.rank(), INFO) << "Code to check kronecker graph generation";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("checks the kronecker graph generation");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("scale", "scale of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("edgefactor", "edgefactor of the graph", ArgvParser::OptionRequiresValue);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  int scale = std::stoi(cmd.optionValue("scale"));
  LOG_IF(!comm.rank(), INFO) << "scale -> " << scale;

  int edgefactor = 16;
  if(cmd.foundOption("edgefactor")) {
    edgefactor = std::stoi(cmd.optionValue("edgefactor"));
  }
  LOG_IF(!comm.rank(), INFO) << "Edgefactor -> " << edgefactor;

  //Object of the graph500 generator class
  conn::graphGen::Graph500Gen g;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<int64_t, int64_t> > edgeList;

  //Populate the edgeList
  g.populateEdgeList(edgeList, scale, edgefactor, comm); 

  //Sum up the edge count across ranks
  auto totalEdgeCount = mxx::reduce(edgeList.size(), 0, comm);
  LOG_IF(!comm.rank(), INFO) << "Total edge count is " << totalEdgeCount;

  //Check if each edge is added in both directions
  bool graphFormatCheck = conn::graphGen::checkEdgeBidirectionality(edgeList, comm);

  LOG_IF(!comm.rank() && graphFormatCheck, INFO) <<  "Graph format check passed";
  LOG_IF(!comm.rank() && !graphFormatCheck, INFO) << "Graph format check failed";

  MPI_Finalize();
  return(0);
}
