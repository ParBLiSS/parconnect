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
 * @file    time_ccl_coloring_graph500.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes connected components in the synthetic undirected kronecker graph
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/graph500/graph500Gen.hpp"
#include "graphGen/common/reduceIds.hpp"
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

  LOG_IF(!comm.rank(), INFO) << "Code computes connected components using coloring in the undirected synthetic graph";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Computes connected components using coloring in the undirected synthetic graph");
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


  //graph500 generator only uses int64_t
  using nodeIdType = int64_t;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<nodeIdType, nodeIdType> > edgeList;

  {
    mxx::section_timer timer(std::cerr, comm);

    //Object of the graph500 generator class
    conn::graphGen::Graph500Gen g;

    //Populate the edgeList
    g.populateEdgeList(edgeList, scale, edgefactor, comm); 

    timer.end_section("Graph generation completed");
  }

  //Count of edges in the graph
  std::size_t nEdges = conn::graphGen::globalSizeOfVector(edgeList, comm);

  LOG_IF(!comm.rank(), INFO) << "Graph size : edges [when every (u,v) has a (v,u) edge] -> " << nEdges;

  {
    mxx::section_timer timer(std::cerr, comm);

    //Compute connected components
    conn::coloring::ccl<nodeIdType> cclInstance(edgeList, comm);
    cclInstance.compute();

    timer.end_section("Coloring completed");

    auto countComponents = cclInstance.computeComponentCount();
    LOG_IF(!comm.rank(), INFO) << "Count of components -> " << countComponents;
  }

  MPI_Finalize();
  return(0);
}
