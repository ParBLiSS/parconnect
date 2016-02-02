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
 * @file    time_ccl_coloring_debruijn_graph.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes connected components in a de Bruijn graph 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/deBruijn/deBruijnGraphGen.hpp"
#include "coloring/labelProp.hpp"
#include "bfs/bfsRunner.hpp"
#include "graphGen/common/reduceIds.hpp"

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

  /**
   * COMMAND LINE ARGUMENTS
   */

  LOG_IF(!comm.rank(), INFO) << "Computing components for de bruijn graph using coloring";

  //Parse command line arguments
  ArgvParser cmd;

  LOG_IF(!comm.rank(), INFO) << "Computing components for deBruijn graph using coloring";

  cmd.setIntroductoryDescription("Computing components for deBruijn graph using coloring");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("file", "fastq sequence file", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Input file
  std::string fileName = cmd.optionValue("file");
  LOG_IF(!comm.rank(), INFO) << "Input file -> " << fileName;

  /**
   * GENERATE GRAPH
   */
  using vertexIdType = int64_t;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  {
    mxx::section_timer timer(std::cerr, comm);

    //Object of the graph generator class
    conn::graphGen::deBruijnGraph g;

    //Populate the edgeList
    g.populateEdgeList(edgeList, fileName, comm); 

    timer.end_section("Graph generation completed");
  }

  //Count of edges in the graph
  std::size_t nEdges = conn::graphGen::globalSizeOfVector(edgeList, comm);

  LOG_IF(!comm.rank(), INFO) << "Graph size : edges ->" << nEdges;

  {
    mxx::section_timer timer(std::cerr, comm);

    //Perform the coloring
    conn::coloring::ccl<vertexIdType> cclInstance(edgeList, comm);
    cclInstance.compute();

    timer.end_section("Coloring completed");

    auto countComponents = cclInstance.computeComponentCount();
    LOG_IF(!comm.rank(), INFO) << "Count of components -> " << countComponents;
  }

  MPI_Finalize();
  return(0);
}
