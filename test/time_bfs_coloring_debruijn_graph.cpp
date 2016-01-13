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
 * @file    time_bfs_coloring_debruijn_graph.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes weakly connected components in the de Bruijn graph
 *          using 1 BFS iteration and coloring 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/deBruijn/deBruijnGraphGen.hpp"
#include "graphGen/common/reduceIds.hpp"
#include "coloring/labelProp.hpp"
#include "bfs/bfsRunner.hpp"

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

  LOG_IF(!comm.rank(), INFO) << "Computing components for deBruijn graph using bfs and coloring";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Computing components for deBruijn graph using bfs and coloring");
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

  //Initialize the distributed vector for saving unique vertices
  std::vector<vertexIdType> uniqueVertexList;

  {
    mxx::section_timer timer(std::cerr, comm);

    //Object of the graph generator class
    conn::graphGen::deBruijnGraph g;

    //Populate the edgeList
    g.populateEdgeList(edgeList, fileName, comm); 

    timer.end_section("Graph generation completed");

    //Call the graph reducer function
    conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, comm);

    timer.end_section("Graph vertices reduction completed for BFS");
  }

  //Count of vertices in the reduced graph
  std::size_t nVertices = conn::graphGen::globalSizeOfVector(uniqueVertexList, comm);

  //Count of edges in the graph
  std::size_t nEdges = conn::graphGen::globalSizeOfVector(edgeList, comm);

  LOG_IF(!comm.rank(), INFO) << "Graph size : vertices -> " << nVertices << ", edges ->" << nEdges;


  /**
   * RUN BFS
   */

  //For saving the size of component discovered using BFS
  std::vector<std::size_t> componentCountsResult;

  {
    mxx::section_timer timer(std::cerr, comm);

    conn::bfs::bfsSupport<vertexIdType> bfsInstance(edgeList, nVertices, comm);

    //Run BFS once
    bfsInstance.runBFSIterations(1, componentCountsResult); 

    timer.end_section("BFS iteration completed");

    //Get the remaining edgeList
    bfsInstance.filterEdgeList();

    timer.end_section("Edgelist filtered for coloring");
  }


  LOG_IF(!comm.rank(), INFO) << "Component size traversed : " << componentCountsResult[0];

  //Sum up the edge count across ranks
  auto totalEdgeCount = mxx::reduce(edgeList.size(), 0, comm);
  LOG_IF(!comm.rank(), INFO) << "Total edge count remaining after BFS : " << totalEdgeCount;

  {
    mxx::section_timer timer(std::cerr, comm);

    //Compute connected components using coloring for the remaining graph
    //Among the subset of ranks which have non-zero count of edges

    comm.with_subset(edgeList.size() > 0, [&](const mxx::comm& comm){
        conn::coloring::ccl<vertexIdType> cclInstance(edgeList, comm);
        cclInstance.compute();

        auto countComponents = cclInstance.computeComponentCount();
        LOG_IF(!comm.rank(), INFO) << "Count of components -> " << countComponents;
    });

    timer.end_section("Coloring completed");
  }

  MPI_Finalize();
  return(0);
}
