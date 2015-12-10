/**
 * @file    time_single_bfs_run.cpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Single bfs run on kroncker graph
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/graph500/graph500Gen.hpp"
#include "graphGen/common/reduceIds.hpp"
#include "bfs/bfsRunner.hpp"

#include "utils/logging.hpp"
#include "utils/prettyprint.hpp"
#include "utils/argvparser.hpp"

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

  LOG_IF(!comm.rank(), INFO) << "Single BFS run on Kronecker graph";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Single BFS run on Kronecker graph");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("scale", "scale of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("edgefactor", "edgefactor of the graph, default = 16", ArgvParser::OptionRequiresValue);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  uint8_t scale = std::stoi(cmd.optionValue("scale"));
  uint8_t edgefactor = 16;
  if(cmd.foundOption("edgefactor")) {
    uint8_t edgefactor = std::stoi(cmd.optionValue("edgefactor"));
  }

  /**
   * GENERATE GRAPH
   */


  //graph500 generator only uses int64_t, which is fine
  using vertexIdType = int64_t;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  //Initialize the distributed vector for saving unique vertices
  std::vector<vertexIdType> uniqueVertexList;

  {
    //Object of the graph500 generator class
    conn::graphGen::graph500Gen g;

    //Populate the edgeList, using undirected mode that includes the edge and its reverse
    g.populateEdgeList(edgeList, scale, edgefactor, conn::graphGen::graph500Gen::UNDIRECTED, comm); 

    //Call the graph reducer function
    conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, comm);
  }

  //Count of vertices in the reduced graph
  std::size_t nVertices = conn::graphGen::globalSizeOfVector(uniqueVertexList, comm);

  //Count of edges in the graph
  std::size_t nEdges = conn::graphGen::globalSizeOfVector(edgeList, comm);

  LOG_IF(!comm.rank(), INFO) << "Graph size : vertices -> " << nVertices << ", edges [when every (u,v) has a (v,u) edge] -> " << nEdges;

  /**
   * RUN BFS
   */

  //For saving the size of component discovered using BFS
  std::vector<std::size_t> componentCountsResult;

  {
    conn::bfs::bfsSupport<vertexIdType> bfsInstance(edgeList, nVertices, comm);

    //Run BFS once
    bfsInstance.runBFSIterations(1, componentCountsResult); 
  }

  LOG_IF(!comm.rank(), INFO) << "Component size traversed :" << componentCountsResult[0];

  MPI_Finalize();
  return(0);
}
