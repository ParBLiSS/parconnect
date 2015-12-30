/**
 * @file    time_bfs_graph500.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   BFS iterations on kroncker graph
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

//External includes
#include "mxx/utils.hpp"
#include "mxx/timer.hpp"

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

  LOG_IF(!comm.rank(), INFO) << "BFS runs on Kronecker graph";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("BFS runs on Kronecker graph");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("scale", "scale of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("iter", "max count of BFS iterations, default = INF", ArgvParser::OptionRequiresValue);
  cmd.defineOption("edgefactor", "edgefactor of the graph, default = 16", ArgvParser::OptionRequiresValue);

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

  std::size_t iterBound = std::numeric_limits<std::size_t>::max();
  if(cmd.foundOption("iter")) {
    std::stringstream sstream(cmd.optionValue("iter"));
    sstream >> iterBound;
    LOG_IF(!comm.rank(), INFO) << "BFS iterations count limit -> " << iterBound;
  }
  else
    LOG_IF(!comm.rank(), INFO) << "BFS iterations count limit -> No limit";


  int edgefactor = 16;
  if(cmd.foundOption("edgefactor")) {
    edgefactor = std::stoi(cmd.optionValue("edgefactor"));
  }
  LOG_IF(!comm.rank(), INFO) << "Edgefactor -> " << edgefactor;

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
    mxx::section_timer timer(std::cerr, comm);

    //Object of the graph500 generator class
    conn::graphGen::Graph500Gen g;

    //Populate the edgeList, using undirected mode that includes the edge and its reverse
    g.populateEdgeList(edgeList, scale, edgefactor, comm); 

    timer.end_section("Graph generation completed");

    //Call the graph reducer function
    conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, comm);

    timer.end_section("Graph vertices reduction completed for BFS");
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
    mxx::section_timer timer(std::cerr, comm);

    conn::bfs::bfsSupport<vertexIdType> bfsInstance(edgeList, nVertices, comm);

    //Run BFS "iterBound" times
    bfsInstance.runBFSIterations(iterBound, componentCountsResult); 

    timer.end_section("BFS iterations completed");
  }


  std::size_t totalVerticesTraversed = std::accumulate(componentCountsResult.begin(), componentCountsResult.end(), 0);

  std::size_t verticesInLargestComponent = *std::max_element(componentCountsResult.begin(), componentCountsResult.end());

  LOG_IF(!comm.rank(), INFO) << "Count of vertices traversed :" << totalVerticesTraversed;
  LOG_IF(!comm.rank(), INFO) << "Count of components :" << componentCountsResult.size();
  LOG_IF(!comm.rank(), INFO) << "Largest Component vertex count :" << verticesInLargestComponent;

  MPI_Finalize();
  return(0);
}
