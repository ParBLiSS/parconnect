/**
 * @file    time_bfs_coloring_fileInput.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes connected components in a graph given as a general file input 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/fileIO/graphReader.hpp"
#include "graphGen/common/reduceIds.hpp"
#include "coloring/labelProp.hpp"
#include "bfs/bfsRunner.hpp"
#include "utils/logging.hpp"
#include "utils/argvparser.hpp"

//External includes
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

  LOG_IF(!comm.rank(), INFO) << "Computing components for general graph file input using bfs and coloring";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Computing components for general graph file input using bfs and coloring");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("file", "input file with edges written along the rows", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("addreverse", "(y/n) y implies reverse of each edge will also be added", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);


  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Input file
  std::string fileName = cmd.optionValue("file");
  bool addReverse = !cmd.optionValue("addreverse").compare("y");  //compare returns 0 iff equal

  LOG_IF(!comm.rank(), INFO) << "Input file -> " << fileName;
  if(addReverse)  LOG_IF(!comm.rank(), INFO) << "Reverse of each edge will be included";
  if(!addReverse)  LOG_IF(!comm.rank(), INFO) << "Reverse of each edge will not be included";


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
    conn::graphGen::GraphFileParser<char *, vertexIdType> g(edgeList, addReverse, comm);

    //Populate the edgeList
    g.populateEdgeList(fileName);

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
