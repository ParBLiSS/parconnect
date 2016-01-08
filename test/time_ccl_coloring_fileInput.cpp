/**
 * @file    time_ccl_coloring_fileInput.cpp
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

  LOG_IF(!comm.rank(), INFO) << "Computing components for general graph file input using coloring";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Computing components for general graph file input using coloring");
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

  {
    mxx::section_timer timer(std::cerr, comm);

    //Object of the graph generator class
    conn::graphGen::GraphFileParser<char *, vertexIdType> g(edgeList, addReverse, comm);

    //Populate the edgeList
    g.populateEdgeList(fileName);

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
