/**
 * @file    time_graphFileIO_gen.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks and times the graph construction from a generic file with edges
 *          Uses parallel IO supported by BLISS 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/fileIO/graphReader.hpp"
#include "graphGen/common/utils.hpp"
#include "utils/logging.hpp"
#include "utils/argvparser.hpp"

//External includes
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

  LOG_IF(!comm.rank(), INFO) << "Code for the graph construction by reading a file";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("constructs the graph by parallel reading of file");
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

  //Graph params
  std::string fileName = cmd.optionValue("file");
  bool addReverse = !cmd.optionValue("addreverse").compare("y");  //compare returns 0 iff equal

  LOG_IF(!comm.rank(), INFO) << "Input file -> " << fileName;
  if(addReverse)  LOG_IF(!comm.rank(), INFO) << "Reverse of each edge will be included";
  if(!addReverse)  LOG_IF(!comm.rank(), INFO) << "Reverse of each edge will not be included";

  typedef int64_t vertexIdType;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  {
    conn::graphGen::GraphFileParser<char *, vertexIdType> g(edgeList, addReverse, fileName, comm);

    g.populateEdgeList();
  }

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
