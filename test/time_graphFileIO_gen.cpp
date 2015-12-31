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

  LOG_IF(!comm.rank(), INFO) << "Code to time graph construction by reading a file";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("constructs the graph by parallel reading of file");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("file", "input file with edges written along the rows", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  std::string fileName = cmd.optionValue("file");
  LOG_IF(!comm.rank(), INFO) << "Input file -> " << fileName;

  typedef std::size_t vertexIdType;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  {
    conn::graphGen::GraphFileParser<char *, vertexIdType> g(edgeList, comm);

    g.populateEdgeList(fileName);
  }

  //Sum up the edge count across ranks
  auto totalEdgeCount = mxx::reduce(edgeList.size(), 0, comm);
  LOG_IF(!comm.rank(), INFO) << "Total edge count is " << totalEdgeCount;

  MPI_Finalize();
  return(0);
}
