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
#include "bfs/bfsSupport.hpp"
#include "utils/logging.hpp"
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

  LOG_IF(!comm.rank(), INFO) << "Single BFS run on Kronecker graph";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Single BFS run on Kronecker graph");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("scale", "scale of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("edgefactor", "edgefactor of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  uint8_t scale = std::stoi(cmd.optionValue("scale"));
  uint8_t edgefactor = std::stoi(cmd.optionValue("edgefactor"));

  runBFSTrial(scale, edgefactor, comm);

  MPI_Finalize();
  return(0);
}
