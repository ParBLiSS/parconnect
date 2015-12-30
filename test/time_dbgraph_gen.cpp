/**
 * @file    time_dbgraph_gen.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks and times the de bruijn graph construction
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/deBruijn/deBruijnGraphGen.hpp"
#include "utils/logging.hpp"
#include "utils/argvparser.hpp"

//External includes
#include "mxx/reduction.hpp"
#include "mxx/utils.hpp"
#include "debruijn/de_bruijn_construct_engine.hpp"
#include "debruijn/de_bruijn_nodes_distributed.hpp"

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

  LOG_IF(!comm.rank(), INFO) << "Code to time de Bruijn graph construction";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("checks and times the de Bruijn graph construction");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("file", "fastq sequence file", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

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

  //Declare a edgeList vector to save edges
  std::vector< std::pair<std::size_t, std::size_t> > edgeList;

  {
    //Object of the graph generator class
    conn::graphGen::deBruijnGraph g;

    //Populate the edgeList
    g.populateEdgeList(edgeList, fileName, comm); 
  }

  //Sum up the edge count across ranks
  auto totalEdgeCount = mxx::reduce(edgeList.size(), 0, comm);
  LOG_IF(!comm.rank(), INFO) << "Total edge count is " << totalEdgeCount;

  MPI_Finalize();
  return(0);
}
