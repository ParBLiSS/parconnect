/**
 * @file    time_ccl_coloring_undirected_chains.cpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks the performance of coloring approach on undirected chain graphs
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */
//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/undirectedChain/undirectedChainGen.hpp"
#include "coloring/labelProp.hpp"
#include "utils/logging.hpp"
#include "utils/argvparser.hpp"

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


  LOG_IF(!comm.rank(), INFO) << "Code computes connected components using coloring on the chains of different lengths";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Code computes connected components using coloring on the chains of different lengths");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("startLength", "length of the smallest chain to run on", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("scaleFactor", "chain length is increased by this factor", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("scaleUpSteps", "Number of times you wish to scale up the graph size, > 0", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);


  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  std::size_t startLength = std::stoi(cmd.optionValue("startLength"));
  uint8_t scaleFactor = std::stoi(cmd.optionValue("scaleFactor"));
  uint8_t scaleUpSteps = std::stoi(cmd.optionValue("scaleUpSteps"));

  //Have atleast 1 graph
  assert(scaleUpSteps > 0);

  //Object of the graph500 generator class
  conn::graphGen::UndirectedChainGen g;

  using nodeIdType = uint64_t;

  std::size_t graphSize = startLength;
  for(int i = 0; i < scaleUpSteps; i++)
  {
    //Declare a edgeList vector to save edges
    std::vector< std::pair<nodeIdType, nodeIdType> > edgeList;

    //Create the edgeList
    g.populateEdgeList(edgeList, graphSize, conn::graphGen::UndirectedChainGen::LOWTOHIGH_IDS, comm); 

    LOG_IF(!comm.rank(), INFO) << "Chain size " << graphSize;

    //Perform the coloring
    conn::coloring::ccl<nodeIdType> cclInstance(edgeList, comm);
    cclInstance.compute();

    //Increase the graph size (or chain length)
    graphSize *= scaleFactor;
  }


  MPI_Finalize();
  return(0);
}


