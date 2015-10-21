/**
 * @file    time_ccl_coloring_undirected.cpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes connected components in the synthetic undirected kronecker graph
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/graph500/graph500Gen.hpp"
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

  //Know the rank and comm size
  int p, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  //Print mpi rank distribution
  mxx::print_node_distribution();

  LOG_IF(!rank, INFO) << "Code computes connected components using coloring in the undirected synthetic graph";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Computes connected components using coloring in the undirected synthetic graph");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("scale", "scale of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("edgefactor", "edgefactor of the graph", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!rank) cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  //Graph params
  uint8_t scale = std::stoi(cmd.optionValue("scale"));
  uint8_t edgefactor = std::stoi(cmd.optionValue("edgefactor"));

  //Object of the graph500 generator class
  conn::graphGen::graph500Gen g;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<int64_t, int64_t> > edgeList;

  //Populate the edgeList
  g.populateEdgeList(edgeList, scale, edgefactor, conn::graphGen::graph500Gen::UNDIRECTED, MPI_COMM_WORLD); 

  //Sum up the edge count across ranks
  auto totalEdgeCount = mxx::reduce(edgeList.size(), 0);
  LOG_IF(!rank, INFO) << "Total edge count is " << totalEdgeCount;

  //Compute connected components
  conn::coloring::ccl<uint32_t, uint64_t> cclInstance(edgeList);
  cclInstance.compute();

  MPI_Finalize();
  return(0);
}
 


