/**
 * @file    check_graph500_gen.cpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Checks if graph500 graph generater works.
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/graph500/graph500Gen.hpp"

int main(int argc, char** argv)
{
  // Initialize the MPI library:
  MPI_Init(&argc, &argv);

  //Graph params
  size_t scale = 9;
  size_t edgefactor = 16;

  //Object of the graph500 generator class
  conn::graphGen::graph500Gen g;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<int64_t, int64_t> > edgeList;

  //Populate the edgeList
  g.populateEdgeList(edgeList, scale, edgefactor, conn::graphGen::graph500Gen::UNDIRECTED); 


  MPI_Finalize();
  return(0);
}
 


