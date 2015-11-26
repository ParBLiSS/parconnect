/**
 * @file    bfsSupport.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Finds connected components using a BFS Run
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */


#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/graph500/graph500Gen.hpp"
#include "CombBLAS/CombBLAS.h"
#include "utils/logging.hpp"

void runBFSTrial(int scale, int edgefactor, mxx::comm &comm)
{
  //Matrix type, to store the adjacency matrix (bool values)
	typedef SpParMat < int64_t, bool, SpDCCols<int64_t,bool> > MatTypeForBool;

  //Matrix type, to store the adjacency matrix (int values)
	typedef SpParMat < int64_t, int64_t, SpDCCols<int64_t,int64_t> > MatTypeForInt;

  // degrees of vertices (including multi-edges and self-loops)
  FullyDistVec<int64_t, int64_t> degrees;	

  //List of edges, distributed in 1D fashion
  DistEdgeList<int64_t> * DEL = new DistEdgeList<int64_t>();

  //Optimization buffer (not sure what it is needed for)
  OptBuf<int32_t, int64_t> optbuf;

  //Copy the edges from our genarator to DEL
  {
    //Object of the graph500 generator class
    conn::graphGen::graph500Gen g;

    //graph500 generator only uses int64_t
    using nodeIdType = int64_t;

    //Declare a edgeList vector to save edges
    std::vector< std::pair<nodeIdType, nodeIdType> > edgeList;

    //Populate the edgeList
    g.populateEdgeList(edgeList, scale, edgefactor, conn::graphGen::graph500Gen::UNDIRECTED, comm); 

    //Copy our edgeList to CombBLAS format of edgeList
    DEL->GenGraph500Data(edgeList, scale);

  }//g is destructed after the scope ends

  //Sparse matrix for adjacency list from list of edges
  //TODO : Avoid copy to DEL and directly use native edgeList
  MatTypeForInt * G = new MatTypeForInt(*DEL, false); 
  delete DEL;	// free memory

  //Calculate degrees
  G->Reduce(degrees, Row, plus<int64_t>(), static_cast<int64_t>(0));	// Identity is 0 

  //Count of vertices
  double nver = (double) degrees.TotalLength();

  //Now represent the adj matrix in the boolean format
  MatTypeForBool A =  MatTypeForBool(*G);			// Convert to Boolean

  A.OptimizeForGraph500(optbuf);		          //Taken from TopDownBFS code, some kind of optimization for graph500 graphs is expected here

  //Frontier
  FullyDistSpVec<int64_t, int64_t> fringe(A.getcommgrid(), A.getncol());	// numerical values are stored 0-based

  //Parent array (acts as a list of vertices in a component for us)
  FullyDistVec<int64_t, int64_t> parents(A.getcommgrid(), A.getncol(), (int64_t) -1);	// numerical values are stored 0-based

  MPI_Barrier(MPI_COMM_WORLD);

  //Source of BFS set to 0
  fringe.SetElement(0 ,0);
  parents.SetElement(0, 0);

  int64_t fringe_size = fringe.getnnz();
  fringe.Apply(myset<int64_t>(1));

  int iterations = 0;
  while(fringe_size > 0) 
  {
    // Top-down
    fringe.setNumToInd();

    //Matrix multiplication
    fringe = SpMV(A, fringe,optbuf);

    //Remove elements from frontier that were already visited before
    fringe = EWiseMult(fringe, parents, true, (int64_t) -1);

    //Update parents array using fringe
    parents.Set(fringe);

    fringe.Apply(myset<int64_t>(1));

    //Global count of elements in the frontier
    fringe_size = fringe.getnnz();

    iterations++;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  FullyDistSpVec<int64_t, int64_t> parentsp = parents.Find(bind2nd(greater<int64_t>(), -1));
  parentsp.Apply(myset<int64_t>(1));
  auto countVertices = parentsp.Reduce(plus<int64_t>(), (int64_t) 0);

  LOG_IF(comm.rank() == 0, INFO) << "Number of iterations: " << iterations;
  LOG_IF(comm.rank() == 0, INFO) << "Number of vertices found: " << countVertices;
}

