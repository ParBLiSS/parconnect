#ifndef BFS_FUNCT2_HPP 
#define BFS_FUNCT2_HPP

template <typename EdgeArrayType>
void bfsComponentFinder2(EdgeArrayType &edgeList, std::size_t vertexCount)
{
  typedef SpParMat < int64_t, bool, SpDCCols<int64_t,bool> > PSpMat_Bool;
  typedef SpParMat < int64_t, int64_t, SpDCCols<int64_t,int64_t> > PSpMat_Int64;

  // Declare objects
  PSpMat_Bool Aeff;	
  OptBuf<int32_t, int64_t> optbuf;	// let indices be 32-bits
  //double initiator[4] = {.57, .19, .19, .05};

  DistEdgeList<int64_t> * DEL = new DistEdgeList<int64_t>();

  //Use the edgelist supplied to build the data structure
  DEL->GenGraphData(edgeList, vertexCount);

  MPI_Barrier(MPI_COMM_WORLD);

  PSpMat_Int64 * G = new PSpMat_Int64(*DEL, false); 
  delete DEL;	// free memory
  MPI_Barrier(MPI_COMM_WORLD);

  Aeff =  PSpMat_Bool(*G);			// Convert to Boolean
  delete G;

  Aeff.OptimizeForGraph500(optbuf);		// Should be called before threading is activated
  MPI_Barrier(MPI_COMM_WORLD);

  //Perform BFS iteration now
  FullyDistVec<int64_t, int64_t> parents ( Aeff.getcommgrid(), Aeff.getncol(), (int64_t) -1);	// identity is -1

  FullyDistSpVec<int64_t, int64_t> fringe(Aeff.getcommgrid(), Aeff.getncol());	// numerical values are stored 0-based

  MPI_Barrier(MPI_COMM_WORLD);

  //Start BFS from vertex 0
  fringe.SetElement(0, 0);

  while(fringe.getnnz() > 0)
  {
    fringe.setNumToInd();
    fringe = SpMV(Aeff, fringe,optbuf);	// SpMV with sparse vector (with indexisvalue flag preset), optimization enabled
    fringe = EWiseMult(fringe, parents, true, (int64_t) -1);	// clean-up vertices that already has parents 
    parents.Set(fringe);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //Report the component size
  FullyDistSpVec<int64_t, int64_t> parentsp = parents.Find(bind2nd(greater<int64_t>(), -1));
  parentsp.Apply(myset<int64_t>(1));


  ostringstream outnew;
  outnew << "Number of vertices found: " << parentsp.Reduce(plus<int64_t>(), (int64_t) 0) << endl; 
  SpParHelper::Print(outnew.str());
}

#endif
