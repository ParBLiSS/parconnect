/**
 * @file    bfsRunner.hpp
 * @ingroup bfs
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Finds connected components using a BFS iterations
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef BFS_ITERATIONS_HPP 
#define BFS_ITERATIONS_HPP

#include <mpi.h>
#include <iostream>
#include <unordered_set>

//Own includes
#include "utils/logging.hpp"
#include "graphGen/common/reduceIds.hpp"

//External includes
#include "graph500-gen/make_graph.h"
#include "CombBLAS/CombBLAS.h"
#include "mxx/comm.hpp"
#include "mxx/distribution.hpp"


namespace conn 
{
  namespace bfs
  {
    /**
     * @class                     conn::bfs::bfsSupport
     * @brief                     supports parallel connected component labeling using BFS iterations
     * @tparam[in]  vertexIdType  type used for vertices in the distributed edge list
     */
    template <typename vertexIdType>
      class bfsSupport
      {
        private:

          //BFS implementation uses signed integer types
          using E = typename std::make_signed<vertexIdType>::type;

          //Distributed set of unvisited vertices
          std::unordered_set<E> unVisitedVertices;

          //Reference to the distributed edge list 
          std::vector< std::pair<E,E> > &edgeList;

          //Matrix type, to store the adjacency matrix (bool values)
          //from combBLAS implementation 
          using booleanMatrixType = SpParMat <E, bool, SpDCCols<E ,bool> >;

          //Matrix type, to store the adjacency matrix (int values)
          //from combBLAS implementation 
          using integerMatrixType = SpParMat <E, E, SpDCCols<E, E> >;

          //Optimization buffer (used as a parameter in combBLAS function calls)
          //TODO: Generalize these types
          OptBuf<int32_t, int64_t> optbuf;

          //combBLAS distributed storage for adjacency matrix
          booleanMatrixType A;

          //This is the communicator which participates for computing the components
          mxx::comm comm;

          //For convenience, define the maximum value 
          E MAX = std::numeric_limits<E>::max();

          //Size of the parents array local to this rank
          std::size_t localDistVecSize;

        public:

        /**
         * @brief                 constructor, builds the adjacency matrix required for BFS
         * @param[in] edgeList    input graph as distributed edgeList
         * @param[in] vertexCount total count of vertices in the graph i.e. (the highest vertex id + 1)
         *                        assuming vertex id begins from 0
         * @param[in] comm        mpi communicator
         *                        TODO : Enable the communicator restriction on all the BFS functions
         */
        bfsSupport(std::vector< std::pair<E,E> > &edgeList, std::size_t vertexCount,
                  const mxx::comm &comm) : edgeList(edgeList), comm(comm.copy())
        {
          //List of edges, distributed in 1D fashion
          DistEdgeList<E> *DEL = new DistEdgeList<E>();

          //Copy our edgeList to CombBLAS format of edgeList
          DEL->GenGraphData(edgeList, vertexCount);

          integerMatrixType *G = new integerMatrixType(*DEL, false); 
          delete DEL;	// free memory

          //Now represent the adj matrix in the boolean format
          A =  booleanMatrixType(*G);			// Convert to Boolean
          delete G;
          //TODO: Do we have to symmetricize the matrix?

          //Copied the statement from TopDownBFS code
          //Some kind of optimization for graph500 graphs is expected here
          A.OptimizeForGraph500(optbuf);		          

          //Helper to initialize the unvisited vertices buffer
          FullyDistVec<E,E> tmp(A.getcommgrid(), A.getncol(), (E)-1);

          for(E i = 0; i < tmp.LocArrSize(); i++)
          {
            //Note that we are saving local id of every vertex
            //This gets convenient when we erase the visited elements later
            unVisitedVertices.emplace(i);
          }
        }

        /**
         * @brief                             runs multiple bfs iterations for graph connectivity
         * @param[in]   noIterations          upper bound on the count of iterations for BFS runs
         *                                    it can also return if the graph is completely visited
         * @param[out]  countComponentSizes   vector of count of vertices visited during each BFS run
         * @details                           Each bfs run begins from unvisited vertex till it traverses
         *                                    that component
         */
        void runBFSIterations(std::size_t noIterations, std::vector<std::size_t> &countComponentSizes)
        {

          //Parent array (acts as a list of vertices in a component for us)
          FullyDistVec<E, E> parents(A.getcommgrid(), A.getncol(), (E) -1);	// numerical values are stored 0-based

          //Record the parents array size
          localDistVecSize = parents.LocArrSize();

          //Exscan of vertex count kept on previous ranks
          E offsetForLocalToGlobal = mxx::exscan(localDistVecSize, comm);

          for(int i = 0; i < noIterations; i++) 
          {
            //Get the source vertex
            E srcPoint = getSource(offsetForLocalToGlobal);

            //If all vertices are visited, then exit
            if(srcPoint == MAX)
              break;

            //Frontier
            FullyDistSpVec<E, E> fringe(A.getcommgrid(), A.getncol());	// numerical values are stored 0-based

            //Barrier
            MPI_Barrier(comm);

            //Set the source of BFS
            fringe.SetElement(srcPoint, srcPoint);
            parents.SetElement(srcPoint, srcPoint);

            //Remove the source vertex from our vertex set
            fringe.removeFromHash(unVisitedVertices);

            //Set to 1 as we include the source
            std::size_t trackCountOfVerticesVisited = 1;

            //Till the frontier is non-empty
            while (fringe.getnnz() > 0)
            {
              // Top-down
              fringe.setNumToInd();

              //Matrix multiplication
              fringe = SpMV(A, fringe, optbuf);

              //Remove elements from frontier that were already visited before
              fringe = EWiseMult(fringe, parents, true, (int64_t) -1);

              //Update parents array using fringe
              parents.Set(fringe);

              //Remove the newly visited elements from our map of vertices
              fringe.removeFromHash(unVisitedVertices);
              trackCountOfVerticesVisited += fringe.getnnz();
            }

            //Keep record of the number of vertices visited
            countComponentSizes.push_back(trackCountOfVerticesVisited);

            MPI_Barrier(comm);
          }

        }

        /**
         * @brief                             Remove the edges corresponding to vertices which have been 
         *                                    covered by BFS
         * @details                           Find the splitters from the unVisitedVertices and split the edges
         *                                    based on their DEST values 
         *                                    Do all2all and remove the edges through linear scan
         * @note                              This function should be called after running the BFS iterations. 
         */
        void filterEdgeList()
        {
          //Exscan of vertex count kept on previous ranks
          E offsetForLocalToGlobal = mxx::exscan(localDistVecSize, comm);

          //Copy all the unvisited elements from set to vector 
          std::vector<E> unVisitedVerticesArray(unVisitedVertices.begin(), unVisitedVertices.end());
          std::for_each(unVisitedVerticesArray.begin(), unVisitedVerticesArray.end(), [&](E& d){d += offsetForLocalToGlobal;});
          std::sort(unVisitedVerticesArray.begin(), unVisitedVerticesArray.end());

          //Now each rank contains the globally sorted list of vertices that were not visited during BFS

          //Define the splitters using the exscan value obtained above
          auto allSplitters = mxx::allgather(offsetForLocalToGlobal);
          allSplitters.erase(allSplitters.begin());

          //New edgelist  
          std::vector< std::pair<E,E> > edgeListNew;

          //Push the unexplored edges to edgeListNew
          {
            const int SRC = 0, DEST = 1;

            //Initialize functor that assigns bucket id to each edge by its source vertex
            conn::graphGen::edgeToBucketAssignment<E, DEST>  edgeRankAssigner(allSplitters);

            mxx::all2all_func(edgeList, edgeRankAssigner, comm);

            //Local sort
            conn::graphGen::edgeComparator<DEST> cmp;
            std::sort(edgeList.begin(), edgeList.end(), cmp);

            //Start traversal over vertex array
            auto it2 = edgeList.begin();
            for(auto it = unVisitedVerticesArray.begin(); it != unVisitedVerticesArray.end(); it++)
            {
              //Edges whose DEST element equals the unvisited vertex element
              auto edgeListRange = conn::utils::findRange(it2, edgeList.end(), *it, cmp); 

              //There ought to be atleast one edge for this vertex
              assert(std::distance(edgeListRange.first , edgeListRange.second) > 0);

              //Insert these edges to our new edgeList
              edgeListNew.insert(edgeListNew.end(), edgeListRange.first, edgeListRange.second);
            }
          }

          //Replace the content of edgeList with the unvisited edges
          edgeList.assign(edgeListNew.begin(), edgeListNew.end());

          //Ensure the block decomposition of edgeList
          mxx::distribute_inplace(edgeList, comm);
        }

        private:

        /**
         * @brief             returns next source to start the BFS iterations
         * @param[in] offset  its the value addition required to convert local
         *                    ids in unVisitedVertices to global vertex ids
         */
        E getSource(E offset)
        {
          //get candidate from this rank
          E firstLocalElement;

          if(unVisitedVertices.empty())
          {
            //Set to MAX if the map is empty
            firstLocalElement = MAX;
          }
          else
          {
            //Convert local to global index if we found valid candidate
            firstLocalElement = *unVisitedVertices.begin() + offset;
          }

          //Find the minimum candidate among all
          //TODO Should we randomize this selection?
          E source = mxx::allreduce(firstLocalElement, mxx::min<E>(), comm);
          return source;
        }
      };

  }
}

#endif
