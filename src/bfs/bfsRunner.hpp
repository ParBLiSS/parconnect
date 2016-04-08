/*
 * Copyright 2016 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
#include "graphGen/common/reduceIds.hpp"
#include "bfs/timer.hpp"
#include "utils/commonfuncs.hpp"

//External includes
#include "extutils/logging.hpp"
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
          static_assert(std::numeric_limits<vertexIdType>::is_signed, "vertexIdType should be a signed type");
          typedef vertexIdType E;

          //Distributed set of unvisited vertices
          std::unordered_set<E> unVisitedVertices;

          //Reference to the distributed edge list 
          std::vector< std::pair<E, E> > &edgeList;

          //Matrix type, to store the adjacency matrix (bool values)
          //from combBLAS implementation 
          using booleanMatrixType = SpParMat <E, bool, SpDCCols<E ,bool> >;

          //Matrix type, to store the adjacency matrix (int values)
          //from combBLAS implementation 
          using integerMatrixType = SpParMat <E, E, SpDCCols<E, E> >;

          //Optimization buffer (used as a parameter in combBLAS function calls)
          //TODO: Generalize these types
          OptBuf<int32_t, int64_t> optbuf;

          //Degrees of each vertex (useful while computing MTEPS score)
          FullyDistVec<E, E> degrees;

          //Record MTEPS score of each iteration
          std::vector<double> MTEPS;

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
        bfsSupport(std::vector< std::pair<E, E> > &edgeList, std::size_t vertexCount,
                  const mxx::comm &comm) : edgeList(edgeList), comm(comm.copy())
        {
          //List of edges, distributed in 1D fashion
          DistEdgeList<E> *DEL = new DistEdgeList<E>();

          //Copy our edgeList to CombBLAS format of edgeList
          DEL->GenGraphData(edgeList, vertexCount);

          comm.barrier();

          integerMatrixType *G = new integerMatrixType(*DEL, false); 
          delete DEL;	// free memory

          comm.barrier();

          //Compute the vertex degrees
          G->Reduce(degrees, Row, plus<E>(), static_cast<E>(0));	// Identity is 0 

          comm.barrier();

          //Now represent the adj matrix in the boolean format
          A =  booleanMatrixType(*G);			// Convert to Boolean
          delete G;

          //Copied the statement from TopDownBFS code
          //Some kind of optimization for graph500 graphs is expected here
          A.OptimizeForGraph500(optbuf);		          

          comm.barrier();

          //Helper to initialize the unvisited vertices buffer
          FullyDistVec<E,E> tmp(A.getcommgrid(), A.getncol(), (E)-1);

          //Record the local array size
          localDistVecSize = tmp.LocArrSize();

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
         * @return                            number of iterations executed by BFS
         */
        std::size_t runBFSIterations(std::size_t noIterations, std::vector<std::size_t> &countComponentSizes)
        {
          //Execute BFS noIterations times
          for(int i = 0; i < noIterations; i++) 
          {
            //Parent array (acts as a list of vertices in a component for us)
            FullyDistVec<E, E> parents(A.getcommgrid(), A.getncol(), (E) -1);	// numerical values are stored 0-based

            //Exscan of vertex count kept on previous ranks
            E offsetForLocalToGlobal = mxx::exscan(localDistVecSize, comm);

            //Get the source vertex
            E srcPoint = getSource(offsetForLocalToGlobal);

            //If all vertices are visited, then exit
            if(srcPoint == MAX)
            {
              LOG_IF(comm.rank() == 0, INFO) << "All vertices already covered, no more BFS iterations required";
              return i;
            }

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

            timePoint t1 = clock::now(); 

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

            comm.barrier();

            FullyDistSpVec<E, E> parentsp = parents.Find(std::bind2nd(std::greater<E>(), -1));
            parentsp.Apply(myset<E>(1));

            //Number of edges traversed
            //std::size_t nEdgesTraversed = EWiseMult(parentsp, degrees, false, 0).Reduce(plus<E>(), 0);
            E nEdgesTraversed = EWiseMult(parentsp, degrees, false, (E) 0).Reduce(plus<E>(), (E) 0);

            //Record the end time of this BFS iteration
            timePoint t2 = clock::now(); 

            //MTEPS = Million edges traversed per second, time is computed in milli seconds
            double MTEPS_Score = static_cast<double>(nEdgesTraversed) / duration(t2 - t1).count() /1000000000.0;

            //Pushing the MTEPS score to our MTEPS vector
            //Take the minimum, although they won't vary much due to barriers
            MTEPS.push_back(mxx::allreduce(MTEPS_Score, mxx::min<double>(), comm));

            comm.barrier();
          }

          return noIterations;

        }

        /**
         * @brief                             Remove the edges corresponding to vertices which have been 
         *                                    covered by BFS
         * @details                           Find the splitters from the sorted edgelist and bucket the
         *                                    unvisited vertices
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

          //Now each rank contains the list of vertices that were not visited during BFS

          //New edgelist  
          std::vector< std::pair<E,E> > edgeListNew;

          //Push the unexplored edges to edgeListNew
          {
            const int SRC = 0, DEST = 1;

            //Globally sort all the edges by SRC layer
            //Should be quick as edgeList was latest sorted by SRC while reducing ids
            if(!mxx::is_sorted(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<SRC>(), comm))
              mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<SRC>(), comm);

            //Define the splitters using the first SRC element of the edge
            auto allSplitters = mxx::allgather(std::get<SRC>(edgeList.front()));
            allSplitters.erase(allSplitters.begin());

            //Initialize functor that assigns rank to all the unique vertices 
            conn::graphGen::vertexToBucketAssignment<E> vertexRankAssigner(allSplitters);

            //All2all to perform the bucketing
            mxx::all2all_func(unVisitedVerticesArray, vertexRankAssigner, comm);

            //Local sort
            std::sort(unVisitedVerticesArray.begin(), unVisitedVerticesArray.end());

            //Do the filtering among the subset of ranks which have non-zero unvisited elements
            //This is required in the cases where BFS traverses all or almost all the edges
            comm.with_subset(unVisitedVerticesArray.size() > 0, [&](const mxx::comm& comm){

                //To resolve the boundary splits across edgelist, we need to do a left shift of first element
                //of this array from ranks 1.. p-1
                E nextProcsFirstUnvisitedVertex = mxx::left_shift(unVisitedVerticesArray.front(), comm);

                //ranks 0.. p-2
                if(comm.rank() < comm.size() - 1)
                  unVisitedVerticesArray.push_back(nextProcsFirstUnvisitedVertex);

                //Start traversal over vertex array
                auto it2 = edgeList.begin();
                for(auto it = unVisitedVerticesArray.begin(); it != unVisitedVerticesArray.end(); it++)
                {
                  //Edges whose SRC element equals the unvisited vertex element
                  auto edgeListRange = conn::utils::findRange(it2, edgeList.end(), *it, conn::utils::TpleComp<SRC>()); 

                  //Insert these edges to our new edgeList
                  edgeListNew.insert(edgeListNew.end(), edgeListRange.first, edgeListRange.second);

                  it2 = edgeListRange.second;
                }
            }); //End of lambda function

          }

          //Replace the content of edgeList with the unvisited edges
          edgeList.assign(edgeListNew.begin(), edgeListNew.end());

            comm.with_subset(edgeList.size() > 0, [&](const mxx::comm& comm){
            //Ensure the block decomposition of edgeList
            mxx::distribute_inplace(edgeList, comm);
          });

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
