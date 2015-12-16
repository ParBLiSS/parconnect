/**
 * @file    reduceIds.hpp
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Given a graph as list of edges, it reduces the vertex ids so that they are contiguous
 *          This process reduces memory usage if graph is represented as an adjacency matrix later.
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef GRAPH_REDUCE_HPP
#define GRAPH_REDUCE_HPP

//Includes
#include <mpi.h>
#include <iostream>
#include <algorithm>

//Own includes
#include "utils/commonfuncs.hpp"

//External includes
#include "mxx/distribution.hpp"
#include "mxx/timer.hpp"
#include "mxx/sort.hpp"
#include "mxx/algos.hpp"
#include "mxx/utils.hpp"

namespace conn 
{
  namespace graphGen
  {

    /**
     * @brief             Functor to assign rank for an edge tuple
     * @param[in] layer   0 or 1 depending upon using the 
     *                    source or dest element of edge
     */
    template <typename E,  uint8_t layer>
      struct edgeToBucketAssignment
      {
        std::vector<E> &splitters;
        
        //Constructor 
        edgeToBucketAssignment(std::vector<E> &splitters) :
          splitters(splitters) {
        }

        int operator() (const std::pair<E,E>& edge)
        {
          if(splitters.size() > 0)
          {
            //Extract the vertex id from edge
            auto valueToBucket = std::get<layer>(edge);

            //Fetch the bucket it belongs to
            auto lower = std::lower_bound(splitters.begin(), splitters.end(), valueToBucket);

            //lower bound will return one less value if valueToBucket is itself one of the splitters
            if(lower != splitters.end() && *lower == valueToBucket)
              lower++;

            int value = std::distance(splitters.begin(), lower);
            
            //value should be between [0, p-1]
            assert(0 <= value && value < splitters.size()+1);

            return value;
          }
          else 
          {
            //Splitter size is 0, implying that comm.size() = 1
            return 0;
          }
        }
      };

    /**
     * @brief             Edge comparator
     * @param[in] layer   0 or 1 depending upon using the 
     *                    source or dest layer of edge
     */
    template <uint8_t layer>
      struct edgeComparator
      {
        //Compare 2 edges
        template <typename E>
          bool operator() (const std::pair<E,E>& e1, const std::pair<E,E>& e2){
            return std::get<layer>(e1) < std::get<layer>(e2);
          }

        //Compare edge and vertex
        template <typename E>
          bool operator() (const std::pair<E,E>& e, const E& v1){
            return std::get<layer>(e) < v1;
          }

        template <typename E>
          bool operator() (const E& v1, const std::pair<E,E>& e){
            return v1 < std::get<layer>(e);
          }


      };

    /**
     * @brief                         Given a graph as list of edges, it updates all the vertex ids so
     *                                that they are contiguos from 0 to |V-1|
     * @param[in]  edgeList           distributed vector of edges
     * @param[out] uniqueVertexList   distributed array of original vertex ids (global index of vertex here is its new id)
     * @details                       Implemented using bucketing and all2all communication
     *                                Due to all2all, this function call will not have good scalability
     */
    template <typename E>
      void reduceVertexIds(std::vector<std::pair<E,E>> &edgeList, std::vector<E> &uniqueVertexList, mxx::comm &comm)
      {
        const int SRC = 0, DEST = 1;

        //Ensure the block decomposition of edgeList
        mxx::distribute_inplace(edgeList, comm);

        ///Phase 1: lets create a unique and sorted list of vertices

        //Populate all the vertices local to this rank
        for(auto &e : edgeList)
        {
          uniqueVertexList.push_back(e.first);
          uniqueVertexList.push_back(e.second);
        }

        //Remove duplicates locally
        auto last = std::unique(uniqueVertexList.begin(), uniqueVertexList.end());
        uniqueVertexList.erase(last, uniqueVertexList.end());

        assert(uniqueVertexList.size() > 0);

        //Sort the vertices, remove the duplicates and enable equal distribution
        mxx::sort(uniqueVertexList.begin(), uniqueVertexList.end(), comm);
        last = mxx::unique(uniqueVertexList.begin(), uniqueVertexList.end(), std::equal_to<E>(), comm);
        uniqueVertexList.erase(last, uniqueVertexList.end());

        mxx::distribute_inplace(uniqueVertexList, comm);
        mxx::sort(uniqueVertexList.begin(), uniqueVertexList.end(), comm); //Sorting is required after distribute

        //Local count of vertices on this rank
        auto localCountOfVertices = uniqueVertexList.size();

        assert(localCountOfVertices > 0);

        //Sum count of vertices on previous ranks
        auto exScanCountOfVertices = mxx::exscan(localCountOfVertices, comm);
        if(comm.rank() == 0)
          exScanCountOfVertices = 0;

        ///Phase 2: Find splitters, bucket and update the edgeTuples 
        auto firstLocalVertexId = uniqueVertexList.front();

        auto allSplitters = mxx::allgather(firstLocalVertexId);
        allSplitters.erase(allSplitters.begin()); //Discard element that came from rank 0

        //Update the SRC layer of all the edges
        {
          //Initialize functor that assigns bucket id to each edge by its source vertex
          edgeToBucketAssignment<E, SRC>  edgeRankAssigner(allSplitters);

          mxx::all2all_func(edgeList, edgeRankAssigner, comm);

          //Local sort
          edgeComparator<SRC> cmp;
          std::sort(edgeList.begin(), edgeList.end(), cmp);

          //Start traversal over vertex array
          auto it2 = edgeList.begin();
          for(auto it = uniqueVertexList.begin(); it != uniqueVertexList.end(); it++)
          {
            //Edges whose src element is equal to current uniqueVertexList element
            auto edgeListRange = conn::utils::findRange(it2, edgeList.end(), *it, cmp); 

            //Id that it should be assigned
            auto trueId = exScanCountOfVertices + std::distance(uniqueVertexList.begin(), it);

            std::for_each(edgeListRange.first, edgeListRange.second, [&](std::pair<E,E> &e){
                std::get<SRC>(e) = trueId;
                });
          }
        }

        //Update the DEST layer of all the edges
        {
          //Initialize functor that assigns bucket id to each edge by its dest vertex
          edgeToBucketAssignment<E, DEST>  edgeRankAssigner(allSplitters);

          //Do all to all
          mxx::all2all_func(edgeList, edgeRankAssigner, comm);

          //Local sort
          edgeComparator<DEST> cmp;
          std::sort(edgeList.begin(), edgeList.end(), cmp);

          //Start traversal over vertex array
          auto it2 = edgeList.begin();
          for(auto it = uniqueVertexList.begin(); it != uniqueVertexList.end(); it++)
          {
            //Edges whose dest element is equal to current uniqueVertexList element
            auto edgeListRange = conn::utils::findRange(it2, edgeList.end(), *it, cmp); 

            //Id that it should be assigned
            auto trueId = exScanCountOfVertices + std::distance(uniqueVertexList.begin(), it);

            std::for_each(edgeListRange.first, edgeListRange.second, [&](std::pair<E,E> &e){
                std::get<DEST>(e) = trueId;
                });
          }
        }

        //Ensure the block decomposition of edgeList
        mxx::distribute_inplace(edgeList, comm);
      }

    /**
     * @brief   returns the global size of the vector
     */
    template <typename vectorType>
      std::size_t globalSizeOfVector(vectorType v, mxx::comm &comm)
      {
        std::size_t localSize = v.size();
        return mxx::allreduce(localSize, std::plus<std::size_t>());
      }
  }
}
 
#endif
