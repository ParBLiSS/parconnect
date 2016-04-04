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
 * @file    reduceIds.hpp
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Given a graph as list of edges, it reduces the vertex ids so that they are contiguous
 *          This process reduces memory usage if graph is represented as an adjacency matrix later.
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef GRAPH_REDUCE_HPP
#define GRAPH_REDUCE_HPP

//Includes
#include <mpi.h>
#include <iostream>
#include <algorithm>

//Own includes
#include "utils/commonfuncs.hpp"
#include "graphGen/common/utils.hpp"

//External includes
#include "mxx/distribution.hpp"
#include "mxx/timer.hpp"
#include "mxx/sort.hpp"
#include "mxx/algos.hpp"
#include "mxx/utils.hpp"
#include "hash/invertible_hash.hpp"

namespace conn 
{
  namespace graphGen
  {

    /**
     * @brief             Functor to assign rank for a given vertex id
     */
    template <typename E>
      struct vertexToBucketAssignment
      {
        std::vector<E> &splitters;
        
        //Constructor 
        vertexToBucketAssignment(std::vector<E> &splitters) :
          splitters(splitters) {
        }

        int operator() (const E& valueToBucket)
        {
          if(splitters.size() > 0)
          {
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
     * @brief   returns the global size of the vector
     */
    template <typename vectorType>
      std::size_t globalSizeOfVector(vectorType &v, mxx::comm &comm)
      {
        std::size_t localSize = v.size();
        return mxx::allreduce(localSize, std::plus<std::size_t>());
      }

    /*
     * @brief   relabes vertex ids using invertible hash function
     */
    template <typename E>
      void permuteVectorIds(std::vector<std::pair<E,E>> &edgeList)
      {
        const int SRC = 0, DEST = 1;

        for(auto &e : edgeList)
        {
          hash_64(std::get<SRC>(e));
          hash_64(std::get<DEST>(e));
        }
      }

    /**
     * @brief                         Given a graph as list of edges, it updates all the vertex ids so
     *                                that they are contiguos from 0 to |V-1|
     * @param[in]  edgeList           distributed vector of edges
     * @param[out] uniqueVertexList   distributed array of original vertex ids (global index of vertex here is its new id)
     * @details                       (u, v) edge in the edgeList is transfromed to (x, y) if u is the xth 
     *                                element in the sorted order of all unique vertices, lly for v is yth 
     *                                Implemented using bucketing and all2all communication
     *                                Due to all2all, this function call will not have good scalability
     */
    template <typename E>
      void reduceVertexIds(std::vector<std::pair<E,E>> &edgeList, std::vector<E> &uniqueVertexList, mxx::comm &comm)
      {
        const int SRC = 0, DEST = 1;

        //Ensure the block decomposition of edgeList
        mxx::distribute_inplace(edgeList, comm);

        ///Phase 1: lets create a unique and sorted list of vertices

        //Reserve the space in the beginning itself
        uniqueVertexList.reserve(2*edgeList.size());

        //Populate all the vertices local to this rank
        for(auto &e : edgeList)
        {
          uniqueVertexList.push_back(e.first);
          uniqueVertexList.push_back(e.second);
        }

        //Remove duplicates locally
        auto last = std::unique(uniqueVertexList.begin(), uniqueVertexList.end());
        uniqueVertexList.erase(last, uniqueVertexList.end());

        //Does the graph have enough unique vertices?
        assert(uniqueVertexList.size() > 0);

        //Sort the vertices, remove the duplicates and enable equal distribution
        mxx::sort(uniqueVertexList.begin(), uniqueVertexList.end(), comm);
        last = mxx::unique(uniqueVertexList.begin(), uniqueVertexList.end(), std::equal_to<E>(), comm);
        uniqueVertexList.erase(last, uniqueVertexList.end());

        //Update the DEST layer of all the edges
        {
          //Globally sort all the edges by DEST layer
          mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<DEST>(), comm);

          auto allSplitters = mxx::allgather(std::get<DEST>(edgeList.front()));
          allSplitters.erase(allSplitters.begin()); //Discard element that came from rank 0

          //Initialize functor that assigns rank to all the unique vertices 
          vertexToBucketAssignment<E> vertexRankAssigner(allSplitters);

          //All2all to perform the bucketing
          mxx::all2all_func(uniqueVertexList, vertexRankAssigner, comm);

          //Sum count of vertices on previous ranks
          auto exScanCountOfVertices = mxx::exscan(uniqueVertexList.size(), comm);
          if(comm.rank() == 0)
            exScanCountOfVertices = 0;

          //Local sort
          std::sort(uniqueVertexList.begin(), uniqueVertexList.end());

          //Start traversal over vertex array
          auto it2 = uniqueVertexList.begin();
          for(auto it = edgeList.begin(); it != edgeList.end();)
          {
            //Edges whose dest element are equal
            auto edgeListRange = conn::utils::findRange(it, edgeList.end(), *it, conn::utils::TpleComp<DEST>()); 

            //Iterator in uniqueVertexList pointing to value equal to dest
            auto srcMatch = std::find(it2, uniqueVertexList.end(), std::get<DEST>(*it));

            //Id that it should be assigned
            auto trueId = exScanCountOfVertices + std::distance(uniqueVertexList.begin(), srcMatch);

            std::for_each(edgeListRange.first, edgeListRange.second, [&](std::pair<E,E> &e){
                std::get<DEST>(e) = trueId;
                });

            //Adavance the loop variables
            it2 = srcMatch;
            it = edgeListRange.second;
          }
        } //DEST layer updated


        //Update the SRC layer of all the edges
        {
          //Globally sort all the edges by SRC layer
          mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<SRC>(), comm);

          auto allSplitters = mxx::allgather(std::get<SRC>(edgeList.front()));
          allSplitters.erase(allSplitters.begin()); //Discard element that came from rank 0

          //Initialize functor that assigns rank to all the unique vertices 
          vertexToBucketAssignment<E> vertexRankAssigner(allSplitters);

          //All2all to perform the bucketing
          mxx::all2all_func(uniqueVertexList, vertexRankAssigner, comm);

          //Sum count of vertices on previous ranks
          auto exScanCountOfVertices = mxx::exscan(uniqueVertexList.size(), comm);
          if(comm.rank() == 0)
            exScanCountOfVertices = 0;

          //Local sort
          std::sort(uniqueVertexList.begin(), uniqueVertexList.end());

          //Start traversal over vertex array
          auto it2 = uniqueVertexList.begin();
          for(auto it = edgeList.begin(); it != edgeList.end();)
          {
            //Edges whose src element are equal
            auto edgeListRange = conn::utils::findRange(it, edgeList.end(), *it, conn::utils::TpleComp<SRC>()); 

            //Iterator in uniqueVertexList pointing to value equal to src
            auto srcMatch = std::find(it2, uniqueVertexList.end(), std::get<SRC>(*it));

            //Id that it should be assigned
            auto trueId = exScanCountOfVertices + std::distance(uniqueVertexList.begin(), srcMatch);

            std::for_each(edgeListRange.first, edgeListRange.second, [&](std::pair<E,E> &e){
                std::get<SRC>(e) = trueId;
                });

            //Adavance the loop variables
            it2 = srcMatch;
            it = edgeListRange.second;
          }
        } //SRC layer updated
      }
  }
}
 
#endif
