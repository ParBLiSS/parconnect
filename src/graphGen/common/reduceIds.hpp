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
     */
    template <typename E>
      void reduceVertexIds(std::vector<std::pair<E,E>> &edgeList, std::size_t &uniqueVertexCount, const mxx::comm &comm)
      {
        const int SRC = 0, DEST = 1;

        //Ensure the block decomposition of edgeList
        mxx::distribute_inplace(edgeList, comm);

        //Update the DEST layer of all the edges
        {
          //Globally sort all the edges by DEST layer
          if(!mxx::is_sorted(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<DEST>(), comm))
            mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<DEST>(), comm);

          E localVertexIndex = 0;

          //my last vertex id
          E lastLocalVertexId = std::get<DEST>(edgeList.back());

          //my first vertex id
          E firstLocalVertexId = std::get<DEST>(edgeList.front());

          //Start traversal over vertex array
          for(auto it = edgeList.begin(); it != edgeList.end();)
          {
            //Edges whose dest element are equal
            auto edgeListRange = conn::utils::findRange(it, edgeList.end(), *it, conn::utils::TpleComp<DEST>()); 

            /*
             * TODO
             * Preserve these mappings in a separate container if 
             * we wish to revert the vertex ids back
             */

            std::for_each(edgeListRange.first, edgeListRange.second, [&](std::pair<E,E> &e){
                std::get<DEST>(e) = localVertexIndex;
                });

            //Adavance the loop variables
            it = edgeListRange.second;
            localVertexIndex++;
          }

          E localUniqueVertexCount = localVertexIndex;

          //Need to compute prefix count of non-shared vertex buckets
          E exScanUniqueVertices = 0;

          E nextProcsFirstVertexId = mxx::left_shift(firstLocalVertexId, comm);
          if(comm.rank() != (comm.size() - 1) && lastLocalVertexId == nextProcsFirstVertexId)
            localUniqueVertexCount--; //Deduct the extra bucket this rank counted

          exScanUniqueVertices = mxx::exscan(localUniqueVertexCount, std::plus<E>(), comm);
          uniqueVertexCount = mxx::allreduce(localUniqueVertexCount, std::plus<E>(), comm);

          if(!comm.rank()) exScanUniqueVertices = 0;

          //Revise all the vertex ids to global index
          std::for_each(edgeList.begin(), edgeList.end(), [&](std::pair<E,E> &e){
              std::get<DEST>(e) += exScanUniqueVertices;
              });
        } //DEST layer updated

        //Update the SRC layer of all the edges
        {
          mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<SRC>(), comm);

          E localVertexIndex = 0;

          //my last vertex id
          E lastLocalVertexId = std::get<SRC>(edgeList.back());

          //my first vertex id
          E firstLocalVertexId = std::get<SRC>(edgeList.front());

          //Start traversal over vertex array
          for(auto it = edgeList.begin(); it != edgeList.end();)
          {
            //Edges whose dest element are equal
            auto edgeListRange = conn::utils::findRange(it, edgeList.end(), *it, conn::utils::TpleComp<SRC>()); 

            std::for_each(edgeListRange.first, edgeListRange.second, [&](std::pair<E,E> &e){
                std::get<SRC>(e) = localVertexIndex;
                });

            //Adavance the loop variables
            it = edgeListRange.second;
            localVertexIndex++;
          }

          E localUniqueVertexCount = localVertexIndex;

          E exScanUniqueVertices = 0;

          E nextProcsFirstVertexId = mxx::left_shift(firstLocalVertexId, comm);
          if(comm.rank() != (comm.size() - 1) && lastLocalVertexId == nextProcsFirstVertexId)
            localUniqueVertexCount--; //Deduct the extra bucket this rank counted

          exScanUniqueVertices = mxx::exscan(localUniqueVertexCount, std::plus<E>(), comm);
          if(!comm.rank()) exScanUniqueVertices = 0;

          //Revise all the vertex ids to global index
          std::for_each(edgeList.begin(), edgeList.end(), [&](std::pair<E,E> &e){
              std::get<SRC>(e) += exScanUniqueVertices;
              });
        }//SRC layer updated
      }
  }
}
 
#endif
