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
 * @file    undirectedChainGen.hpp
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Generates distributed chain graph in parallel
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef UNDIRECTED_CHAIN_GEN_HPP
#define UNDIRECTED_CHAIN_GEN_HPP

//Includes
#include <mpi.h>
#include <iostream>
#include <vector>

//External includes
#include "mxx/timer.hpp"
#include "mxx/partition.hpp"

namespace conn 
{
  namespace graphGen
  {
    /**
     * @class     conn::graphGen::UndirectedChainGen
     * @brief     Parallel graph generator for undirected chains
     */
    class UndirectedChainGen
    {
      public:

        static const uint8_t LOWTOHIGH_IDS = 0;

        //TODO: Low Priority task, enable chain construction with random ids
        //static const uint8_t RANDOM_IDS = 1;
        
        /**
         * @brief                 populates the edge list vector 
         * @param[in] chainLength length of the graph i.e. the count of nodes (not the edges)
         * @param[in] mode        choose between LOWTOHIGH_IDS, RANDOM_IDS
         *                        like 0-1-2...100 (low to high ids)
         *                        or 45-72-21..34 (random ids)
         * @param[out] edgeList   input vector to fill up
         */
        template <typename T>
        void populateEdgeList( std::vector< std::pair<T, T> > &edgeList, 
            uint64_t chainLength, 
            uint8_t mode,
            const mxx::comm &comm = mxx::comm())
        {
          mxx::section_timer timer;

          //sanity check
          assert(mode == LOWTOHIGH_IDS);

          //Know the portion of graph to be generated locally 
          //Lets divide the chainLength value by the number of processes
          mxx::partition::block_decomposition<T> part(chainLength, comm.size(), comm.rank());

          if(mode == LOWTOHIGH_IDS)
          {
            //First node (vertex) id on this MPI node
            T beginNodeId = part.excl_prefix_size(); 
            T lastNodeId = part.excl_prefix_size() + part.local_size() -1 ;

            if(part.local_size() > 0)
            {
              for(int i = beginNodeId; i < lastNodeId; i++)
              {
                edgeList.emplace_back(i, i + 1);
                edgeList.emplace_back(i + 1, i);
              }

              //If not last rank, attach edges from last node of this rank to first node of next rank
              if(part.prefix_size() < chainLength)
              {
                edgeList.emplace_back(lastNodeId, lastNodeId+1);
                edgeList.emplace_back(lastNodeId+1, lastNodeId);
              }
            }
          }

          timer.end_section("graph generation completed");
        }

    };
  }
}

#endif
