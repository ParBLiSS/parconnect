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
 * @file    utils.hpp
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Common utility functions related to edgelist and graph representation.
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef GRAPH_UTILS_HPP
#define GRAPH_UTILS_HPP

//Own inlcudes
#include "utils/commonfuncs.hpp"

//external includes
#include "extutils/logging.hpp"
#include "mxx/reduction.hpp"
#include "mxx/sort.hpp"

namespace conn 
{
  namespace graphGen
  {

    /**
     * @details   Helper function to print the edgelist distribution
     *            across ranks during the algorithm's exection. Prints the min, mean and 
     *            max count of the active tuples
     */
    template <typename Iterator, typename T = std::size_t>
      void printEdgeListDistribution(Iterator begin, Iterator end, mxx::comm &comm)
      {
        T localWorkLoad = std::distance(begin, end);

        T maxLoad  = mxx::reduce(localWorkLoad, 0, mxx::max<T>() , comm);
        T minLoad  = mxx::reduce(localWorkLoad, 0, mxx::min<T>() , comm);
        T meanLoad = mxx::reduce(localWorkLoad, 0, std::plus<T>(), comm)/ comm.size();

        auto sep = ",";
        LOG_IF(comm.rank() == 0, INFO) << "Distribution of edge list; min-mean-max : " << minLoad << sep << meanLoad << sep << maxLoad;
      }

    /**
     * @brief     Helper function to confirm edges are represented in both directions
     * @note      Use while debugging or testing only
     * @details   let E1 be the order of edgeList sorted by <SRC, DEST> and E2 be the 
     *            copy of edgeList sorted by <DEST, SRC>. Then ith edge in E1 should
     *            be equal to flip of ith edge in E2
     * @return    true if the check passes on all the MPI ranks
     */
    template <typename E>
      bool checkEdgeBidirectionality( std::vector< std::pair<E,E> > &edgeList1, mxx::comm &comm)
      {
        const int SRC = 0, DEST = 1;

        assert(std::distance(edgeList1.begin(), edgeList1.end()) > 0);

        //Create a copy
        auto edgeList2 = edgeList1;

        //Sort edgeList1 by SRC, DEST
        mxx::sort(edgeList1.begin(), edgeList1.end(), conn::utils::TpleComp2Layers<SRC, DEST>(), comm);

        //Sort edgeList2 by DEST, SRC
        mxx::sort(edgeList2.begin(), edgeList2.end(), conn::utils::TpleComp2Layers<DEST, SRC>(), comm);
        
        bool localCheck = true; 

        //Check equality 
        for(int i = 0; i < edgeList1.size(); i++)
        {
          localCheck = localCheck && (edgeList1[i].first == edgeList2[i].second);
          localCheck = localCheck && (edgeList1[i].second == edgeList2[i].first);
        }

        //Check if the condition passed on all ranks
        unsigned char localCheckVal = (unsigned char) localCheck;
        unsigned char globalCheckVal = mxx::allreduce(localCheckVal, mxx::min<unsigned char>(), comm); 
        
        if(globalCheckVal == 0)
          return false;
        else
          return true;
      }

  }
}

#endif
