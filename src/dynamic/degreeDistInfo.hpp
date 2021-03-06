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
 * @file    degreeDistInfo.hpp   
 * @ingroup dynamic
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Helps in choosing between coloring and BFS by using the degree distribution
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef GRAPH_DEGREE_DIST_HPP
#define GRAPH_DEGREE_DIST_HPP

//Includes
#include <mpi.h>
#include <iostream>
#include <algorithm>

//Own includes
#include "utils/commonfuncs.hpp"

//External includes
#include "mxx/distribution.hpp"
#include "mxx/reduction.hpp"
#include "mxx/sort.hpp"
#include "plfit/src/plfit.h"
#include "mxx/timer.hpp"

namespace conn 
{
  namespace dynamic
  {

    template <typename VecType>
      double fitCurve(VecType &data)
      {
        /* construct the plfit options */
        plfit_discrete_options_t plfit_discrete_options;
        plfit_discrete_options_init(&plfit_discrete_options);
        plfit_discrete_options.finite_size_correction = 0;
        plfit_discrete_options.p_value_method = PLFIT_P_VALUE_SKIP;
        plfit_discrete_options.alpha_method = PLFIT_LBFGS;

        plfit_result_t result;

        plfit_discrete(&data[0], data.size(), &plfit_discrete_options, &result);

        //fprintf(stderr,"\n");
        //fprintf(stderr,"\talpha = %12.5lf\n", result.alpha);
        //fprintf(stderr,"\txmin  = %12.5lf\n", result.xmin );
        //fprintf(stderr,"\tL     = %12.5lf\n", result.L    );
        //fprintf(stderr,"\tD     = %12.5lf\n", result.D    );
        //fprintf(stderr,"\n");

        return result.D;

      }

    /**
     * @brief                   Decides if its optimal to run BFS iteration based on the degree distribution
     * @param[in]  edgeList     distributed vector of edges
     * @return                  true if BFS should be executed, false otherwise
     * @NOTE                    Assumes each edge is present both ways in the edgeList vector         
     */
    template <typename E>
      bool runBFSDecision(std::vector<std::pair<E,E>> &edgeList, mxx::comm &comm)
      {
#ifdef BENCHMARK_CONN
        mxx::section_timer timer(std::cerr, comm);
#endif

        const int SRC = 1, DEST = 0;  //Reverse the layers to avoid sorting during relabeling vertices

        const int sampler = 11;

        //Ensure the block decomposition of edgeList
        mxx::distribute_inplace(edgeList, comm);

        //Sort by source, dest vertex
        mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp2Layers<SRC,DEST>(), comm);

        //Map to hold degree frequency
        std::unordered_map<std::size_t, std::size_t> degreeCountMap;

        std::size_t maxDegree = 0;

        //Vector to hold boundary vertex degrees
        std::vector<std::pair<E,E>> boundaryVertexDegrees;

        for(auto it = edgeList.begin(); it != edgeList.end();)
        {
          auto equalSrcRange = conn::utils::findRange(it, edgeList.end(), *it, conn::utils::TpleComp<SRC>()); 

          //A vertex may have duplicate destination vertices, ignore them
          auto uniqueDestVerticeIter = std::unique(equalSrcRange.first, equalSrcRange.second);

          auto currentDegree = std::distance(equalSrcRange.first, uniqueDestVerticeIter);

          if(equalSrcRange.first == edgeList.begin())   //First bucket
            boundaryVertexDegrees.emplace_back(std::get<SRC>(*it), currentDegree);
          else if(equalSrcRange.second == edgeList.end() && equalSrcRange.first != edgeList.begin())   //Last bucket (and different from first bucket)
            boundaryVertexDegrees.emplace_back(std::get<SRC>(*it), currentDegree);
          else
          {
            //This bucket is completely local to this rank
            degreeCountMap[currentDegree]++;

            if(currentDegree > maxDegree)
              maxDegree = currentDegree;
          }

          it = equalSrcRange.second;
        }

        auto globalBoundaryVertexDegrees = mxx::gatherv(boundaryVertexDegrees, 0, comm);

        if(comm.rank() == 0)
        {
          const int SRC = 0, COUNT = 1;

          for(auto it = boundaryVertexDegrees.begin(); it != boundaryVertexDegrees.end();)
          {
            auto equalSrcRange = conn::utils::findRange(it, globalBoundaryVertexDegrees.end(), *it, conn::utils::TpleComp<SRC>());
          
            auto currentDegree = std::accumulate(equalSrcRange.first, equalSrcRange.second, 0, [&](const E &p1, const std::pair<E,E> &p2){
                return p1 + std::get<COUNT>(p2);
                });

            if(currentDegree > maxDegree)
              maxDegree = currentDegree;

            degreeCountMap[currentDegree]++;

            it = equalSrcRange.second;
          }
        }
        
        maxDegree = mxx::allreduce(maxDegree, mxx::max<std::size_t>(), comm);

        //Convert to vector
        //Holds the frequency of degrees [1,2..maxDegree]
        std::vector<double> degreeHolder(maxDegree);
        for(size_t i = 0; i < degreeHolder.size(); i++)
          degreeHolder[i] = degreeCountMap[i+1];

        //Reduce degree vector to rank 0
        auto globalDegreeCounts = mxx::reduce(degreeHolder, 0, comm);

        //Add 1 to each element (for stable curve fitting)
        std::for_each(globalDegreeCounts.begin(), globalDegreeCounts.end(), [](double& d) { d+=1.0;});

        int decision = 0;   //false

#ifdef BENCHMARK_CONN
        timer.end_section("Degree distribution completed");
#endif

        if(!comm.rank())
        {
          auto statisticVal = fitCurve(globalDegreeCounts);
          if(statisticVal < 0.05)
          {
            LOG_IF(!comm.rank(), INFO) << "Kolmogorov-Smirnov statistic " << statisticVal << " (below 0.05)"; 
            decision = 1; //true
          }
          else
          {
            LOG_IF(!comm.rank(), INFO) << "Kolmogorov-Smirnov statistic " << statisticVal << " (above 0.05)"; 
          }

        }

#ifdef BENCHMARK_CONN
        timer.end_section("K-S test completed");
#endif

        auto gbDecision = mxx::allreduce(decision, mxx::max<int>(), comm);

        return gbDecision == 1;
      }

  }
}
 
#endif
