/**
 * @file    graph500Gen.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Wrapper around code files for parallel graph500 graph generator.
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef GRAPH500_GEN_HPP
#define GRAPH500_GEN_HPP

//Includes
#include <mpi.h>
#include <iostream>
#include <vector>

//Own includes
#include "graphGen/graph500/include/make_graph.h"
#include "mxx/timer.hpp"

namespace conn 
{
  namespace graphGen
  {
    /**
     * @class     conn::graphGen::graph500Gen
     * @brief     Wraps around the c code for graph500 generation
     */
    class graph500Gen
    {
      public:
        //edge id is restricted to type int64_t in the c code, so we use the same 
        using T = int64_t; 
        
        static const uint8_t DIRECTED = 0;
        static const uint8_t UNDIRECTED = 1;

      private:
        static const double initiator[4];

      public:

        /**
         * @brief                 populates the edge list vector 
         * @param[out] edgeList   input vector to fill up
         * @param[in] mode        choose between DIRECTED or UNDIRECTED
         * @param[in] scale       scale of the graph
         * @param[in] edgeFactor  edgeFactor of the graph
         */
        void populateEdgeList( std::vector< std::pair<int64_t, int64_t> > &edgeList, 
            uint8_t scale, 
            uint8_t edgeFactor, 
            uint8_t mode,
            const mxx::comm &comm = mxx::comm())
        {
          //seeds to use
          int64_t seeds[2] = {1,2};

          int64_t desired_nedges = edgeFactor * (1UL << scale);

          int64_t nedges;

          T *edges;

          mxx::section_timer timer;

          //Use the internal function to populate the edges
          make_graph(scale, desired_nedges, seeds[0], seeds[1], initiator, &nedges, &edges);

          //sanity check
          assert(mode == DIRECTED || mode == UNDIRECTED);

          for (int i = 0; i < nedges; ++i) 
          {
            T src = edges[2*i];
            T dest = edges[2*i+1];

            if ((src >= 0) && (dest >= 0)) 
            { 
              // valid edge
              edgeList.emplace_back(src, dest);

              //Insert the reverse edge if the mode is undirected
              if(mode == UNDIRECTED)
                edgeList.emplace_back(dest, src);
            }
          }

          //Free temporary memory
          free(edges);

          timer.end_section("graph generation completed");
        }

    };

    const double graph500Gen::initiator[4] =  {.57, .19, .19, .05};
  }
}

#endif
