/**
 * @file    deBruijnGraphGen.hpp 
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Builds the edgelist for de bruijn graph using BLISS library.
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef DE_BRUIJN_GEN_HPP
#define DE_BRUIJN_GEN_HPP

//Includes
#include <mpi.h>
#include <iostream>
#include <vector>

//Own includes
#include "graphGen/common/timer.hpp"

//External includes
#include "debruijn/de_bruijn_node_trait.hpp"

namespace conn 
{
  namespace graphGen
  {
    /**
     * @class                     conn::graphGen::deBruijnGraph
     * @brief                     Builds the edgelist of de Bruijn graph 
     * @tparam[in]  NodeMapType   Type of map to use for storing DBG
     * @tparam[in]  SeqParser     Parser type, depends on the sequence file format
     */
    template <typename NodeMapType, template <typename> class SeqParser>
    class deBruijnGraph
    {
      public:

        /** 
         * @brief                 populates the edge list vector 
         * @param[in]   fileName
         * @param[out]  edgelist
         */
        template <typename E>
        void populateEdgeList( std::vector< std::pair<E, E> > &edgeList, 
            std::string &fileName,
            const mxx::comm &comm)
        {
          Timer timer;

          //Initialize the map
          NodeMapType idx(comm, comm.size());

          //Build the de Bruijn graph as distributed map
          idx.template build<SeqParser>(fileName, comm);

          auto it = idx.cbegin();

          //Deriving data type of de Bruijn graph storage container
          using mapPairType = typename std::iterator_traits<decltype(it)>::value_type;
          using constkmerType =  typename std::tuple_element<0, graphNodeType>::type; 
          using kmerType = typename std::remove_const<constNodeInfoType>::type; //Remove const from nodetype
          using edgeCountInfoType = typename std::tuple_element<1, graphNodeType>::type;

          //Temporary storage for each kmer's neighbors in the graph
          std::vector<kmerType> tmpNeighborVector1;
          std::vector<kmerType> tmpNeighborVector2;

          //Read the index and populate the edges inside edgeList
          for(; it != idx.cend(); it++)
          {
            auto sourceKmer = it->first;

            //Get incoming neighbors
            bliss::de_bruijn::node::node_utils<nodeInfoType, edgeCountInfoType>::get_in_neighbors(sourceKmer, it->second, tmpNeighborVector1);

            //Get outgoing neigbors
            bliss::de_bruijn::node::node_utils<nodeInfoType, edgeCountInfoType>::get_out_neighbors(sourceKmer, it->second, tmpNeighborVector2);

            //Push the edges to our edgeList
            for(auto &e : tmpNeighborVector1)
              edgeList.push_back(std::min(sourceKmer, sourceKmer.reverse_complement()), 
                                  std::min(e, e.reverse_complement()));

            for(auto &e : tmpNeighborVector2)
              edgeList.push_back(std::min(sourceKmer, sourceKmer.reverse_complement()), 
                                  std::min(e, e.reverse_complement()));

            //TODO: Is each edge pushed in both directions ? 
          }

          timer.end_section("graph generation completed");
        }

    };
  }
}

#endif
