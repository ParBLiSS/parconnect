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
#include "debruijn/de_bruijn_construct_engine.hpp"
#include "debruijn/de_bruijn_nodes_distributed.hpp"

namespace conn 
{
  namespace graphGen
  {
    /**
     * @class                     conn::graphGen::deBruijnGraph
     * @brief                     Builds the edgelist of de Bruijn graph 
     * @details                   Sequences are expected in the FASTQ format
     *                            Restrict the alphabets of DNA to {A,C,G,T} 
     */
    class deBruijnGraph
    {
      public:

        //Kmer size set to 31, and alphabets set to 4 nucleotides
        using Alphabet = bliss::common::DNA;
        using KmerType = bliss::common::Kmer<31, Alphabet>;


        //BLISS internal data structure for storing de bruijn graph
        template <typename EdgeEnc>
          using NodeMapType = typename bliss::de_bruijn::de_bruijn_nodes_distributed<
          KmerType, bliss::de_bruijn::node::edge_exists<EdgeEnc>, int,
          bliss::kmer::transform::lex_less,
          bliss::kmer::hash::farm>;


        //Parser type, depends on the sequence file format
        //We restrict the usage to FASTQ format
        template <typename baseIter>
          using SeqParser = typename bliss::io::FASTQParser<baseIter>;


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
          bliss::de_bruijn::de_bruijn_engine<NodeMapType> idx(comm);

          //Build the de Bruijn graph as distributed map
          idx.template build<SeqParser>(fileName, comm);

          auto it = idx.cbegin();

          //Deriving data type of de Bruijn graph storage container
          using mapPairType = typename std::iterator_traits<decltype(it)>::value_type;
          using constkmerType =  typename std::tuple_element<0, mapPairType>::type; 
          using kmerType = typename std::remove_const<constkmerType>::type; //Remove const from nodetype
          using edgeCountInfoType = typename std::tuple_element<1, mapPairType>::type;

          //Temporary storage for each kmer's neighbors in the graph
          std::vector<kmerType> tmpNeighborVector1;
          std::vector<kmerType> tmpNeighborVector2;

          bliss::kmer::transform::lex_less<KmerType> minKmer;

          static_assert(std::is_same<typename kmerType::KmerWordType, uint64_t>::value, "Kmer word type should be set to uint64_t");

          //Read the index and populate the edges inside edgeList
          for(; it != idx.cend(); it++)
          {
            auto sourceKmer = it->first;

            //Get incoming neighbors
            bliss::de_bruijn::node::node_utils<kmerType, edgeCountInfoType>::get_in_neighbors(sourceKmer, it->second, tmpNeighborVector1);

            //Get outgoing neigbors
            bliss::de_bruijn::node::node_utils<kmerType, edgeCountInfoType>::get_out_neighbors(sourceKmer, it->second, tmpNeighborVector2);

            //typename kmerType::KmerWordType* sourceVertexData = minKmer(sourceKmer).getData();
            auto s = minKmer(sourceKmer).getData()[0];

            //Push the edges to our edgeList
            for(auto &e : tmpNeighborVector1)
            {

              auto d = minKmer(e).getData()[0];
              edgeList.emplace_back(s, d);
            }

            //Same procedure for the outgoing edges
            for(auto &e : tmpNeighborVector2)
            {
              auto d = minKmer(e).getData()[0];
              edgeList.emplace_back(s, d);
            }
          }

          timer.end_section("graph generation completed");
        }

    };
  }
}

#endif
