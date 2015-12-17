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
#include "utils/logging.hpp"

//external includes
#include "mxx/reduction.hpp"

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

  }
}

#endif
