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
 * @file    commonfuncs.hpp
 * @ingroup utils
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Common utility functions
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef COMMON_FUNCTIONS_HPP 
#define COMMON_FUNCTIONS_HPP

//Includes
#include <iostream>
#include <fstream>

//External includes
#include "mxx/comm.hpp"


namespace conn 
{
  namespace utils
  {

    //Implements std::equal_range using sequential scan
    template<typename ForwardIterator, typename T, class Compare>
      std::pair<ForwardIterator, ForwardIterator> findRange(ForwardIterator first, ForwardIterator second, const T& val, Compare comp)
      {
        //Default return values
        ForwardIterator it1 = second;
        ForwardIterator it2 = second;

        ForwardIterator start = first;

        for(; start != second ; start++)
        {
          //Check equivalence to val
          if(!comp(*start, val) && !comp(val, *start))
          {
            it1 = start;
            break;
          }
        }

        for(; start != second ; start++)
        {
          //Check where the equivalency is broken
          if(! (!comp(*start, val) && !comp(val, *start)))
          {
            it2 = start;
            break;
          }
        }

        return std::make_pair(it1, it2);
      }

    /**
     * @brief           Functor for comparing tuples by single index layer 
     * @tparam layer    Tuple's index which is used for comparison
     * @tparam op       comparator, default as std::less
     */
    template <size_t layer, template<typename> class op = std::less>
      struct TpleComp
      {
        //Compare two tuples using their values 
        template<typename T>
          bool operator() (T const &t1, T const &t2)
          {
            return op<typename std::tuple_element<layer, T>::type>() (std::get<layer>(t1), std::get<layer>(t2));
          }

        //Compare tuple directly against a value
        template<typename T>
          bool operator() (T const &t1, typename std::tuple_element<layer, T>::type const &val)
          {
            return op<typename std::tuple_element<layer, T>::type>() (std::get<layer>(t1), val);
          }

        template<typename T>
          bool operator() (typename std::tuple_element<layer, T>::type const &val, T const &t2)
          {
            return op<typename std::tuple_element<layer, T>::type>() (val, std::get<layer>(t2));
          }
      };
   
    /**
     * @brief           Functor for comparing tuples by two index layers
     * @tparam layer1   Tuple's primary index which is used for comparison
     * @tparam layer2   Tuple's secondary index which is used for comparison
     * @tparam op1      layer1 comparator, default as std::less
     * @tparam op2      layer2 comparator, default as std::less
     */
    template <size_t layer1, size_t layer2, template<typename> class op1 = std::less, template<typename> class op2 = std::less>
      struct TpleComp2Layers
      {
        template<typename T>
          bool operator() (T const &t1, T const &t2)
          {
            return 
              op1<typename std::tuple_element<layer1, T>::type>() (std::get<layer1>(t1), std::get<layer1>(t2)) ||
              (
               std::get<layer1>(t1) == std::get<layer1>(t2) &&
               op2<typename std::tuple_element<layer2, T>::type>() (std::get<layer2>(t1), std::get<layer2>(t2))
              );
          }
      };

    /**
     * @brief           Reduces two tuples by using comparator on single layer
     * @tparam layer    Tuple's index which is used for comparison
     * @tparam op       comparator, default as std::less
     */
    template <size_t layer, template<typename> class op = std::less>
      struct TpleReduce
      {
        template<typename T>
          T operator() (T const &t1, T const &t2)
          {
            if(TpleComp<layer, op>() (t1, t2))
              return t1;
            else
              return t2;
          }
      };
   
    /**
     * @brief           Reduces two tuples by using comparators on two layers
     * @tparam layer1   Tuple's primary index which is used for comparison
     * @tparam layer2   Tuple's secondary index which is used for comparison
     * @tparam op1      layer1 comparator, default as std::less
     * @tparam op2      layer2 comparator, default as std::less
     */
    template <size_t layer1, size_t layer2, template<typename> class op1 = std::less, template<typename> class op2 = std::less>
      struct TpleReduce2Layers
      {
        template<typename T>
          T operator() (T const &t1, T const &t2)
          {
            if(TpleComp2Layers< layer1, layer2, op1, op2 >() (t1, t2))
              return t1;
            else
              return t2;
          }
      };

    //Write all the edges to a single file
    //Could be slow since it writes sequentially
    template <typename EdgeVectorType>
      void writeEdgesToFile(EdgeVectorType &edges, std::string &outFile, const mxx::comm &comm)
      {
        //Gather complete edgeList on rank 0
        auto fullEdgeList = mxx::gatherv(edges, 0, comm);

        if(!comm.rank())
        {

          std::ofstream f;
          f.open(outFile);

          for(auto &e : fullEdgeList)
            f << e.first << " " << e.second << std::endl;

          f.close();

        }
      }
  }
}

#endif
