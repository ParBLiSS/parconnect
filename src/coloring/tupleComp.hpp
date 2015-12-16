/**
 * @file    tupleComp.hpp
 * @ingroup coloring
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Custom tuple comparators
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef TUPLE_COMP_HPP 
#define TUPLE_COMP_HPP

//Includes
#include <iostream>
#include <tuple>

//Own includes

namespace conn 
{
  namespace coloring
  {
    /*
     * COMPARISON
     */

    /**
     * @brief           Functor for comparing tuples by single index layer 
     * @tparam layer    Tuple's index which is used for comparison
     * @tparam op       comparator, default as std::less
     */
    template <size_t layer, template<typename> class op = std::less>
      struct TpleComp
      {
        template<typename T>
          bool operator() (T const &t1, T const &t2)
          {
            return op<typename std::tuple_element<layer, T>::type>() (std::get<layer>(t1), std::get<layer>(t2));
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

    /*
     * REDUCTION
     */
 
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

  }
}

#endif
