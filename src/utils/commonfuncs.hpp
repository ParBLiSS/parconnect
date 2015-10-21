/**
 * @file    commonfuncs.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Common utility functions
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef COMMON_FUNCTIONS_HPP 
#define COMMON_FUNCTIONS_HPP

//Includes
#include <iostream>



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
  }
}

#endif
