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
 * @file    sort.hpp
 * @ingroup mxx_extra
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   augments extra functions to mxx 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef MXX2_SORT_HPP
#define MXX2_SORT_HPP

//External includes
#include "mxx/comm_fwd.hpp"
#include "mxx/samplesort.hpp"

namespace mxx {

  /**
   * @brief           count unique elements in the global range
   * @param[in] cmp   comparator which was used for sorting the values 
   *                  Example: std::less if integers are sorted in the ascending order
   */
  template <typename Iterator, typename BinaryPredicate>
    std::size_t uniqueCount(Iterator begin, Iterator end, BinaryPredicate cmp, const mxx::comm& comm = mxx::comm()) {
      typedef typename std::iterator_traits<Iterator>::value_type T;

      //Initialize the count to 0
      size_t uniqueCountLocal = 0;

      mxx::comm c = comm.split(begin != end);
      comm.with_subset(begin != end, [&](const mxx::comm& c) {
          size_t n = std::distance(begin, end);

          // send last item to next processor
          T last = *(begin + (n-1));
          T prev = mxx::right_shift(last, c);

          // skip elements which are equal to the last one on the previous processor
          if (c.rank() > 0)
            while (begin != end && !cmp(prev, *begin))
              ++begin;

          if (begin == end)
            return;
          else
            uniqueCountLocal = 1;   //First element we have is unique

          Iterator beginPrev = begin;

          //count unique elements (first entry would be marked as duplicate)
          while (++begin != end)
            if (cmp(*beginPrev, *begin))
            {
              uniqueCountLocal++;
              beginPrev = begin;
            }
      });

      size_t uniqueCountGlobal;
      mxx::allreduce(&uniqueCountLocal, 1, &uniqueCountGlobal ,comm);
      return uniqueCountGlobal;
    }

  template <typename Iterator>
    std::size_t uniqueCount(Iterator begin, Iterator end, const mxx::comm& comm = mxx::comm()) {
      return uniqueCount(begin, end, std::equal_to<typename std::iterator_traits<Iterator>::value_type>(), comm);
    }

#include "mxx/comm_def.hpp"

} // namespace mxx

#endif // MXX_SORT_HPP
