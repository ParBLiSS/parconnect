/**
 * @file    sort.hpp
 * @brief   augments extra functions to mxx 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef MXX2_SORT_HPP
#define MXX2_SORT_HPP

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
