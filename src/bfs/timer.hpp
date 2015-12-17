/**
 * @file    timer.hpp
 * @ingroup bfs
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Timer log switch used in the bfsRunner
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef BFS_TIMER_HPP 
#define BFS_TIMER_HPP

//external includes
#include "mxx/timer.hpp"

//Switch to 1 if verbose time log is required during the 
//label propagation, else keep it 0
#define BFS_ENABLE_TIMER 0

namespace conn 
{
  namespace bfs
  {

#if COLORING_ENABLE_TIMER
    using Timer = mxx::section_timer_impl<std::chrono::duration<double, std::milli> >;
#else
    using Timer = mxx::empty_section_timer_impl;
#endif

    //Time moment 
    typedef typename std::chrono::steady_clock::time_point timePoint;

    //Clock 
    typedef typename std::chrono::steady_clock clock;

    //Duration type in milliseconds
    typedef typename std::chrono::duration<double, std::milli> duration;

  }
}

#endif
