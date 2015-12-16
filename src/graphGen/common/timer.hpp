/**
 * @file    timer.hpp
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Timer log switch for graph generation methods
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef GRAPHGEN_TIMER_HPP 
#define GRAPHGEN_TIMER_HPP

//external includes
#include "mxx/timer.hpp"

//Switch to 1 if verbose time log is required during the 
//graph generation, else keep it 0
#define ENABLE_TIMER 0

namespace conn 
{
  namespace graphGen
  {

#if ENABLE_TIMER
    using Timer = mxx::section_timer_impl<std::chrono::duration<double, std::milli> >;
#else
    using Timer = mxx::empty_section_timer_impl;
#endif

  }
}

#endif
