/**
 * @file    timer.hpp
 * @ingroup coloring
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Timer log switch for label propagation method
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef LABEL_PROPAGATION_TIMER_HPP 
#define LABEL_PROPAGATION_TIMER_HPP

//external includes
#include "mxx/timer.hpp"

//Switch to 1 if verbose time log is required during the 
//label propagation, else keep it 0
#define COLORING_ENABLE_TIMER 0

namespace conn 
{
  namespace coloring
  {

#if COLORING_ENABLE_TIMER
    using Timer = mxx::section_timer_impl<std::chrono::duration<double, std::milli> >;
#else
    using Timer = mxx::empty_section_timer_impl;
#endif

  }
}

#endif
