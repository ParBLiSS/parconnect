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

#ifdef BENCHMARK_CONN
#define COLORING_ENABLE_TIMER 1
#else
#define COLORING_ENABLE_TIMER 0
#endif

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
