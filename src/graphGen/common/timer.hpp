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
