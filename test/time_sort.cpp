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
 * @file    time_sort.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Measures sorting time of a random vector 
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//External includes
#include "extutils/logging.hpp"
#include "extutils/argvparser.hpp"
#include "mxx/collective.hpp"
#include "mxx/timer.hpp"
#include "mxx/utils.hpp"
#include "mxx/sort.hpp"

INITIALIZE_EASYLOGGINGPP
using namespace std;
using namespace CommandLineProcessing;

int main(int argc, char** argv)
{
  // Initialize the MPI library:
  MPI_Init(&argc, &argv);

  //Initialize the communicator
  mxx::comm comm;

  //Print mpi rank distribution
  mxx::print_node_distribution();

  /**
   * COMMAND LINE ARGUMENTS
   */

  LOG_IF(!comm.rank(), INFO) << "Computing sorting benchmark timings";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Computes sorting benchmark timings");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("count", "total count of elements to sort ", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) std::cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  std::istringstream istr(cmd.optionValue("count"));
  std::size_t count;
  istr >> count;

  LOG_IF(!comm.rank(), INFO) << "mxx sort of " << count << " elements";

  //Adjust and make count equal on all ranks
  count = count - count % comm.size();

  std::size_t localcount = count/comm.size();

  using valueType = int64_t;

  //Declare a edgeList vector to save edges
  std::vector< valueType > buffer(localcount);

  LOG_IF(!comm.rank(), INFO) << "Generating vector";

  //Generating the vector here
  std::generate(buffer.begin(), buffer.end(), std::rand);

  LOG_IF(!comm.rank(), INFO) << "Begin sorting";
  mxx::section_timer timer(std::cerr, comm);

  //Do sorting here
  mxx::sort(buffer.begin(), buffer.end(), comm);

  timer.end_section("sorting completed");

  MPI_Finalize();
  return(0);
}
