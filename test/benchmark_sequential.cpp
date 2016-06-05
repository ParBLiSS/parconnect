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
 * @file    benchmark_sequential.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes connected components in a graph given as a general file input using sequential REM algorithm 
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

//Includes
#include <mpi.h>
#include <iostream>

//Own includes
#include "graphGen/fileIO/graphReader.hpp"
#include "graphGen/deBruijn/deBruijnGraphGen.hpp"
#include "graphGen/graph500/graph500Gen.hpp"
#include "graphGen/undirectedChain/undirectedChainGen.hpp" 
#include "graphGen/common/reduceIds.hpp"

//External includes
#include "extutils/logging.hpp"
#include "extutils/argvparser.hpp"
#include "mxx/reduction.hpp"
#include "mxx/utils.hpp"
#include "mxx/timer.hpp"

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

  LOG_IF(!comm.rank(), INFO) << "Starting sequential implementation to compute connectivity";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Benchmark for sequential implementation to compute connectivity");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("input", "dbg or kronecker or generic or chain", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("file", "input file (if input = dbg or generic)", ArgvParser::OptionRequiresValue);
  cmd.defineOption("scale", "scale of the graph (if input = kronecker)", ArgvParser::OptionRequiresValue);
  cmd.defineOption("chainLength", "length of undirected chain graph (if input = chain)", ArgvParser::OptionRequiresValue);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) std::cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  if (comm.size() > 1)
  {
    if (!comm.rank()) std::cout << "Run sequential benchmark using single process only"  << "\n";
    exit(1);
  }

  /**
   * GENERATE GRAPH
   */
  using vertexIdType = int64_t;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  std::size_t nVertices;

  LOG_IF(!comm.rank(), INFO) << "Generating graph";

#ifdef BENCHMARK_CONN
  mxx::section_timer timer(std::cerr, comm);
#endif

  //Construct graph based on the given input mode
  if(cmd.optionValue("input") == "generic")
  {
    //Input file
    std::string fileName;
    if(cmd.foundOption("file"))
      fileName = cmd.optionValue("file");
    else
    {
      std::cout << "Required option missing: '--file'\n";
      exit(1);
    }

    LOG_IF(!comm.rank(), INFO) << "Input file -> " << fileName;

    //Add reverse of the edges
    bool addReverse = true;

    //Object of the graph generator class
    conn::graphGen::GraphFileParser<char *, vertexIdType> g(edgeList, addReverse, fileName, comm);

    //Populate the edgeList
    g.populateEdgeList();
  }
  else if(cmd.optionValue("input") == "dbg")
  {
    //Input file
    std::string fileName;
    if(cmd.foundOption("file"))
      fileName = cmd.optionValue("file");
    else
    {
      std::cout << "Required option missing: '--file'\n";
      exit(1);
    }

    LOG_IF(!comm.rank(), INFO) << "Input file -> " << fileName;

    //Object of the graph generator class
    conn::graphGen::deBruijnGraph g;

    //Populate the edgeList
    g.populateEdgeList(edgeList, fileName, comm); 
  }
  else if(cmd.optionValue("input") == "kronecker")
  {
    int scale;
    if(cmd.foundOption("scale"))
      scale = std::stoi(cmd.optionValue("scale"));
    else
    {
      std::cout << "Required option missing: '--scale'\n";
      exit(1);
    }

    LOG_IF(!comm.rank(), INFO) << "Scale -> " << scale;

    int edgefactor = 16;

    //Object of the graph500 generator class
    conn::graphGen::Graph500Gen g;

    //Populate the edgeList
    g.populateEdgeList(edgeList, scale, edgefactor, comm); 
  }
  else if(cmd.optionValue("input") == "chain")
  {
    int chainLength;
    if(cmd.foundOption("chainLength"))
      chainLength = std::stoi(cmd.optionValue("chainLength"));
    else
    {
      std::cout << "Required option missing: '--chainLength'\n";
      exit(1);
    }

    LOG_IF(!comm.rank(), INFO) << "Chain length -> " << chainLength;

    //Object of the chain generator class
    conn::graphGen::UndirectedChainGen g;

    //Populate the edgeList 
    g.populateEdgeList(edgeList, chainLength, comm);
  }
  else
  {
    std::cout << "Wrong input value given" << std::endl;
    exit(1);
  }

  //Count of edges in the graph
  std::size_t nEdges = conn::graphGen::globalSizeOfVector(edgeList, comm);

#ifdef BENCHMARK_CONN
  timer.end_section("Graph construction completed");
#endif

  /**
   * COMPUTE CONNECTIVITY
   */

  auto start = std::chrono::steady_clock::now();

  LOG_IF(!comm.rank(), INFO) << "Beginning computation, benchmark timer started";

  //Relable the ids
  conn::graphGen::permuteVectorIds(edgeList);
  LOG_IF(!comm.rank(), INFO) << "Vertex ids permuted";

  conn::graphGen::reduceVertexIds(edgeList, nVertices, comm);
  LOG_IF(!comm.rank(), INFO) << "Ids compacted for REM algorithm";

  LOG_IF(!comm.rank(), INFO) << "Graph size : vertices -> " << nVertices << ", edges -> " << nEdges/2  << " (x2)";

  /**
   * Execute REM algorithm
   */

  std::size_t num_comp = nVertices;

  std::vector<vertexIdType> p(nVertices);
  for(size_t i = 0; i < nVertices; i++)
    p[i] = i;


  for(auto &e : edgeList)
  {
    vertexIdType rx = e.first;
    vertexIdType ry = e.second;

    while (p[rx] != p[ry]) {   // Check if rx and ry have the same parent
      if (p[rx] < p[ry]) {     // Find the one with the smaller parent
        if (rx == p[rx]) {     // If x is a root we can link
          p[rx] = p[ry];
          num_comp--;
          break;
        }
        vertexIdType z = p[rx];         // Splicing
        p[rx] = p[ry];
        rx = z;
      }
      else {
        if (ry == p[ry]) {      // If y is a root we can link
          p[ry] = p[rx];
          num_comp--;
          break;
        }
        vertexIdType z = p[ry];         // Splicing
        p[ry] = p[rx];
        ry = z;
      }
    } //while loop 
  } //End of loop over edges

  LOG_IF(!comm.rank(), INFO) << "Count of components -> " << num_comp;

  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count(); 

  LOG_IF(!comm.rank(), INFO) << "Time excluding graph construction (ms) -> " << elapsed_time;

  MPI_Finalize();
  return(0);
}
