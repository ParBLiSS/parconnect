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
 * @file    benchmark_parconnect_auto.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Computes connected components in a graph given as a general file input 
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
#include "graphGen/common/reduceIds.hpp"
#include "coloring/labelProp.hpp"
#include "bfs/bfsRunner.hpp"
#include "dynamic/degreeDistInfo.hpp"

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

  LOG_IF(!comm.rank(), INFO) << "Starting executable for benchmarking in the Student Cluster Competition";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("Benchmark for computing connectivity in the Student Cluster Competition");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("input", "dbg or kronecker or generic", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("file", "input file (if input = dbg or generic)", ArgvParser::OptionRequiresValue);
  cmd.defineOption("scale", "scale of the graph (if input = kronecker)", ArgvParser::OptionRequiresValue);
  cmd.defineOption("pointerDouble", "set to y/n to control pointer doubling during coloring", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError)
  {
    if (!comm.rank()) std::cout << cmd.parseErrorDescription(result) << "\n";
    exit(1);
  }

  /**
   * GENERATE GRAPH
   */
  using vertexIdType = int64_t;

  //Declare a edgeList vector to save edges
  std::vector< std::pair<vertexIdType, vertexIdType> > edgeList;

  //Initialize the distributed vector for saving unique vertices
  std::vector<vertexIdType> uniqueVertexList;

  LOG_IF(!comm.rank(), INFO) << "Generating graph";

  //Fetch the pointer doubling choice
  bool pointerDouble;
  if(cmd.optionValue("pointerDouble") == "y")
    pointerDouble = true;
  else
    pointerDouble = false;


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
  else
  {
    std::cout << "Wrong input value given" << std::endl;
    exit(1);
  }

#ifdef BENCHMARK_CONN
  timer.end_section("Graph construction completed");
#endif

  /**
   * COMPUTE CONNECTIVITY
   */

  comm.barrier();
  auto start = std::chrono::steady_clock::now();

  LOG_IF(!comm.rank(), INFO) << "Beginning computation, benchmark timer started";

  //Relable the ids
  conn::graphGen::permuteVectorIds(edgeList);
  LOG_IF(!comm.rank(), INFO) << "Vertex ids permuted";

#ifdef BENCHMARK_CONN
  timer.end_section("Vertex Ids permuted");
#endif

  bool runBFS = conn::dynamic::runBFSDecision(edgeList, comm);

#ifdef BENCHMARK_CONN
    timer.end_section("Graph fit stastistics calculated");
#endif

  //Call the graph reducer function
  if(runBFS) 
  {
    conn::graphGen::reduceVertexIds(edgeList, uniqueVertexList, comm);
    LOG_IF(!comm.rank(), INFO) << "Ids compacted for BFS run";

#ifdef BENCHMARK_CONN
    timer.end_section("Vertex Ids relabeled (contiguous)");
#endif
  }


  //Count of vertices, edges in the reduced graph
  std::size_t nVertices = conn::graphGen::globalSizeOfVector(uniqueVertexList, comm);
  std::size_t nEdges = conn::graphGen::globalSizeOfVector(edgeList, comm);

  if(runBFS) 
    LOG_IF(!comm.rank(), INFO) << "Graph size : vertices -> " << nVertices << ", edges -> " << nEdges/2;

  if(!runBFS)
    LOG_IF(!comm.rank(), INFO) << "Graph size : edges -> " << nEdges/2;

  //For saving the size of component discovered using BFS
  std::vector<std::size_t> componentCountsResult;

  std::size_t noBFSIterationsExecuted = 0;

  if(runBFS)
  {
    conn::bfs::bfsSupport<vertexIdType> bfsInstance(edgeList, nVertices, comm);

    //Run BFS once
    noBFSIterationsExecuted = bfsInstance.runBFSIterations(1, componentCountsResult); 

#ifdef BENCHMARK_CONN
    timer.end_section("BFS iterations executed");
#endif

    LOG_IF(!comm.rank(), INFO) << "Number of vertices visited by 1st BFS iteration -> " << componentCountsResult[0];

    //Get the remaining edgeList
    bfsInstance.filterEdgeList();

#ifdef BENCHMARK_CONN
    timer.end_section("Remaining graph filtered out");
#endif
  }

  std::size_t countComponents = noBFSIterationsExecuted;

  LOG_IF(!comm.rank(), INFO) << noBFSIterationsExecuted << " BFS iterations executed";

  if(pointerDouble)
  {
    comm.with_subset(edgeList.size() > 0, [&](const mxx::comm& comm){
        conn::coloring::ccl<vertexIdType, conn::coloring::lever::ON> cclInstance(edgeList, comm);
        cclInstance.compute();

        countComponents += cclInstance.computeComponentCount();
    });
  }
  else
  {
     comm.with_subset(edgeList.size() > 0, [&](const mxx::comm& comm){
        conn::coloring::ccl<vertexIdType, conn::coloring::lever::OFF> cclInstance(edgeList, comm);
        cclInstance.compute();

        countComponents += cclInstance.computeComponentCount();
    });
  }

#ifdef BENCHMARK_CONN
    timer.end_section("Coloring completed");
#endif

  countComponents = mxx::allreduce(countComponents, mxx::max<std::size_t>());
  LOG_IF(!comm.rank(), INFO) << "Count of components -> " << countComponents;

  comm.barrier();
  auto end = std::chrono::steady_clock::now();
  auto elapsed_time  = std::chrono::duration<double, std::milli>(end - start).count(); 

  LOG_IF(!comm.rank(), INFO) << "Time excluding graph construction (ms) -> " << elapsed_time;

  MPI_Finalize();
  return(0);
}
