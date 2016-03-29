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
 * @file    utils_exportGraphDotFormat.cpp
 * @ingroup 
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Writes the edgeList into dot files for running graphviz 
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
#include "graphGen/common/dotExport.hpp"

//External includes
#include "extutils/logging.hpp"
#include "extutils/argvparser.hpp"


INITIALIZE_EASYLOGGINGPP
using namespace std;
using namespace CommandLineProcessing;

int main(int argc, char** argv)
{
  // Initialize the MPI library:
  MPI_Init(&argc, &argv);

  //Initialize the communicator
  mxx::comm comm;

  /**
   * COMMAND LINE ARGUMENTS
   */

  LOG_IF(!comm.rank(), INFO) << "This executable exports graph into dot files.";

  //Parse command line arguments
  ArgvParser cmd;

  cmd.setIntroductoryDescription("This executable exports graph into dot files");
  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption("input", "dbg or kronecker or generic", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  cmd.defineOption("file", "input file (if input = dbg or generic)", ArgvParser::OptionRequiresValue);
  cmd.defineOption("scale", "scale of the graph (if input = kronecker)", ArgvParser::OptionRequiresValue);
  cmd.defineOption("outputPath", "path to the directory where files will be written", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

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

  LOG_IF(!comm.rank(), INFO) << "Generating graph";

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


  std::string outputPath =  cmd.optionValue("outputPath");

  conn::graphGen::writeEdgeListDotFormat(edgeList, outputPath, comm);

  comm.barrier();

  LOG_IF(!comm.rank(), INFO) << "Files written to folder " << outputPath;
  LOG_IF(!comm.rank(), INFO) << "You can join these files in the rank order for using as input to graphViz";

  MPI_Finalize();
  return(0);
}
