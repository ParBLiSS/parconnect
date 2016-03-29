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
 * @file    dotExport.hpp
 * @ingroup graphGen
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   takes the edgeList and writes the graph in graphviz dot format
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */


#ifndef DOTFORMAT_EXPORT_HPP
#define DOTFORMAT_EXPORT_HPP

//Includes
#include <iostream>
#include <algorithm>

//External includes
#include "mxx/distribution.hpp"

namespace conn 
{
  namespace graphGen
  {
    /**
     * @brief             writes the edgeList in the dot format
     *                    Each rank writes a file, need to concatenate them in order
     * @NOTE              Assumes each each is represented twice (u-v) plus (v-u) in the edgeList
     * @param outputPath  path to output directory where this function creates the files
     */
    template <typename E>
      void writeEdgeListDotFormat(std::vector<std::pair<E,E>> &edgeList, std::string &outputPath, mxx::comm &comm)
      {
        const int SRC = 0, DEST = 1;

        //Ensure the block decomposition of edgeList
        mxx::distribute_inplace(edgeList, comm);

        std::string fileName = outputPath + "/" + "graph." + std::to_string(comm.rank()) + ".dot";
        std::ofstream out(fileName);

        //Rank 0 writes initial header
        if(comm.rank() == 0)
          out << "graph G {\n";

        for(auto &e : edgeList) {
          if(std::get<SRC>(e) < std::get<DEST>(e)) {
            out << std::to_string(std::get<SRC>(e)) + " -- " + std::to_string(std::get<DEST>(e)) + ";\n";
          }
        }

        //Last ranks writes closing bracket
        if(comm.rank() == comm.size() - 1)
          out << "}\n";

        out.close();
      }


  }
}



#endif
