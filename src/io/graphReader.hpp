/**
 * @file    graphReader.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Parallel graph reader using MPI file IO
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef PAR_GRAPH_READER_HPP 
#define PAR_GRAPH_READER_HPP

//Includes
#include <iostream>

//External includes
#include "io/file_loader.hpp"

namespace conn 
{
  namespace io
  {

    /**
     * @class     conn::graphGen::graphFileParser
     * @brief     Enables parallel reading of the edgelist file
     */
    template<typename Iterator>
    class GraphFileParser
    {

      protected:

        /// constant representing the EOL character
        static const typename std::iterator_traits<Iterator>::value_type eol = '\n';

        /// alias for this type, for use inside this class.
        using type = GraphFileParser<Iterator>;

      public:

        /**
         * @brief                   populates the edge list vector 
         * @param[in]   filename    name of the file to read
         * @param[out]  edgeList    input vector to fill up
         */
        void populateEdgeList( std::vector< std::pair<E, E> > &edgeList, 
            std::string &filename,
            const mxx::comm &comm)
        {
          Timer timer;

          //Define file loader type
          typedef bliss::io::FileLoader<CharType> FileLoaderType;  

          //==== create file Loader
          FileLoaderType loader(comm, filename);

          //====  now process the file, one L1 block partition per MPI Rank 
          typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();

          //Range of the file partition this MPI rank reads
          auto localFileRange = partition.getRange();

          //Iterator over data in the file
          auto dataIter = partition.begin();

          // initialize the byte offset counter
          std::size_t i = localFileRange.start;

          while (i < localFileRange.end) {

            bool EOLfound = findEOL(dataIter, t.end, i);
            std::cout << comm.rank() << " EOL " << EOLfound << "\n";

            //if (!EOLfound)
              //lastBrokenHeader = true;

            ++dataIter;
            ++i;
          }

          timer.end_section("File IO completed, graph built");
        }

      private:

        /**
         *  @brief                    search for first EOL character in a iterator.
         *  @param[in/out]  iter      iterator to advance
         *  @param[in]      end       position to stop the traversal
         *  @param[in/out]  offset    the global offset within the file.
         *  @return                   True if EOL char is found, and false if not found.
         **/
        inline bool findEOL(Iterator& iter, const size_t end, size_t &offset) const{
          while (*iter != type::eol) {
            if (offset == end) { 
              return false;
            }
            ++iter;
            ++offset;
          }
          return true;
        }

    }


  }
}

#endif
