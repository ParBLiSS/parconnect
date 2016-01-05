/**
 * @file    graphReader.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Parallel graph reader using BLISS parallel IO support
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef PAR_GRAPH_READER_HPP 
#define PAR_GRAPH_READER_HPP

//Includes
#include <iostream>

//Own includes
#include "graphGen/common/timer.hpp"
#include "utils/logging.hpp"

//External includes
#include "io/file_loader.hpp"
#include "common/base_types.hpp"
#include "mxx/comm.hpp"


namespace conn 
{
  namespace graphGen
  {

    /**
     * @class     conn::graphGen::graphFileParser
     * @brief     Enables parallel reading of the edgelist file
     */
    template<typename Iterator, typename E>
    class GraphFileParser : public bliss::io::BaseFileParser<Iterator >
    {
      private:

        //Type of inherited class
        using baseType = typename bliss::io::BaseFileParser<Iterator >;

        //MPI communicator
        mxx::comm comm;

        //Switch to determine if reverse of each edge should be included as well
        bool addReverseEdge;

        //Reference to the distributed edge list 
        std::vector< std::pair<E,E> > &edgeList;

      public:

        /**
         * @brief                 constructor for this class
         * @param[in] edgeList    Edgelist to build
         * @param[in] comm        mpi communicator
         */
        template <typename vID>
          GraphFileParser(std::vector< std::pair<vID, vID> > &edgeList, bool addReverseEdge,
              const mxx::comm &comm) :  edgeList(edgeList), 
                                        addReverseEdge(addReverseEdge),
                                        comm(comm.copy())
          {
            comm.barrier();
            LOG_IF(!comm.rank(), INFO) << "GraphFileParser constructor called";
          }


        /**
         * @brief                   populates the edge list vector 
         * @param[in]   filename    name of the file to read
         */
        void populateEdgeList(std::string &filename)
        {
          Timer timer;

          //Value type over which Iterator is defined
          typedef typename std::iterator_traits<Iterator>::value_type IteratorValueType;

          //Define file loader type
          typedef bliss::io::FileLoader<IteratorValueType> FileLoaderType;  

          comm.barrier();
          LOG_IF(!comm.rank(), INFO) << "Defining loader";

          //==== create file Loader
          FileLoaderType loader(filename, comm);

          comm.barrier();
          LOG_IF(!comm.rank(), INFO) << "Fetching L1 block";

          //====  now process the file, one L1 block partition per MPI Rank 
          typename FileLoaderType::L1BlockType partition = loader.getNextL1Block();

          LOG_IF(!comm.rank(), INFO) << "Global file range [" << loader.getFileRange().start << "," << loader.getFileRange().end << "]"; 

          comm.barrier();
          LOG_IF(!comm.rank(), INFO) << "Obtaining range";

          //Range of the file partition this MPI rank reads
          auto localFileRange = partition.getRange();

          LOG(INFO) << "Rank " << comm.rank() << "'s local file range [" << localFileRange.start << "," << localFileRange.end << "]"; 

          //Iterator over data in the file
          typename FileLoaderType::L1BlockType::iterator dataIter = partition.begin();

          //Initialize the byte offset counter over the range of this rank
          std::size_t i = localFileRange.start;

          comm.barrier();
          LOG_IF(!comm.rank(), INFO) << "Adjusting begin iterator ";

          adjustInitialIterPosition(dataIter, partition.end(), i);

          //Backup iterators and offset for back-tracking
          typename FileLoaderType::L1BlockType::iterator backUpdataIter;
          std::size_t backupi;

          bool lastEdgeRead;

          comm.barrier();
          LOG_IF(!comm.rank(), INFO) << "Beginning parsing of the file ";

          //Begin parsing the contents 
          while (dataIter != partition.end()) {

            backUpdataIter = dataIter;
            backupi = i;
            lastEdgeRead = readAnEdge(dataIter, partition.end(), i);

          }

          comm.barrier();
          LOG_IF(!comm.rank(), INFO) << "Reading last edge";


          //All but last rank read one more edge
          if(comm.rank() != comm.size() - 1)
          {
            if(lastEdgeRead == true)
            {
              //Read edge from next record
              LOG(INFO) << "Rank " << comm.rank() << " Current offset[t] " << i; 
              this->findNonEOL(dataIter, dataIter + std::min(loader.getFileRange().end - i, localFileRange.end - localFileRange.start), i);
              readAnEdge(dataIter, dataIter + std::min(loader.getFileRange().end - i, localFileRange.end - localFileRange.start), i);
            }
            else
            {
              //Read full record again
              LOG(INFO) << "Rank " << comm.rank() << " Current offset[f] " << backupi; 
              this->findNonEOL(backUpdataIter, backUpdataIter + std::min(loader.getFileRange().end - backupi, localFileRange.end - localFileRange.start) , backupi);
              readAnEdge(backUpdataIter, backUpdataIter + std::min(loader.getFileRange().end - backupi, localFileRange.end - localFileRange.start) , backupi);
            }
          }

          LOG(INFO) << "Rank " << comm.rank() << " read the last edge"; 

          timer.end_section("File IO completed, graph built");
        }

      private:

        template <typename Iter>
          void adjustInitialIterPosition(Iter& curr, Iter end, std::size_t& i)
          {
            //every rank except 0 skips initial partial or full record
            if(comm.rank() != 0)
            {
              this->findEOL(curr, end, i);
              this->findNonEOL(curr, end, i);
            }

            //skip initial sentences beginning with '%'
            while(*curr == '%')
            {
              if(curr == end)
                return;

              this->findEOL(curr, end, i);
              this->findNonEOL(curr, end, i);
            }
          }

        /**
         * @brief             reads an edge assuming iterator points to 
         *                    beginning of a valid record
         * @param[in] curr    current iterator position
         * @return            true if edge is read successfully,
         *                    false if partition boundary is encountered 
         */
        template <typename Iter>
          bool readAnEdge(Iter& curr, Iter end, std::size_t& offset) 
          {
            //make sure we point to non EOL value
            this->findNonEOL(curr, end, offset);

            //string to save the contents of a line
            std::string readLine;

            //till we see end of line
            while( ! (*curr == baseType::eol || *curr == baseType::cr)) 
            {
              //return if we are crossing the boundary
              if(curr == end)
                return false;

              //keep pushing the character to string
              readLine.append(curr, 1);

              //advance iterator
              curr++; offset++;

              //return if we are crossing the boundary
              if(curr == end)
                return false;
            }

            parseStringForEdge(readLine);

            return true;
          }


        /**
         * @brief             assumes string with two integers as input,
         *                    parse the integers and insert to edgeList
         * @param[in] record  string with 2 integers separated by space
         */
        inline void parseStringForEdge(std::string &record)
        {
          std::stringstream stream(record);

          size_t n = std::count(record.begin(), record.end(), ' ');

          E vertex1, vertex2;

          stream >> vertex1;
          stream >> vertex2;

          if(n == 1)
          {
            edgeList.emplace_back(vertex1, vertex2);
            if(addReverseEdge)
              edgeList.emplace_back(vertex2, vertex1);
          }
        }
    };
  }
}

#endif
