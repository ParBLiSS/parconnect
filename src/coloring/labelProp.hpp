/**
 * @file    labelProp.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Custom tuple comparators
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef LABEL_PROPAGATION_HPP 
#define LABEL_PROPAGATION_HPP

//Includes
#include <iostream>

//Own includes
#include "coloring/tupleComp.hpp"
#include "utils/commonfuncs.hpp"
#include "utils/logging.hpp"
#include "utils/prettyprint.hpp"

#include "mxx/sort.hpp"
#include "mxx_extra/sort.hpp"
#include "mxx/timer.hpp"
#include "mxx/comm.hpp"

namespace conn 
{
  namespace coloring
  {

    /**
     * @brief     enum to declare tuple format used during coloring
     */
    enum cclTupleIds
    {
      Pc,         //Pc layer, current partition id
      Pn,         //Pn layer, candidate parition ids
      nId      //Node id , remains unchanged
    };

    /**
     * @brief     pair format of each edge in the edgeList
     */
    enum edgeListTIds
    {
      src,        //source vertex id for edge
      dst         //destination vertex id for edge
    };

    /**
     * @class     conn::coloring::ccl
     * @brief     supports parallel connected component labeling using label propagation technique
     */
    template<typename pIdtype = uint32_t, typename nIdType = uint64_t>
    class ccl 
    {
      public:
        //Type for saving parition ids
        using partitionIdtype = pIdtype;

        //Type for saving node ids
        using nodeIdType = nIdType;

      private:

        using T = std::tuple<pIdtype, pIdtype, nodeIdType>;
        std::vector<T> tupleVector;

        pIdtype MAX = std::numeric_limits<pIdtype>::max();

      public:
        //Constructor
        template <typename edgeListPairsType>
        ccl(edgeListPairsType &edgeList) 
        {
          //Parse the edgeList
          convertEdgeListforCCL(edgeList);

          //Re-distribute the tuples uniformly across the ranks
          mxx::distribute_inplace(tupleVector, mxx::comm());
        }

        //Compute the connected component labels
        void compute(mxx::comm comm = mxx::comm())
        {
          //Size of vector should be >= 0
          assert(tupleVector.begin() != tupleVector.end());

          runConnectedComponentLabeling();
        }

        /**
         * @brief     count the components in the graph after ccl (useful for debugging/testing)
         * @return    count 
         * @note      should be called after computing connected components. 
         *            assumes vector is sorted by Pc
         */
        std::size_t getComponentCount()
        {
          //Vector should be sorted by Pc
          assert(is_sorted(tupleVector.begin(), tupleVector.end(), TpleComp<cclTupleIds::Pc>()));

          //Count unique Pc values
          return mxx::uniqueCount(tupleVector.begin(), tupleVector.end(),  TpleComp<cclTupleIds::Pc>());
        }

      private:

        /**
         * @brief     converts the edgelist to vector of tuples needed for ccl
         * @details   For the bucket in the edgeList ...<(u, v1), (u,v2)>...
         *            we append <(u, ~, u), (u, ~, v1), (u, ~, v2)> to our tupleVector
         *            We ignore the bucket splits across ranks here, because that shouldn't affect the correctness and complexity
         */
        template <typename edgeListPairsType>
        void convertEdgeListforCCL(edgeListPairsType &edgeList, mxx::comm comm = mxx::comm())
        {
          mxx::section_timer timer;

          //Sort the edgeList by src id of each edge
          mxx::sort(edgeList.begin(), edgeList.end(), TpleComp<edgeListTIds::src>(), comm); 

          //Reserve the approximate required space in our vector
          tupleVector.reserve(edgeList.size());

          for(auto it = edgeList.begin(); it != edgeList.end(); )
          {
            //Range of edges with same source vertex
            auto equalRange = conn::utils::findRange(it, edgeList.end(), *it, TpleComp<edgeListTIds::src>());

            //Range would include atleast 1 element
            assert(std::distance(equalRange.first, equalRange.second) >= 1);

            //Insert the self loop
            tupleVector.emplace_back(std::get<edgeListTIds::src>(*it), MAX, std::get<edgeListTIds::src>(*it)); 

            //Insert other vertex members in this partition 
            for(auto it2 = equalRange.first; it2 != equalRange.second; it2++)
              tupleVector.emplace_back(std::get<edgeListTIds::src>(*it2), MAX, std::get<edgeListTIds::dst>(*it2));;

            it = equalRange.second;
          }

          timer.end_section("vector of tuples initialized for ccl");

          //Log the total count of tuples 
          auto totalTupleCount = mxx::reduce(tupleVector.size(), 0);

          LOG_IF(comm.rank() == 0, INFO) << "Total tuple count is " << totalTupleCount;
        }

        /**
         * @brief     run the iterative algorithm for ccl
         */
        void runConnectedComponentLabeling(mxx::comm comm = mxx::comm())
        {
          bool converged = false;

          int iterCount = 0;

          mxx::section_timer timer;

          while(!converged)
          {
            updatePn();
            converged = updatePc();
            iterCount ++;
          }

          timer.end_section("Coloring done");

          LOG_IF(comm.rank() == 0, INFO) << "Iteration count " << iterCount;
        }

        /**
         * @brief     update the Pn layer by sorting the tuples using node ids
         */
        void updatePn(mxx::comm comm = mxx::comm())
        {
          //Sort by nid
          mxx::sort(tupleVector.begin(), tupleVector.end(), TpleComp<cclTupleIds::nId>(), comm); 

          //We need to run a exscan over first element of last bucket to resolve boundary splits
          //Find the element with max node id and min Pc locally 
          auto firstElementOfLastBucket = mxx::local_reduce(tupleVector.begin(), tupleVector.end(), TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>());

          //Result of exscan, again look for max nodeid and min Pc on previous ranks
          auto prevMinPc = mxx::exscan(firstElementOfLastBucket, TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>(), comm);  

          //Now we can update the Pn layer of all the buckets locally
          for(auto it = tupleVector.begin(); it !=  tupleVector.end();)
          {
            //Range of tuples with the same node id
            auto equalRange = conn::utils::findRange(it, tupleVector.end(), *it, TpleComp<cclTupleIds::nId>());

            //Range would include atleast 1 element
            assert(std::distance(equalRange.first, equalRange.second) >= 1);
            
            //Minimum from local bucket
            auto thisBucketsMinPc = mxx::local_reduce(equalRange.first, equalRange.second, TpleReduce<cclTupleIds::Pc>());

            //Compare against the exscanned value (choose exscan value if it has equal node id and lesser pn)
            //This check will be meaningful only for the first local bucket
            auto getLowerTupleValue = TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>() (prevMinPc, thisBucketsMinPc);

            //Update this bucket
            std::for_each(equalRange.first, equalRange.second, [&getLowerTupleValue](T &e){
                std::get<cclTupleIds::Pn>(e) = std::get<cclTupleIds::Pc>(getLowerTupleValue);
                });

            //Advance the loop pointer
            it = equalRange.second;
          }
        }

        /**
         * @brief     update the Pc layer by choosing min Pn
         * @return    bool value, true if the algorithm is converged
         */
        bool updatePc(mxx::comm comm = mxx::comm())
        {
          //converged yet
          uint8_t converged = 1;    // 1 means true, we will update it below

          //Sort by Pc
          mxx::sort(tupleVector.begin(), tupleVector.end(), TpleComp<cclTupleIds::Pc>(), comm); 

          //We need to run a exscan over first element of last bucket to resolve boundary splits
          //Find the element with max Pc and min Pn locally 
          auto firstElementOfLastBucket = mxx::local_reduce(tupleVector.begin(), tupleVector.end(), TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>());

          //Result of exscan, again look for max Pc and min Pn on previous ranks
          auto prevMinPc = mxx::exscan(firstElementOfLastBucket, TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>(), comm);  

          //Now we can update the Pn layer of all the buckets locally
          for(auto it = tupleVector.begin(); it !=  tupleVector.end();)
          {
            //Range of tuples with the same Pc
            auto equalRange = conn::utils::findRange(it, tupleVector.end(), *it, TpleComp<cclTupleIds::Pc>());

            //Range would include atleast 1 element
            assert(std::distance(equalRange.first, equalRange.second) >= 1);
            
            //Minimum from local bucket
            auto thisBucketsMinPc = mxx::local_reduce(equalRange.first, equalRange.second, TpleReduce<cclTupleIds::Pn>());

            //Compare against the exscanned value (choose exscan value if it has equal node id and lesser pn)
            //This check will be meaningful only for the first local bucket
            auto getLowerTupleValue = TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>() (prevMinPc, thisBucketsMinPc);

            //Update this bucket (merge the partition Pc with new min Pn)
            std::for_each(equalRange.first, equalRange.second, [&getLowerTupleValue, &converged](T &e){

                //Check for termination
                auto oldPc = std::get<cclTupleIds::Pc>(e);
                auto newPc = std::get<cclTupleIds::Pn>(getLowerTupleValue);
                if(oldPc != newPc) converged = 0;

                //Update the Pc value
                std::get<cclTupleIds::Pc>(e) = std::get<cclTupleIds::Pn>(getLowerTupleValue);
                });

            //Advance the loop pointer
            it = equalRange.second;
          }

          //Know convergence of all the ranks
          uint8_t allConverged;
          mxx::allreduce(&converged, 1, &allConverged, mxx::min<uint8_t>());

          return (allConverged == 1  ? true : false);
        }
    };

  }
}

#endif
