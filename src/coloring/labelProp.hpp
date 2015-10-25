/**
 * @file    labelProp.hpp
 * @ingroup group
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Connected component labeling using label propagation (or coloring) approach
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef LABEL_PROPAGATION_HPP 
#define LABEL_PROPAGATION_HPP

//Includes
#include <iostream>

//Own includes
#include "coloring/tupleComp.hpp"
#include "coloring/labelProp_utils.hpp"
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
     * @class                     conn::coloring::ccl
     * @brief                     supports parallel connected component labeling using label propagation technique
     * @tparam[in]  pIdtype       type used for partition ids
     * @tparam[in]  nIdType       type used for node id
     * @tparam[in]  OPTIMIZATION  optimization level for benchmarking, use loadbalanced for the best version 
     */
    template<typename pIdtype = uint32_t, typename nIdType = uint64_t, uint8_t OPTIMIZATION = opt_level::loadbalanced>
    class ccl 
    {
      public:
        //Type for saving parition ids
        using partitionIdtype = pIdtype;

        //Type for saving node ids
        using nodeIdType = nIdType;

        //This is the communicator which participates for computing the components
        mxx::comm comm;

      private:

        using T = std::tuple<pIdtype, pIdtype, nodeIdType>;
        std::vector<T> tupleVector;

        //Used during initialization of <Pn>
        //Also used to mark partitions as stable
        pIdtype MAX = std::numeric_limits<pIdtype>::max();

        //Used to mark tuples as stable 
        //partitions would become stable if all its tuples are stable
        pIdtype MAX2 = std::numeric_limits<pIdtype>::max() -1 ;

      public:
        /**
         * @brief                 public constructor
         * @param[in] edgeList    distributed vector of edges
         * @param[in] c           mpi communicator for the execution 
         */
        template <typename edgeListPairsType>
        ccl(edgeListPairsType &edgeList, const mxx::comm &c) : comm(c.copy()) 
        {
          //Parse the edgeList
          convertEdgeListforCCL(edgeList);

          //Re-distribute the tuples uniformly across the ranks
          mxx::distribute_inplace(tupleVector, comm);
        }

        //Compute the connected component labels
        void compute()
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
          if(!mxx::is_sorted(tupleVector.begin(), tupleVector.end(), TpleComp<cclTupleIds::Pc>(), comm))
            mxx::sort(tupleVector.begin(), tupleVector.end(), TpleComp<cclTupleIds::Pc>(), comm);

          //Count unique Pc values
          return mxx::uniqueCount(tupleVector.begin(), tupleVector.end(),  TpleComp<cclTupleIds::Pc>(), comm);
        }

        /**
         * @brief     Free the communicator
         * @note      Its mandatory to call this function 
         */
        void free_comm()
        {
          comm.~comm();
        }

      private:

        /**
         * @brief     converts the edgelist to vector of tuples needed for ccl
         * @details   For the bucket in the edgeList ...<(u, v1), (u,v2)>...
         *            we append <(u, ~, u), (u, ~, v1), (u, ~, v2)> to our tupleVector
         *            We ignore the bucket splits across ranks here, because that shouldn't affect the correctness and complexity
         */
        template <typename edgeListPairsType>
        void convertEdgeListforCCL(edgeListPairsType &edgeList)
        {
          mxx::section_timer timer(std::cerr, comm);

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
          auto totalTupleCount = mxx::reduce(tupleVector.size(), 0, comm);

          LOG_IF(comm.rank() == 0, INFO) << "Total tuple count is " << totalTupleCount;
        }

        /**
         * @brief     run the iterative algorithm for ccl
         */
        void runConnectedComponentLabeling()
        {
          bool converged = false;

          int iterCount = 0;

          mxx::section_timer timer(std::cerr, comm);

          auto begin = tupleVector.begin();
          auto end = tupleVector.end();
          auto mid = begin;   //range mid-end represents active partitions

          while(!converged)
          {
            updatePn(mid, end);

            converged = updatePc(mid, end);

            //parition the dataset into stable and active paritions, if optimization is enabled
            if(!converged && (OPTIMIZATION == opt_level::stable_partition_removed || OPTIMIZATION == opt_level::loadbalanced))
            {
              //use std::partition to move stable tuples to the left
              //After this step, mid-end would represent active tuples
              mid = partitionStableTuples( mid, end );

              if(OPTIMIZATION == opt_level::loadbalanced)
                mid = mxx::block_decompose_partitions(begin, mid, end, comm);
                //Re distributed the tuples to balance the load across the ranks
            }


            iterCount ++;
          }

          timer.end_section("Coloring done");

          LOG_IF(comm.rank() == 0, INFO) << "Iteration count " << iterCount;
        }

        /**
         * @brief             update the Pn layer by sorting the tuples using node ids
         * @param[in] begin   To iterate over the vector of tuples, marks the range of active tuples 
         * @param[in] end     end iterator
         */
        template <typename Iterator>
        void updatePn(Iterator begin, Iterator end)
        {
          //Sort by nid,Pc
          mxx::sort(begin, end, TpleComp2Layers<cclTupleIds::nId, cclTupleIds::Pc>(), comm); 

          //Resolve last and first bucket's boundary splits
          
          //First, find the element with max node id and min Pc locally 
          //Or in other words, get the min Pc of the last bucket
          auto minPcOfLastBucket = mxx::local_reduce(begin, end, TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>());

          //Second, do exscan, look for max nodeid and min Pc on previous ranks
          auto prevMinPc = mxx::exscan(minPcOfLastBucket, TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>(), comm);  

          //We also need to know max Pc of the first bucket on the next rank (to check for stability)
          auto maxPcOfFirstBucket = mxx::local_reduce(begin, end, TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::less, std::greater>());

          //reverse exscan, look for min nodeid and max Pc on forward ranks
          auto nextMaxPc = mxx::exscan(maxPcOfFirstBucket,  TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::less, std::greater>(), comm.reverse()); 

          //Now we can update the Pn layer of all the buckets locally
          for(auto it = begin; it !=  end;)
          {
            //Range of tuples with the same node id
            auto equalRange = conn::utils::findRange(it, end, *it, TpleComp<cclTupleIds::nId>());

            //Range would include atleast 1 element
            assert(std::distance(equalRange.first, equalRange.second) >= 1);

            //Minimum Pc from local bucket
            auto thisBucketsMinPcLocal = mxx::local_reduce(equalRange.first, equalRange.second, TpleReduce<cclTupleIds::Pc>());

            //Maximum Pc from local bucket
            auto thisBucketsMaxPcLocal = mxx::local_reduce(equalRange.first, equalRange.second, TpleReduce<cclTupleIds::Pc, std::greater>());

            //For now, mark global minimum as local
            auto thisBucketsMaxPcGlobal = thisBucketsMaxPcLocal;
            auto thisBucketsMinPcGlobal = thisBucketsMinPcLocal;

            //Treat first, last buckets as special cases
            if(equalRange.first == begin)
            {
              //Use value from previous rank
              thisBucketsMinPcGlobal =  comm.rank() == 0 ? thisBucketsMinPcLocal : TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>() (prevMinPc, thisBucketsMinPcLocal);
            }

            if(equalRange.second == end)
            {
              //Use value from next rank
              thisBucketsMaxPcGlobal = comm.rank() == comm.size() - 1 ? thisBucketsMaxPcLocal : TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::less, std::greater>() (nextMaxPc, thisBucketsMaxPcLocal);

            }

            //If min Pc < max Pc for this bucket, update Pn or else mark them as stable
            if(TpleComp<cclTupleIds::Pc>()(thisBucketsMinPcGlobal, thisBucketsMaxPcGlobal))
              std::for_each(equalRange.first, equalRange.second, [&](T &e){
                  std::get<cclTupleIds::Pn>(e) = std::get<cclTupleIds::Pc>(thisBucketsMinPcGlobal);
                  });
            else
              std::for_each(equalRange.first, equalRange.second, [&](T &e){
                  std::get<cclTupleIds::Pn>(e) = MAX2;
                  });

            //Advance the loop pointer
            it = equalRange.second;
          }
        }

        /**
         * @brief     update the Pc layer by choosing min Pn
         * @param[in] begin   To iterate over the vector of tuples, marks the range of active tuples 
         * @param[in] end     end iterator
         * @return    bool value, true if the algorithm is converged
         */
        template <typename Iterator>
          bool updatePc(Iterator begin, Iterator end)
          {
            //converged yet
            uint8_t converged = 1;    // 1 means true, we will update it below

            //Sort by Pc, Pn
            mxx::sort(begin, end, TpleComp2Layers<cclTupleIds::Pc, cclTupleIds::Pn>(), comm); 

            //Resolve last bucket's boundary split
            
            //First, find the element with max Pc and min Pn locally 
            //Or in other words, get the min Pn of the last bucket
            auto minPnOfLastBucket = mxx::local_reduce(begin, end, TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>());

            //Result of exscan, again look for max Pc and min Pn on previous ranks
            auto prevMinPn = mxx::exscan(minPnOfLastBucket, TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>(), comm);  

            //Now we can update the Pc layer of all the buckets locally
            for(auto it = begin; it !=  end;)
            {
              //Range of tuples with the same Pc
              auto equalRange = conn::utils::findRange(it, end, *it, TpleComp<cclTupleIds::Pc>());

              //Range would include atleast 1 element
              assert(std::distance(equalRange.first, equalRange.second) >= 1);


              //Minimum Pn from local bucket
              auto thisBucketsMinPnLocal = mxx::local_reduce(equalRange.first, equalRange.second, TpleReduce<cclTupleIds::Pn>());

              //For now, mark global minimum as local
              auto thisBucketsMinPnGlobal = thisBucketsMinPnLocal;

              //Treat first, last buckets as special cases
              if(equalRange.first == begin)
              {
                //Use value from previous rank
                thisBucketsMinPnGlobal =  comm.rank() == 0 ? thisBucketsMinPnLocal : TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>() (prevMinPn, thisBucketsMinPnLocal);
              }

              //If min Pn < MAX2 for this bucket, update the Pc to new value or else mark the partition as stable
              if(std::get<cclTupleIds::Pn>(thisBucketsMinPnGlobal) < MAX2) {
                converged = 0;
                std::for_each(equalRange.first, equalRange.second, [&](T &e){
                    std::get<cclTupleIds::Pc>(e) = std::get<cclTupleIds::Pn>(thisBucketsMinPnGlobal);
                    });
              }
              else
                std::for_each(equalRange.first, equalRange.second, [&](T &e){
                    std::get<cclTupleIds::Pn>(e) = MAX;
                    });

              //Advance the loop pointer
              it = equalRange.second;
            }

            //Know convergence of all the ranks
            uint8_t allConverged;
            mxx::allreduce(&converged, 1, &allConverged, mxx::min<uint8_t>(), comm);

            return (allConverged == 1  ? true : false);
          }

        template <typename Iterator>
          Iterator partitionStableTuples(Iterator begin, Iterator end)
          {
            return std::partition(begin, end, [&](const T &e){
                return std::get<cclTupleIds::Pn>(e) == MAX;
                });
          }
    };

  }
}

#endif
