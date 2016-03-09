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
 * @file    labelProp.hpp
 * @ingroup coloring
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
#include "coloring/labelProp_utils.hpp"
#include "coloring/timer.hpp" //Timer switch 
#include "utils/commonfuncs.hpp"

//external includes
#include "mxx/sort.hpp"
#include "mxx_extra/sort.hpp"
#include "mxx/comm.hpp"
#include "extutils/logging.hpp"

namespace conn 
{
  namespace coloring
  {

    /**
     * @class                     conn::coloring::ccl
     * @brief                     supports parallel connected component labeling using label propagation technique
     * @tparam[in]  nIdType       type used for node id
     * @tparam[in]  DOUBLING      controls whether pointer doubling would be executed or not, 'ON' by default.
     * @tparam[in]  OPTIMIZATION  optimization level for benchmarking, use loadbalanced for the best version 
     */
    template<typename nIdType = uint64_t, uint8_t DOUBLING = lever::ON, uint8_t OPTIMIZATION = opt_level::loadbalanced>
    class ccl 
    {
      public:
        //Type for saving parition ids
        using pIdtype = nIdType;

        //Type for saving node ids
        using nodeIdType = nIdType;

      private:

        //This is the communicator which participates for computing the components
        mxx::comm comm;

      private:

        using T = std::tuple<pIdtype, pIdtype, nodeIdType>;
        std::vector<T> tupleVector;

        //Used during initialization of <Pn>
        //Also used to mark partitions as stable
        pIdtype MAX_PID = std::numeric_limits<pIdtype>::max();

        //Used to mark tuples as stable 
        //partitions would become stable if all its tuples are stable
        pIdtype MAX_PID2 = std::numeric_limits<pIdtype>::max() -1 ;

        //Used to mark the special tuples used during doubling
        nodeIdType MAX_NID = std::numeric_limits<nodeIdType>::max();

      public:
        /**
         * @brief                 public constructor
         * @param[in] edgeList    distributed vector of edges
         * @param[in] c           mpi communicator for the execution 
         */
        template <typename E>
        ccl(std::vector<std::pair<E,E>> &edgeList, const mxx::comm &c) : comm(c.copy()) 
        {
          //nodeIdType and E should match
          //If they don't, modify the class type or the edgeList type
          static_assert(std::is_same<E, nodeIdType>::value, "types must match");

          //Parse the edgeList
          convertEdgeListforCCL(edgeList);

          //Re-distribute the tuples uniformly across the ranks
          mxx::distribute_inplace(tupleVector, comm);
        }

        /**
         * @brief   Compute the connected component labels
         * @note    Note that the communicator is freed after the computation
         */
        void compute()
        {
          //Size of vector should be >= 0
          assert(tupleVector.begin() != tupleVector.end());

          runConnectedComponentLabeling();
        }

        /**
         * @brief     count the components in the graph after ccl (useful for debugging/testing)
         * @note      should be called after computing connected components. 
         */
        std::size_t computeComponentCount()
        {
          std::size_t componentCount;

          //Vector should be sorted by Pc
          comm.with_subset(tupleVector.begin() !=  tupleVector.end() , [&](const mxx::comm& comm){

            if(!mxx::is_sorted(tupleVector.begin(), tupleVector.end(), conn::utils::TpleComp<cclTupleIds::Pc>(), comm))
              mxx::sort(tupleVector.begin(), tupleVector.end(), conn::utils::TpleComp<cclTupleIds::Pc>(), comm);

            //Count unique Pc values
            componentCount =  mxx::uniqueCount(tupleVector.begin(), tupleVector.end(),  conn::utils::TpleComp<cclTupleIds::Pc>(), comm);
          });

          componentCount = mxx::allreduce(componentCount, mxx::max<std::size_t>(), comm);

          return componentCount;
        }

      private:

        /**
         * @brief     Free the communicator
         * @note      Required to make sure that the communicator is freed before MPI_Finalize 
         */
        void free_comm()
        {
          comm.~comm();
        }

        /**
         * @brief     converts the edgelist to vector of tuples needed for ccl
         * @details   For the bucket in the edgeList ...<(u, v1), (u,v2)>...
         *            we append <(u, ~, u), (u, ~, v1), (u, ~, v2)> to our tupleVector
         *            We ignore the bucket splits across ranks here, because that shouldn't affect the correctness and complexity
         */
        template <typename edgeListPairsType>
        void convertEdgeListforCCL(edgeListPairsType &edgeList)
        {
          Timer timer(std::cerr, comm);

          //Sort the edgeList by src id of each edge
          mxx::sort(edgeList.begin(), edgeList.end(), conn::utils::TpleComp<edgeListTIds::src>(), comm); 

          //Reserve the approximate required space in our vector
          tupleVector.reserve(edgeList.size());

          for(auto it = edgeList.begin(); it != edgeList.end(); )
          {
            //Range of edges with same source vertex
            auto equalRange = conn::utils::findRange(it, edgeList.end(), *it, conn::utils::TpleComp<edgeListTIds::src>());

            //Range would include atleast 1 element
            assert(std::distance(equalRange.first, equalRange.second) > 0);

            //Insert other vertex members in this partition 
            for(auto it2 = equalRange.first; it2 != equalRange.second; it2++)
              tupleVector.emplace_back(std::get<edgeListTIds::src>(*it2), MAX_PID, std::get<edgeListTIds::dst>(*it2));;

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
          //variable to track convergence
          bool converged = false;

          //counting iterations
          int iterCount = 0;

          //range [tupleVector.begin() -- tupleVector.begin() + distance_begin_mid) marks the set of stable partitions in the vector
          //range [tupleVector.begin() + distance_begin_mid -- tupleVector.end()) denotes the active tuples 
          //Initially all the tuples are active, therefore we set distance_begin_mid to 0
          std::size_t distance_begin_mid = 0;

          while(!converged)
          {

            LOG_IF(comm.rank() == 0, INFO) << "Iteration #" << iterCount + 1;
            Timer timer(std::cerr, comm);

            //Temporary storage for extra tuples needed for doubling
            std::vector<T> parentRequestTupleVector;

            //Define the iterators over tupleVector
            auto begin = tupleVector.begin();
            auto mid = tupleVector.begin() + distance_begin_mid;
            auto end = tupleVector.end();

            //Log the min, mean and max count of active tuples across ranks
            //printWorkLoad(mid, end, comm);

            //Update Pn layer (Explore neighbors of a node and find potential partition candidates
            updatePn(mid, tupleVector.end());

            timer.end_section("\tPn update done");
            
            //Update the Pc layer, choose the best candidate
            converged = updatePc(mid, tupleVector.end(), parentRequestTupleVector);

            timer.end_section("\tPc update done");

            //Perform pointer doubling if enabled
            if(DOUBLING)
            {
              doPointerDoubling(distance_begin_mid, parentRequestTupleVector);

              //Due to insertion and deletion of elements, block decomposed property is lost during 
              //the pointer doubling, so redo it
              mxx::distribute_inplace(tupleVector, comm);

              timer.end_section("\tPointer doubling done");
            }



            //IMPORTANT : iterators over tupleVector are invalid and need to be redefined
            //Why? Because vector could undergo reallocation during pointer doubling
            
            begin = tupleVector.begin();
            mid = tupleVector.begin() + distance_begin_mid;
            end = tupleVector.end();

            //parition the dataset into stable and active paritions, if optimization is enabled
            if(!converged && (OPTIMIZATION == opt_level::stable_partition_removed || OPTIMIZATION == opt_level::loadbalanced))
            {
              //use std::partition to move stable tuples to the left
              mid = partitionStableTuples<cclTupleIds::Pn>(mid, end);

              timer.end_section("\tStable partitons placed aside");

              if(OPTIMIZATION == opt_level::loadbalanced)
              {
                mid = mxx::block_decompose_partitions(begin, mid, end, comm);
                //Re distributed the tuples to balance the load across the ranks
              
                timer.end_section("\tLoad balanced");
              }
            }
            distance_begin_mid = std::distance(begin, mid);

            iterCount ++;
          }

          LOG_IF(comm.rank() == 0, INFO) << "Algorithm took " << iterCount << " iterations";
        }

        /**
         * @brief             update the Pn layer by sorting the tuples using node ids
         * @param[in] begin   To iterate over the vector of tuples, marks the range of active tuples 
         * @param[in] end     end iterator
         */
        template <typename Iterator>
        void updatePn(Iterator begin, Iterator end)
        {
          comm.with_subset(begin != end, [&](const mxx::comm& com)
          {

              //Sort by nid,Pc
              mxx::sort(begin, end, conn::utils::TpleComp2Layers<cclTupleIds::nId, cclTupleIds::Pc>(), com); 

              //Resolve last and first bucket's boundary splits

              //First, find the element with max node id and min Pc locally 
              //Or in other words, get the min Pc of the last bucket
              auto minPcOfLastBucket = mxx::local_reduce(begin, end, conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>());

              //Second, do exscan, look for max nodeid and min Pc on previous ranks
              auto prevMinPc = mxx::exscan(minPcOfLastBucket, conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>(), com);  

              //We also need to know max Pc of the first bucket on the next rank (to check for stability)
              auto maxPcOfFirstBucket = mxx::local_reduce(begin, end, conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::less, std::greater>());

              //reverse exscan, look for min nodeid and max Pc on forward ranks
              auto nextMaxPc = mxx::exscan(maxPcOfFirstBucket,  conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::less, std::greater>(), com.reverse()); 

              //Now we can update the Pn layer of all the buckets locally
              for(auto it = begin; it !=  end;)
              {
                //Range of tuples with the same node id
                auto equalRange = conn::utils::findRange(it, end, *it, conn::utils::TpleComp<cclTupleIds::nId>());

                //Range would include atleast 1 element
                assert(std::distance(equalRange.first, equalRange.second) > 0);

                //Minimum Pc from local bucket
                auto thisBucketsMinPcLocal = mxx::local_reduce(equalRange.first, equalRange.second, conn::utils::TpleReduce<cclTupleIds::Pc>());

                //Maximum Pc from local bucket
                auto thisBucketsMaxPcLocal = mxx::local_reduce(equalRange.first, equalRange.second, conn::utils::TpleReduce<cclTupleIds::Pc, std::greater>());

                //For now, mark global minimum as local
                auto thisBucketsMaxPcGlobal = thisBucketsMaxPcLocal;
                auto thisBucketsMinPcGlobal = thisBucketsMinPcLocal;

                //Treat first, last buckets as special cases
                if(equalRange.first == begin)
                {
                  //Use value from previous rank
                  thisBucketsMinPcGlobal =  com.rank() == 0 ? thisBucketsMinPcLocal : conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>() (prevMinPc, thisBucketsMinPcLocal);
                }

                if(equalRange.second == end)
                {
                  //Use value from next rank
                  thisBucketsMaxPcGlobal = com.rank() == com.size() - 1 ? thisBucketsMaxPcLocal : conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::less, std::greater>() (nextMaxPc, thisBucketsMaxPcLocal);

                }

                auto maxPcValue = std::get<cclTupleIds::Pc>(thisBucketsMaxPcGlobal);
                auto minPcValue = std::min(std::get<cclTupleIds::Pc>(thisBucketsMinPcGlobal), std::get<cclTupleIds::nId>(*equalRange.first));

                //If min Pc < max Pc for this bucket, update Pn or else mark them as stable
                if(minPcValue < maxPcValue)
                  std::for_each(equalRange.first, equalRange.second, [&](T &e){
                      std::get<cclTupleIds::Pn>(e) = minPcValue;
                      });
                else
                  std::for_each(equalRange.first, equalRange.second, [&](T &e){
                      std::get<cclTupleIds::Pn>(e) = MAX_PID2;
                      });

                //Advance the loop pointer
                it = equalRange.second;
              }
          });
        }

        /**
         * @brief                             update the Pc layer by choosing min Pn
         * @param[in] begin                   To iterate over the vector of tuples, marks the range of active tuples 
         * @param[in] end                     end iterator
         * @param[in] partitionStableTuples   storate to keep 'parentRequest' tuples for doubling
         * @return                            bool value, true if the algorithm is converged
         */
        template <typename Iterator>
          bool updatePc(Iterator begin, Iterator end, std::vector<T>& parentRequestTupleVector)
          {
            //converged yet
            uint8_t converged = 1;    // 1 means true, we will update it below

            //Work only among ranks which have non-zero tuples left
            comm.with_subset(begin != end, [&](const mxx::comm& com)
            {
                //Sort by Pc, Pn
                mxx::sort(begin, end, conn::utils::TpleComp2Layers<cclTupleIds::Pc, cclTupleIds::Pn>(), com); 

                //Resolve last bucket's boundary split

                //First, find the element with max Pc and min Pn locally 
                //Or in other words, get the min Pn of the last bucket
                auto minPnOfLastBucket = mxx::local_reduce(begin, end, conn::utils::TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>());

                //Result of exscan, again look for max Pc and min Pn on previous ranks
                auto prevMinPn = mxx::exscan(minPnOfLastBucket, conn::utils::TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>(), com);  


                //Now we can update the Pc layer of all the buckets locally
                for(auto it = begin; it !=  end;)
                {
                  //Range of tuples with the same Pc
                  auto equalRange = conn::utils::findRange(it, end, *it, conn::utils::TpleComp<cclTupleIds::Pc>());

                  //Range would include atleast 1 element
                  assert(std::distance(equalRange.first, equalRange.second) > 0);

                  //Minimum Pn from local bucket
                  auto thisBucketsMinPnLocal = mxx::local_reduce(equalRange.first, equalRange.second, conn::utils::TpleReduce<cclTupleIds::Pn>());

                  //For now, mark global minimum as local
                  auto thisBucketsMinPnGlobal = thisBucketsMinPnLocal;

                  //Treat first, last buckets as special cases
                  if(equalRange.first == begin)
                  {
                    //Use value from previous rank
                    thisBucketsMinPnGlobal =  com.rank() == 0 ?  thisBucketsMinPnLocal : 
                      conn::utils::TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>() (prevMinPn, thisBucketsMinPnLocal);
                  }

                  //If min Pn < MAX_PID2 for this bucket, update the Pc to new value or else mark the partition as stable
                  if(std::get<cclTupleIds::Pn>(thisBucketsMinPnGlobal) < MAX_PID2) 
                  {

                    //Algorithm not converged yet because we found an active partition
                    converged = 0;

                    //Update Pc
                    std::for_each(equalRange.first, equalRange.second, [&](T &e){
                        std::get<cclTupleIds::Pc>(e) = std::get<cclTupleIds::Pn>(thisBucketsMinPnGlobal);
                        });

                    //Insert a 'parentRequest' tuple in the vector for doubling
                    if(DOUBLING)
                      parentRequestTupleVector.emplace_back(MAX_PID, MAX_PID, std::get<cclTupleIds::Pn>(thisBucketsMinPnGlobal));
                  }
                  else
                  {
                    //stable
                    std::for_each(equalRange.first, equalRange.second, [&](T &e){
                        std::get<cclTupleIds::Pn>(e) = MAX_PID;
                        });
                  }

                  //Advance the loop pointer
                  it = equalRange.second;
                }

            });

            //Know convergence of all the ranks
            uint8_t allConverged;
            mxx::allreduce(&converged, 1, &allConverged, mxx::min<uint8_t>(), comm);

            return (allConverged == 1  ? true : false);
          }

        /**
         * @brief                               Function that performs the pointer doubling
         *
         * @details                             
         *                                      'parentRequest' tuples serve the purpose of fetching parent of a partition
         *                                      In the beginning, they have the format <MAX_PID, MAX_PID, newPc>
         *                                      Goal is to 
         *                                        1.  Update the Pn layer of this tuple with the partition id of node newPc
         *                                            Since there could be multiple partition ids of a node 
         *                                            (because a node can belong to multiple partitions at an instant),
         *                                            we pick the minimum of them
         *                                            This would require cut-pasting all the tuples from parentRequestTupleVector to tupleVector
         *                                        2.  Flip the 'parentRequest' tuples and update the partition newPc to the value obtained above
         *                                        3.  Make sure to delete the 'parentRequest' tuples from tupleVector
         *
         * @param[in] beginOffset               tupleVector.begin() + beginOffset would denote the begin iterator for active tuples
         * @param[in] parentRequestTupleVector  All the 'parentRequest' tuples
         */
        void doPointerDoubling(std::size_t beginOffset, std::vector<T>& parentRequestTupleVector)
        {
          //Copy the tuples from parentRequestTupleVector to tupleVector 
          tupleVector.insert(tupleVector.end(), parentRequestTupleVector.begin(), parentRequestTupleVector.end());

          //Range of active tuples in tupleVector needs to be updated 
          auto begin = tupleVector.begin() + beginOffset;
          auto end = tupleVector.end();

          //Work among ranks with non-zero count of tuples
          comm.with_subset(begin != end, [&](const mxx::comm& com){

              //1. Repeat the procedure of updatePn, but just modify the 'parentRequest' tuples
              //   We can distinguish the 'parentRequest' tuples as they have Pc = MAX_PID

              //Same code as updatePn()
              mxx::sort(begin, end, conn::utils::TpleComp2Layers<cclTupleIds::nId, cclTupleIds::Pc>(), com); 
              auto minPcOfLastBucket = mxx::local_reduce(begin, end, conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>());
              auto prevMinPc = mxx::exscan(minPcOfLastBucket, conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>(), com);  
              for(auto it = begin; it !=  end;)
              {
                auto equalRange = conn::utils::findRange(it, end, *it, conn::utils::TpleComp<cclTupleIds::nId>());
                auto thisBucketsMinPcLocal = mxx::local_reduce(equalRange.first, equalRange.second, conn::utils::TpleReduce<cclTupleIds::Pc>());
                auto thisBucketsMinPcGlobal = thisBucketsMinPcLocal;
                if(equalRange.first == begin)
                {
                  thisBucketsMinPcGlobal =  com.rank() == 0 ? thisBucketsMinPcLocal : 
                  conn::utils::TpleReduce2Layers<cclTupleIds::nId, cclTupleIds::Pc, std::greater, std::less>() (prevMinPc, thisBucketsMinPcLocal);
                }

                std::for_each(equalRange.first, equalRange.second, [&](T &e)
                {
                  if(std::get<cclTupleIds::Pc>(e) == MAX_PID)
                  {
                    std::get<cclTupleIds::Pn>(e) = std::get<cclTupleIds::Pc>(thisBucketsMinPcGlobal);

                    //flip this 'parentRequest' tuple
                    std::get<cclTupleIds::Pc>(e) =  std::get<cclTupleIds::nId>(e);
                    std::get<cclTupleIds::nId>(e) =  MAX_NID;
                  }
                });

                it = equalRange.second;
              }

              //2. Now repeat the procedure of updatePc()
              mxx::sort(begin, end, conn::utils::TpleComp2Layers<cclTupleIds::Pc, cclTupleIds::Pn>(), com); 
              auto minPnOfLastBucket = mxx::local_reduce(begin, end, conn::utils::TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>());
              auto prevMinPn = mxx::exscan(minPnOfLastBucket, conn::utils::TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>(), com);  
              for(auto it = begin; it !=  end;)
              {
                auto equalRange = conn::utils::findRange(it, end, *it, conn::utils::TpleComp<cclTupleIds::Pc>());
                auto thisBucketsMinPnLocal = mxx::local_reduce(equalRange.first, equalRange.second, conn::utils::TpleReduce<cclTupleIds::Pn>());
                auto thisBucketsMinPnGlobal = thisBucketsMinPnLocal;
                if(equalRange.first == begin)
                {
                  thisBucketsMinPnGlobal =  com.rank() == 0 ?  thisBucketsMinPnLocal : 
                    conn::utils::TpleReduce2Layers<cclTupleIds::Pc, cclTupleIds::Pn, std::greater, std::less>() (prevMinPn, thisBucketsMinPnLocal);

                }

                //update the Pc for pointer jumping
                //Ignore the stable partitions
                if(std::get<cclTupleIds::Pn>(*equalRange.first) != MAX_PID)
                  std::for_each(equalRange.first, equalRange.second, [&](T &e){
                      std::get<cclTupleIds::Pc>(e) = std::get<cclTupleIds::Pn>(thisBucketsMinPnGlobal);
                      });

                it = equalRange.second;
              }
          });


          //3.  Now remove the 'parentRequest' tuples from tupleVector
          //    We can distinguish the 'parentRequest' tuples as they have nId = MAX_NID

          //operator should be '!=' because we want 'parentRequest' tuples to move towards right for deletion later
          auto mid = partitionStableTuples<cclTupleIds::nId, std::not_equal_to>(begin,end);

          //Erase the 'parentRequest' tuples
          tupleVector.erase(mid, end);
        }

        /*
         * @brief               partition the tuple array
         * @tparam[in]  layer   the layer of tuple used to partition them
         * @tparam[in]  op      binary operator that takes tuple_value and max
         *
         * @details             if op(value, max) is true, element is moved to left half                
         */
        template <uint8_t layer, template<typename> class op = std::equal_to, typename Iterator>
          Iterator partitionStableTuples(Iterator begin, Iterator end)
          {
            //type of the tuple_element at index layer
            using eleType = typename std::tuple_element<layer, T>::type;

            //max value
            auto max = std::numeric_limits<eleType>::max();

            //do the partition
            return std::partition(begin, end, [&](const T &e){
                return op<eleType>()(std::get<layer>(e), max);
                });
          }

        /**
         * @details   Helper function to print the load distribution i.e. the active tuples
         *            across ranks during the algorithm's exection. Prints the min, mean and 
         *            max count of the active tuples
         */
        template <typename Iterator, typename T = std::size_t>
          void printWorkLoad(Iterator begin, Iterator end, const mxx::comm &comm)
          {
            T localWorkLoad = std::distance(begin, end);

            T maxLoad  = mxx::reduce(localWorkLoad, 0, mxx::max<T>() , comm);
            T minLoad  = mxx::reduce(localWorkLoad, 0, mxx::min<T>() , comm);
            T meanLoad = mxx::reduce(localWorkLoad, 0, std::plus<T>(), comm)/ comm.size();

            auto sep = ",";
            LOG_IF(comm.rank() == 0, INFO) << "Load distribution of active tuples min-mean-max : " << minLoad << sep << meanLoad << sep << maxLoad;
          }

        /**
         * @brief     Print verbose log of tuple counts on all the ranks (both active and inactive)
         * @note      Use only while debugging
         */
        template <typename Iterator, typename T = std::size_t>
          void printVerboseTupleCounts(Iterator begin, Iterator mid, Iterator end)
          {
            T inactiveTupleCount = std::distance(begin, mid);
            T activeTupleCount = std::distance(mid, end);

            std::pair<T,T> tupleCounts = std::make_pair(inactiveTupleCount, activeTupleCount);

            //Gather to rank 0
            auto gatherValues = mxx::gather(tupleCounts, 0, comm);

            //Print the pairs 
            if(!comm.rank()) std::cerr << gatherValues << std::endl; 
          }

    };

  }
}

#endif
