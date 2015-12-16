/**
 * @file    labelProp_utils.hpp
 * @ingroup coloring
 * @author  Chirag Jain <cjain7@gatech.edu>
 * @brief   Utils for labelProp.hpp
 *
 * Copyright (c) 2015 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef LABEL_PROPAGATION_UTILS_HPP 
#define LABEL_PROPAGATION_UTILS_HPP



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
      nId         //Node id , remains unchanged
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
     * @brief     optimization level for the coloring algorithm 
     */
    enum opt_level
    {
      naive,          
      stable_partition_removed,  //removes stable partitions from the working set
      loadbalanced    //enables load balance, recommended setting, used by default
    };

    /**
     * @brief   On and off switch
     */
    enum lever
    {
      OFF,
      ON
    };
  }
}

#endif

