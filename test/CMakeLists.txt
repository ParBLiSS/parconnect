cmake_minimum_required(VERSION 2.6)

# project settings
project(parConnect-test)

##### add definitions to find path to data files
add_definitions( -DPROJECT_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#################
## PERFORMANCE ##
#################

#add_executable(benchmark_parconnect benchmark_parconnect.cpp)
#target_link_libraries(benchmark_parconnect ${EXTRA_LIBS} MPITypelib CommGridlib)

add_executable(parconnect benchmark_parconnect_auto.cpp)
target_link_libraries(parconnect ${EXTRA_LIBS} MPITypelib CommGridlib plfit0)

#add_executable(benchmark_sequential benchmark_sequential.cpp)
#target_link_libraries(benchmark_sequential ${EXTRA_LIBS})

#add_executable(time_all2all time_all2all.cpp)
#target_link_libraries(time_all2all ${EXTRA_LIBS})

#add_executable(time_sort time_sort.cpp)
#target_link_libraries(time_sort ${EXTRA_LIBS})

#################
## CORRECTNESS ##
#################

OPTION(BUILD_CONN_UTILS "Build the utilities utils_*.cpp" OFF)

if(BUILD_CONN_UTILS)
  add_executable(exportGraphDotFormat utils_exportGraphDotFormat.cpp)
  target_link_libraries(exportGraphDotFormat ${EXTRA_LIBS})

  add_executable(exportBinaryFormat utils_exportBinaryFormat.cpp)
  target_link_libraries(exportBinaryFormat ${EXTRA_LIBS})
endif(BUILD_CONN_UTILS)

OPTION(BUILD_CONN_TESTS "Build the test code for checking parallel connectivity" OFF)

if(BUILD_CONN_TESTS)
  add_executable(test-coloring test_ccl_coloring.cpp)
  target_link_libraries(test-coloring mxx-gtest-main)

  add_executable(test-graphgen test_graphgen.cpp)
  target_link_libraries(test-graphgen mxx-gtest-main)

  add_executable(test-bfsRunner test_bfsRunner.cpp)
  target_link_libraries(test-bfsRunner mxx-gtest-main MPITypelib CommGridlib)
endif(BUILD_CONN_TESTS)
