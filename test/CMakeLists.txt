cmake_minimum_required(VERSION 2.6)

# project settings
project(parConnect-test)

##### add definitions to find path to data files
add_definitions( -DPROJECT_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#################
## PERFORMANCE ##
#################

add_executable(benchmark_parconnect benchmark_parconnect.cpp)
target_link_libraries(benchmark_parconnect ${EXTRA_LIBS} MPITypelib CommGridlib)

add_executable(time_all2all time_all2all.cpp)
target_link_libraries(time_all2all ${EXTRA_LIBS})

#################
## CORRECTNESS ##
#################

OPTION(BUILD_CONN_TESTS "Build the test code for checking parallel connectivity" OFF)

if(BUILD_CONN_TESTS)
  add_executable(test-coloring test_ccl_coloring.cpp)
  target_link_libraries(test-coloring mxx-gtest-main)

  add_executable(test-graphgen test_graphgen.cpp)
  target_link_libraries(test-graphgen mxx-gtest-main)

  add_executable(test-bfsRunner test_bfsRunner.cpp)
  target_link_libraries(test-bfsRunner mxx-gtest-main MPITypelib CommGridlib)
endif(BUILD_CONN_TESTS)
