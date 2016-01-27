cmake_minimum_required(VERSION 2.6)

# project settings
project(parConnect-test)

##### add definitions to find path to data files
add_definitions( -DPROJECT_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

#################
## PERFORMANCE ##
#################

#add_executable(time_bfs_coloring_graph500 time_bfs_coloring_graph500.cpp)
#target_link_libraries(time_bfs_coloring_graph500  ${EXTRA_LIBS} MPITypelib CommGridlib)

#add_executable(time_bfs_coloring_debruijn_graph time_bfs_coloring_debruijn_graph.cpp)
#target_link_libraries(time_bfs_coloring_debruijn_graph ${EXTRA_LIBS} MPITypelib CommGridlib)

#add_executable(time_bfs_coloring_fileInput time_bfs_coloring_fileInput.cpp)
#target_link_libraries(time_bfs_coloring_fileInput ${EXTRA_LIBS} MPITypelib CommGridlib)

add_executable(benchmark_parconnect benchmark_parconnect.cpp)
target_link_libraries(benchmark_parconnect ${EXTRA_LIBS} MPITypelib CommGridlib)

#################
## CORRECTNESS ##
#################

add_executable(test-coloring test_ccl_coloring.cpp)
target_link_libraries(test-coloring mxx-gtest-main)

add_executable(test-graphgen test_graphgen.cpp)
target_link_libraries(test-graphgen mxx-gtest-main)

add_executable(test-bfsRunner test_bfsRunner.cpp)
target_link_libraries(test-bfsRunner mxx-gtest-main MPITypelib CommGridlib)
