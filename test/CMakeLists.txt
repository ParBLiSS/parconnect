cmake_minimum_required(VERSION 2.6)

# project settings
project(parConnect-test)

#################
## PERFORMANCE ##
#################

add_executable(time_graph500_gen time_graph500_gen.cpp)
target_link_libraries(time_graph500_gen ${EXTRA_LIBS})

add_executable(time_ccl_coloring_graph500  time_ccl_coloring_graph500.cpp)
target_link_libraries(time_ccl_coloring_graph500 ${EXTRA_LIBS})

add_executable(time_ccl_coloring_chains time_ccl_coloring_chains.cpp)
target_link_libraries(time_ccl_coloring_chains ${EXTRA_LIBS})

add_executable(time_bfs_graph500 time_bfs_graph500.cpp)
target_link_libraries(time_bfs_graph500 ${EXTRA_LIBS} MPITypelib CommGridlib)

add_executable(time_bfs_coloring_graph500 time_bfs_coloring_graph500.cpp)
target_link_libraries(time_bfs_coloring_graph500  ${EXTRA_LIBS} MPITypelib CommGridlib)

#################
## CORRECTNESS ##
#################

add_executable(test-coloring test_ccl_coloring.cpp)
target_link_libraries(test-coloring mxx-gtest-main)

add_executable(test-graphgen test_graphgen.cpp)
target_link_libraries(test-graphgen mxx-gtest-main)

add_executable(test-bfsRunner test_bfsRunner.cpp)
target_link_libraries(test-bfsRunner mxx-gtest-main MPITypelib CommGridlib)
