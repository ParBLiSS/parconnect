cmake_minimum_required(VERSION 2.6)

# project settings
project(parConnect-test)

#################
## PERFORMANCE ##
#################

add_executable(time_graph500_gen time_graph500_gen.cpp)
target_link_libraries(time_graph500_gen ${EXTRA_LIBS})

add_executable(time_ccl_coloring_undirected_graph500  time_ccl_coloring_undirected_graph500.cpp)
target_link_libraries(time_ccl_coloring_undirected_graph500 ${EXTRA_LIBS})

#################
## CORRECTNESS ##
#################

add_executable(test-coloring test_ccl_coloring.cpp)
target_link_libraries(test-coloring mxx-gtest-main)
