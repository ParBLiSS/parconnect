cmake_minimum_required(VERSION 2.6)

# project settings
project(parConnect-test)

add_executable(check_graph500_gen check_graph500_gen.cpp)
target_link_libraries(check_graph500_gen ${EXTRA_LIBS})
