cmake_minimum_required(VERSION 3.7)
project(solver-generic)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

add_library(solver-generic SHARED src/SolverGeneric.cpp)
target_include_directories(solver-generic PUBLIC include/)
target_link_libraries(solver-generic mbsolve-lib)
