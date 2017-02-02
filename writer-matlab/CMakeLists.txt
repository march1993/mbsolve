cmake_minimum_required(VERSION 3.7)
project(writer-matlab)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${Matlab_INCLUDE_DIRS} include/)

add_library(writer-matlab SHARED src/WriterMATLAB.cpp)
target_include_directories(writer-matlab PUBLIC include/)
target_link_libraries(writer-matlab mbsolve-lib)
target_link_libraries(writer-matlab ${Matlab_LIBRARIES})
