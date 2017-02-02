cmake_minimum_required(VERSION 3.7)
project(writer-matlab)

include_directories(${mbsolve-lib_INCLUDE_DIR} ${Matlab_INCLUDE_DIRS} include/
  ${writer-matlab_BINARY_DIR}/exports/)

add_library(writer-matlab SHARED src/WriterMATLAB.cpp)

include(GenerateExportHeader)
generate_export_header(writer-matlab
  BASE_NAME writer-matlab
  EXPORT_MACRO_NAME WRITER_MATLAB_EXPORT
  EXPORT_FILE_NAME exports/writer-matlab_EXPORTS.hpp
  STATIC_DEFINE WRITER_MATLAB_BUILT_AS_STATIC)

target_include_directories(writer-matlab PUBLIC include/
  ${writer-matlab_BINARY_DIR}/exports/)
target_link_libraries(writer-matlab mbsolve-lib)
target_link_libraries(writer-matlab ${Matlab_LIBRARIES})
