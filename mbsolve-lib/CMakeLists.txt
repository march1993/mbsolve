cmake_minimum_required(VERSION 3.7)
project(mbsolve-lib)

include_directories(include/ ${mbsolve-lib_BINARY_DIR}/exports/)

add_library(mbsolve-lib SHARED
  src/Quantity.cpp src/Solver.cpp src/Writer.cpp
  include/Device.hpp include/Quantity.hpp
  include/Scenario.hpp include/Solver.hpp
  include/Types.hpp include/Writer.hpp)
target_include_directories(mbsolve-lib PUBLIC include/
  ${mbsolve-lib_BINARY_DIR}/exports/)
include(GenerateExportHeader)
generate_export_header(mbsolve-lib
  BASE_NAME mbsolve-lib
  EXPORT_MACRO_NAME MBSOLVE_LIB_EXPORT
  EXPORT_FILE_NAME exports/mbsolve-lib_EXPORTS.hpp
  STATIC_DEFINE MBSOLVE_LIB_BUILT_AS_STATIC)
