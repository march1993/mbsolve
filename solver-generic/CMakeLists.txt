cmake_minimum_required(VERSION 3.7)
project(solver-generic)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/
  ${solver_generic_BINARY_DIR}/exports/)

add_library(solver-generic SHARED src/SolverGeneric.cpp)

include(GenerateExportHeader)
generate_export_header(solver-generic
  BASE_NAME solver-generic
  EXPORT_MACRO_NAME SOLVER_GENERIC_EXPORT
  EXPORT_FILE_NAME exports/solver-generic_EXPORTS.hpp
  STATIC_DEFINE SOLVER_GENERIC_BUILT_AS_STATIC)

target_include_directories(solver-generic PUBLIC include/
    ${solver_generic_BINARY_DIR}/exports/)
target_link_libraries(solver-generic mbsolve-lib)
