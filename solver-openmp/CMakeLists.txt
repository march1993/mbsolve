cmake_minimum_required(VERSION 3.7)
project(solver-openmp)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/
  ${solver-openmp_BINARY_DIR}/exports/)

add_library(solver-openmp SHARED
  src/solver_2lvl_pc.cpp
  src/solver_3lvl_pc.cpp
  src/DensityMatrixData.cpp)

include(GenerateExportHeader)
generate_export_header(solver-openmp
  BASE_NAME solver-openmp
  EXPORT_MACRO_NAME SOLVER_OPENMP_EXPORT
  EXPORT_FILE_NAME exports/solver-openmp_EXPORTS.hpp
  STATIC_DEFINE SOLVER_OPENMP_BUILT_AS_STATIC)

target_compile_options(solver-openmp PRIVATE ${OpenMP_CXX_FLAGS})

target_include_directories(solver-openmp PUBLIC include/
  ${solver-openmp_BINARY_DIR}/exports/)
target_link_libraries(solver-openmp ${OpenMP_CXX_FLAGS})
target_link_libraries(solver-openmp mbsolve-lib)
