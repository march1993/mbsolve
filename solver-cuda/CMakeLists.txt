cmake_minimum_required(VERSION 3.7)
project(solver-cuda)

set(CUDA_SEPARABLE_COMPILATION ON)

cuda_add_library(solver-cuda SHARED src/CUDADensityMatrix.cu
  src/SolverCUDA2lvl.cu include/SolverCUDA2lvl.hpp
  include/CUDADensityMatrix.hpp include/CUDACommon.hpp)

target_include_directories(solver-cuda PUBLIC include/)
target_link_libraries(solver-cuda mbsolve-lib)
