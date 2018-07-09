cmake_minimum_required(VERSION 3.6)
project(solver-mpich)

# specify matrix exponential method
set(MAT_EXP_METHOD "DIAG" CACHE STRING "specify matrix exponential method")
set_property(CACHE MAT_EXP_METHOD PROPERTY STRINGS "PADE" "DIAG" "RODR")

string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOLVER_MPICH_SRCS
  src/solver_mpich_clvl_os_red.cpp
)

add_library(solver-mpich SHARED ${SOLVER_MPICH_SRCS})

set_property(TARGET solver-mpich PROPERTY CXX_STANDARD 11)

if(ARCH STREQUAL "AVX2")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-mpich PRIVATE "-mavx2")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-mpich PRIVATE "-xAVX2")
  else()
  endif()
elseif(ARCH STREQUAL "SSE3")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-mpich PRIVATE "-msse3")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-mpich PRIVATE "-xSSE3")
  else()
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-mpich PRIVATE "-march=native")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-mpich PRIVATE "-xHost")
  else()
  endif()
endif()

if(MAT_EXP_METHOD STREQUAL "DIAG")
  target_compile_definitions(solver-mpich PRIVATE EXP_METHOD=1)
elseif(MAT_EXP_METHOD STREQUAL "RODR")
  target_compile_definitions(solver-mpich PRIVATE EXP_METHOD=2)
else()
endif()

target_compile_options(solver-mpich PRIVATE ${MPICH})

target_include_directories(solver-mpich PUBLIC include/)
target_link_libraries(solver-mpich ${MPICH_CXX_FLAGS})
target_link_libraries(solver-mpich mbsolve-lib)

target_include_directories(solver-mpich PUBLIC ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(solver-mpich ${MPI_CXX_LIBRARIES})
