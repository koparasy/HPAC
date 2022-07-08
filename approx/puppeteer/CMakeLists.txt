cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(APPROX_BENCH LANGUAGES C CXX)

option(WITH_VISUALIZE "Plot raw data to the respectibe paper figures" OFF)
option(WITH_SRC_CODE "Create build directories of puppetteer" ON)

IF( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE Release FORCE )
ENDIF()

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ standard to conform to")
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(CMAKE_UTIL_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LIBAPPROX_LIB_FILE "ApproxUtil")

if (WITH_SRC_CODE)
  add_subdirectory(libraries)
  add_subdirectory(benchmarks)
endif()

if (WITH_VISUALIZE)
  add_subdirectory(data)
endif()

