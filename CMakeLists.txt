cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(cluster-mc CXX)
include(cmake/postfix.cmake)

include_directories(${PROJECT_SOURCE_DIR})

include(eigen3)
include(standards)
include(lattice)
FetchContent_MakeAvailable(${FetchContents})
include_directories(${FetchContent_includes})

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

# IO Test
include(add_iotest)
enable_testing()

option(BUILD_STANDALONE_APPLICATIONS "Build standalone applications" ON)
option(BUILD_ALPS_APPLICATIONS "Build ALPS applications" ON)

add_subdirectory(test)

if(BUILD_STANDALONE_APPLICATIONS)
  add_subdirectory(standalone)
endif(BUILD_STANDALONE_APPLICATIONS)

if(ALPS_VERSION AND BUILD_ALPS_APPLICATIONS)
  add_subdirectory(alps)
endif(ALPS_VERSION AND BUILD_ALPS_APPLICATIONS)
