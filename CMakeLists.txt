#
# Cluster-MC: Cluster Algorithm Monte Carlo Methods
#
# Copyright (C) 1997-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(cluster_mc NONE)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
include(core)

# Eigen3
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

# IO Test
include(add_iotest)
enable_testing()

option(BUILD_STANDALONE_APPLICATIONS "Build standalone applications" ON)
option(BUILD_ALPS_APPLICATIONS "Build ALPS applications" ON)

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/tool/standards ${PROJECT_SOURCE_DIR}/tool/lattice)

if(BUILD_STANDALONE_APPLICATIONS)
  add_subdirectory(standalone)
endif(BUILD_STANDALONE_APPLICATIONS)

if(ALPS_VERSION AND BUILD_ALPS_APPLICATIONS)
  add_subdirectory(alps)
endif(ALPS_VERSION AND BUILD_ALPS_APPLICATIONS)
