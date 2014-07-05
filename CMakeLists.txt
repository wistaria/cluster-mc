#
# Cluster-MC: Cluster Algorithm Monte Carlo Methods
#
# Copyright (C) 1997-2014 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(inifinite_range NONE)

option(BUILD_STANDALONE_APPLICATIONS "Build standalone applications" ON)
option(BUILD_ALPS_APPLICATIONS "Build ALPS applications" ON)

set(ALPS_HOME_ENV $ENV{ALPS_HOME})
if((ALPS_ROOT_DIR OR ALPS_HOME_ENV) AND BUILD_ALPS_APPLICATIONS)
  find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME}
    NO_SYSTEM_ENVIRONMENT_PATH)
  message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
  include(${ALPS_USE_FILE})
else((ALPS_ROOT_DIR OR ALPS_HOME_ENV) AND BUILD_ALPS_APPLICATIONS)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
  find_package(BoostSrc)
  include_directories(${Boost_ROOT_DIR})
  include(add_alps_test)
endif((ALPS_ROOT_DIR OR ALPS_HOME_ENV) AND BUILD_ALPS_APPLICATIONS)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# enable C and C++ compilers
enable_language(C CXX)

enable_testing()
include(CTest)

include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

if(BUILD_STANDALONE_APPLICATIONS)
  add_subdirectory(standalone)
endif(BUILD_STANDALONE_APPLICATIONS)

if(BUILD_ALPS_APPLICATIONS AND ALPS_LIBRARIES)
  add_subdirectory(alps)
endif(BUILD_ALPS_APPLICATIONS AND ALPS_LIBRARIES)
