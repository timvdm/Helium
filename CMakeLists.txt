cmake_minimum_required(VERSION 2.6)
project(Helium)

# set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

############################################################
#
# Find dependencies
#
############################################################

# Setup 3th party packages (included with source)
include_directories(${CMAKE_SOURCE_DIR}/src/thirdparty)
include_directories(${CMAKE_SOURCE_DIR}/src/thirdparty/jsoncpp)

# Find and setup boost (required)
find_package(Boost 1.48 COMPONENTS timer system iostreams REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

# Find and setup Eigen3 (optional)
find_package(Eigen3)
if (EIGEN3_FOUND)
  message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Find OpenCL (optional)
option(ENABLE_OPENCL "Enable OpenCL support" OFF)
set(HAVE_OPENCL 0) # used by configure_file()
if (ENABLE_OPENCL)
  find_package(OpenCL)
  if (OPENCL_FOUND)
    message(STATUS "OpenCL found: enabling OpenCL features...")
    include_directories(${OPENCL_INCLUDE_DIR})
    add_definitions(-DHAVE_OPENCL)
    set(HAVE_OPENCL 1) # used by configure_file()
  endif()
endif()

############################################################
#
# Compiler configuration
#
############################################################

# compiler options
set(CMAKE_CXX_FLAGS "-g -pedantic -Wall -Wno-long-long -Wno-sign-compare ${CMAKE_CXX_FLAGS}")

# C++11 (optional) [only seems to work with clang on mac for now]
option(ENABLE_CPP11 "Enable C++11 features" OFF)
if (ENABLE_CPP11)
  if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}")
  endif()
  add_definitions(-DHAVE_CPP11)
  set(HAVE_CPP11 1) # used by configure_file()
else()
  set(HAVE_CPP11 0) # used by configure_file()
endif()

if (ENABLE_COVERAGE)
  set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage ${CMAKE_CXX_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS "-O2 ${CMAKE_CXX_FLAGS}")
endif()

############################################################
#
# Platform checks
#
############################################################

# platform check for POPCNT instruction
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("int main() { __builtin_popcountll(0); }" HAVE_POPCNT)
if (HAVE_POPCNT)
  add_definitions(-DHAVE_POPCNT)
endif()

############################################################
#
# The Helium library
#
############################################################

configure_file(${CMAKE_SOURCE_DIR}/src/config.h.cmake
  ${CMAKE_BINARY_DIR}/include/Helium/config.h)

include_directories(${CMAKE_BINARY_DIR}/include)
add_subdirectory(src)

set(Helium_LIBRARIES
  helium
  pthread
)

############################################################
#
# Unit Testing & Tools
#
############################################################


# unit tests
enable_testing()
include(CTest)
set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
add_subdirectory(test)

add_subdirectory(tools)

############################################################
#
# Doxygen Documentation
#
############################################################

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif()

############################################################
#
# Scripting Bindings
#
############################################################

# Python
option(ENABLE_PYTHON "Enable python bindings" OFF)
if (ENABLE_PYTHON)
  add_subdirectory(bindings/python)
endif()
