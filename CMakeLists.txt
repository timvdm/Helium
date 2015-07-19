cmake_minimum_required(VERSION 2.6)
project(Helium)

set(HELIUM_MAJOR 0)
set(HELIUM_MINOR 2)
set(HELIUM_PATCH 0)

# set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# define HEDATADIR for tests that need input files
set(HEDATADIR "${CMAKE_SOURCE_DIR}/data/")

############################################################
#
# Find dependencies
#
############################################################

include_directories(${CMAKE_BINARY_DIR}/include)

# Setup 3th party packages (included with source)
include_directories(${CMAKE_SOURCE_DIR}/src/thirdparty)

# Find and setup boost (required)
if (MSVC)
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-DBOOST_PYTHON_STATIC_LIB)
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_MULTITHREADED ON)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()
find_package(Boost 1.48 COMPONENTS chrono timer system iostreams filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Find and setup Eigen3 (required)
if (MSVC)
  add_definitions(-DEIGEN_DONT_ALIGN_STATICALLY=1)
endif()
find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

option(ENABLE_OPENBABEL "Enable OpenBabel support" OFF)
if (ENABLE_OPENBABEL)
  find_package(OpenBabel2)
endif()

# Find OpenCL (optional)
option(ENABLE_OPENCL "Enable OpenCL support" OFF)
set(HAVE_OPENCL 0) # used by configure_file()
if (ENABLE_OPENCL)
  find_package(OpenCL)
  if (OPENCL_FOUND)
    message(STATUS "OpenCL found: enabling OpenCL features...")
    include_directories(${OPENCL_INCLUDE_DIR})
    set(HAVE_OPENCL 1) # used by configure_file()
  endif()
endif()

############################################################
#
# Compiler configuration
#
############################################################

# compiler options
if (MSVC)
  set(CompilerFlags
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE
  )
  foreach(CompilerFlag ${CompilerFlags})
    string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
  endforeach()
else()
  set(CMAKE_CXX_FLAGS "-g -pedantic -Wall -Wno-long-long -Wno-sign-compare -Wno-deprecated-register  ${CMAKE_CXX_FLAGS}")
endif()

# C++11
if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
elseif(UNIX)
  # asume clang++
  set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}")
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

add_subdirectory(src)

set(Helium_LIBRARIES
  helium
  ${Boost_LIBRARIES}
)

if (UNIX)
  set(Helium_LIBRARIES
    ${Helium_LIBRARIES}
    pthread
  )
endif()

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

# Build service
option(ENABLE_SERVICE "Enable the RESTful service" OFF)
if (ENABLE_SERVICE)
  add_subdirectory(service)
endif()

############################################################
#
# CPack
#
############################################################

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Helium")
set(CPACK_PACKAGE_VENDOR "Tim Vandermeersch")
set(CPACK_PACKAGE_VERSION_MAJOR ${HELIUM_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${HELIUM_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${HELIUM_PATCH})
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "helium-${HELIUM_MAJOR}.${HELIUM_MINOR}.${HELIUM_PATCH}")

set(ignore_files
  .git
  .gitignore
  .DS_Store
  build
  build-*
  devel
  docs
  examples
  data/100K.smi
  data/100K.hel
  data/1M.smi
  data/1M.hel
)

set(CPACK_SOURCE_IGNORE_FILES "${ignore_files}")

include(CPack)
