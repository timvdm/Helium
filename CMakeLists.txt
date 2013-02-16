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
find_package(Boost 1.48 COMPONENTS timer system REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

# Find and setup Eigen3 (optional)
find_package(Eigen3)
if (EIGEN3_FOUND)
  message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Find OpenCL (optional)
option(ENABLE_OPENCL "Enable OpenCL support" OFF)
if (ENABLE_OPENCL)
  find_package(OpenCL)
  if (OPENCL_FOUND)
    message(STATUS "OpenCL found: enabling OpenCL features...")
    add_definitions(-DHAVE_OPENCL)
    include_directories(${OPENCL_INCLUDE_DIR})
  endif()
endif()

############################################################
#
# Compiler configuration
#
############################################################

# compiler options
set(CMAKE_CXX_FLAGS "-g -pedantic -Wall -O2 -Wno-sign-compare ${CMAKE_CXX_FLAGS}")

# C++11 (optional) [only seems to work with clang on mac for now]
option(ENABLE_CPP11 "Enable C++11 features" OFF)
if (ENABLE_CPP11)
  if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}")
  endif()
  add_definitions(-DHAVE_CPP11)
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

set(Helium_HDRS
  src/bitvec.h
  src/canonical.h
  src/components.h
  src/enumeratepaths.h
  src/enumeratesubgraphs.h
  src/extendedconnectivities.h
  src/fileio.h
  src/fingerprints.h
  src/invariants.h
  src/molecule.h
  src/substructure.h
  src/tie.h
  src/timeout.h
  src/util.h
  src/hemol.h
)

set(Helium_util_HDRS
  src/util/fileio.h
  src/util/functor.h
  src/util/string.h
  src/util/typetraits.h
  src/util/vector.h
)

set(Helium_fileio_HDRS
  src/fileio/fingerprints.h
)

set(Helium_SRCS
  src/hemol.cpp
  src/fileio/file.cpp
  src/thirdparty/jsoncpp/jsoncpp.cpp
)

# Helium library
add_library(helium SHARED ${Helium_SRCS})
target_link_libraries(helium ${Boost_LIBRARIES} ${OPENCL_LIBRARIES})

set(Helium_LIBRARIES
  helium
  pthread
)

# install target for header files
set(INCLUDE_INSTALL_DIR  "include" CACHE STRING "Install dir for include files")
install(FILES ${Helium_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium")
install(FILES ${Helium_util_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/util")
install(FILES ${Helium_fileio_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/fileio")

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
