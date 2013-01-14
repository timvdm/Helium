cmake_minimum_required(VERSION 2.6)
project(Helium)

# set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# find and setup OpenBabel
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

#find and setup JsonCpp
set(JsonCpp_DIR "/usr/local/lib")
find_package(JsonCpp)
include_directories(${JsonCpp_INCLUDE_DIRS})

# FIXME
include_directories("/usr/local/include/eigen3")

set(Helium_LIBRARIES
  ${OPENBABEL2_LIBRARIES}
  ${JsonCpp_LIBRARIES}
  boost_timer
  boost_system
  pthread
)


# compiler options
set(CMAKE_CXX_FLAGS "-g -pedantic -Wall -O2 -Wno-sign-compare ${CMAKE_CXX_FLAGS}")


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

# install target for smiley header file
set(INCLUDE_INSTALL_DIR  "include" CACHE STRING "Install dir for include files")
install(FILES ${Helium_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium")
install(FILES ${Helium_util_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/util")

# unit test
enable_testing()
set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
add_subdirectory(test)

add_subdirectory(tools)


