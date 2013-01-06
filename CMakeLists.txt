cmake_minimum_required(VERSION 2.6)
project(Helium)

# set the module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# find and setup OpenBabel
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

# compiler options
set(CMAKE_CXX_FLAGS "-g -pedantic -Wall -O2 -Wno-sign-compare ${CMAKE_CXX_FLAGS}")


set(Helium_HDRS
  src/molecule.h
  )

# install target for smiley header file
install(FILES Helium_HDRS DESTINATION include/Helium)


# unit test
enable_testing()
set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
add_subdirectory(test)

add_subdirectory(tools)
