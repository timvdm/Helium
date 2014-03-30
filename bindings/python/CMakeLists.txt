# Find and setup python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Find and setup boost (required)
find_package(Boost 1.48 COMPONENTS python REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

message("python executable: ${PYTHON_EXECUTABLE}")
message("python library: ${PYTHON_LIBRARY}")

set(python_SRCS
  main.cpp
  common.cpp
  dfs.cpp
  element.cpp
  error.cpp
  molecule.cpp
  rings.cpp
  smarts.cpp
  smiles.cpp
  smirks.cpp
  substructure.cpp
  ../../src/smarts.cpp
  ../../src/smirks.cpp
)

set(CMAKE_CXX_FLAGS "-O2")

add_library(PyHelium SHARED ${python_SRCS})
set_target_properties(PyHelium PROPERTIES OUTPUT_NAME helium)
target_link_libraries(PyHelium ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
