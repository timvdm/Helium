# Find and setup python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Find and setup boost (required)
find_package(Boost 1.48 COMPONENTS python iostreams REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

message("python executable: ${PYTHON_EXECUTABLE}")
message("python library: ${PYTHON_LIBRARY}")

set(python_SRCS
  main.cpp
  bitvec.cpp
  common.cpp
  diagram.cpp
  distancematrix.cpp
  element.cpp
  error.cpp
  molecule.cpp
  rings.cpp
  smarts.cpp
  smiles.cpp
  smirks.cpp
  algorithms/bfs.cpp
  algorithms/canonical.cpp
  algorithms/dfs.cpp
  algorithms/dijkstra.cpp
  algorithms/enumeratepaths.cpp
  algorithms/enumeratesubgraphs.cpp
  algorithms/extendedconnectivities.cpp
  algorithms/floydwarshall.cpp
  algorithms/components.cpp
  algorithms/invariants.cpp
  fileio/file.cpp
  fileio/fingerprints.cpp
  fileio/fps.cpp
  fileio/moleculefile.cpp
  fingerprints/fingerprints.cpp
  fingerprints/similarity.cpp
  depict/depict.cpp
  ../../src/smarts.cpp
  ../../src/smirks.cpp
  ../../src/fileio/file.cpp
  ../../src/jsoncpp.cpp
  ../../src/diagram.cpp
  ../../src/depict/svgpainter.cpp
  ../../src/depict/depict.cpp
)

set(CMAKE_CXX_FLAGS "-O2")

add_library(PyHelium SHARED ${python_SRCS})
set_target_properties(PyHelium PROPERTIES
  PREFIX ""
  SUFFIX ".so"
  OUTPUT_NAME helium)
target_link_libraries(PyHelium ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
if (NOT WIN32)
  execute_process(COMMAND
    ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE PYTHON_INSTDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  install(TARGETS PyHelium
    LIBRARY DESTINATION ${PYTHON_INSTDIR}
    COMPONENT PyHelium)
endif()

