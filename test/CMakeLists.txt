set(tests
  molecule
  fileio
  enumeratepaths
  enumeratesubgraphs
  canonical
  util
  components
  fingerprints
  #maximalmatching
  isomorphism
  cycles
  smiles
  bitvec
  distancematrix
  dijkstra
  floydwarshall
  razinger
  dfs
  bfs
  smarts
  smirks
  smartcomponents
  )

foreach(test ${tests})
  add_executable(test_${test} ${test}.cpp)
  target_link_libraries(test_${test} ${Helium_LIBRARIES})
  add_test(${test}_Test ${TEST_PATH}/test_${test})
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "FAIL")
endforeach(test ${tests})

find_package(OpenBabel2)
if (OPENBABEL2_FOUND)
  message(STATUS "OpenBabel found: building openbabel example")
  include_directories(${OPENBABEL2_INCLUDE_DIR})
  add_executable(test_openbabel openbabel.cpp)
  target_link_libraries(test_openbabel ${Helium_LIBRARIES} ${OPENBABEL2_LIBRARIES})
endif()

if (ENABLE_PYTHON)
  add_subdirectory(python)
endif()
