set(tests
  molecule
  fileio
  enumeratepaths
  enumeratesubgraphs
  canonical
  substructure
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
  )

foreach(test ${tests})
  add_executable(test_${test} ${test}.cpp)
  target_link_libraries(test_${test} ${Helium_LIBRARIES})
  add_test(${test}_Test ${TEST_PATH}/test_${test})
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "FAIL")
endforeach(test ${tests})


