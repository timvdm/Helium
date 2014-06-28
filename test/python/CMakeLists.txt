include(UsePythonTest)

set(python_tests
  bfs
  bitvec
  canonical
  components
  depict
  dfs
  dijkstra
  element
  enumeratepaths
  enumeratesubgraphs
  extendedconnectivities
  file
  fingerprintfiles
  fingerprints
  floydwarshall
  fps
  invariants
  molecule
  moleculefile
  rings
  similarity
  smarts
  smiles
  smirks
)

foreach(pytest ${python_tests})
  set_source_files_properties(${pytest}.py PROPERTIES
    PYTHONPATH "${CMAKE_BINARY_DIR}/bindings/python"
    HEDATADIR "${CMAKE_SOURCE_DIR}/data")
  add_python_test(pytest_${pytest} ${pytest}.py)
  set_tests_properties(pytest_${pytest} PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
    ENVIRONMENT "HEDATADIR=${CMAKE_SOURCE_DIR}/data")
endforeach()
