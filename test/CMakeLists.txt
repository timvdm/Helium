# define DATADIR for tests that need input files
add_definitions(-DDATADIR="${CMAKE_SOURCE_DIR}/data/")

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
  )

foreach(test ${tests})
  add_executable(test_${test} ${test}.cpp)
  target_link_libraries(test_${test}  ${Helium_LIBRARIES} ${OPENBABEL2_LIBRARIES})
  add_test(${test}_Test ${TEST_PATH}/test_${test})
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "FAIL")
endforeach(test ${tests})


