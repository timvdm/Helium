include(UsePythonTest)

set(python_tests
  element
  molecule
  rings
  smarts
  smiles
  smirks
)
    
foreach(pytest ${python_tests})
  add_python_test(pytest_${pytest} ${pytest}.py)
  set_tests_properties(pytest_${pytest} PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
endforeach()
