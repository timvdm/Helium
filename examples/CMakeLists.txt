set(examples
  dfs1
  dfs2
  bfs1
  components
  molecule1
  molecule2
  molecule3
  molecule4
  moleculefile
  substructure
  smarts
  smirks
  cycles
  enumeratesubgraphs
  smartmol
  )

foreach(example ${examples})
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} ${Helium_LIBRARIES})
endforeach(example ${examples})


find_package(OpenBabel2)
if (OPENBABEL2_FOUND)
  message(STATUS "OpenBabel found: building openbabel example")
  include_directories(${OPENBABEL2_INCLUDE_DIR})
  add_executable(opsmirks openbabel.cpp)
  target_link_libraries(opsmirks ${Helium_LIBRARIES} ${OPENBABEL2_LIBRARIES})
endif()

