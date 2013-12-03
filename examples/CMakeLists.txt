set(examples
  dfs1
  dfs2
  bfs1
  components
  molecule1
  molecule2
  molecule3
  molecule4
  substructure
  smirks

  enumeratesubgraphs
  )

foreach(example ${examples})
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} ${Helium_LIBRARIES})
endforeach(example ${examples})


