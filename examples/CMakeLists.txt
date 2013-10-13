set(examples
  dfs1
  dfs2
  bfs1
  )

foreach(example ${examples})
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} ${Helium_LIBRARIES})
endforeach(example ${examples})


