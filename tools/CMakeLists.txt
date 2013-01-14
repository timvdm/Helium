# define DATADIR for tests that need input files
add_definitions(-DDATADIR="${CMAKE_SOURCE_DIR}/data/")

set(tools
  convert
  index
  index2invindex
  benchmark
  query
  filter
  )


foreach(tool ${tools})
  add_executable(${tool} ${tool}.cpp)
  target_link_libraries(${tool} ${Helium_LIBRARIES})
  #  install(TARGETS ${tool}
  #                  RUNTIME DESTINATION bin
  #                  LIBRARY DESTINATION lib
  #                  ARCHIVE DESTINATION lib)
endforeach(tool ${tools})


