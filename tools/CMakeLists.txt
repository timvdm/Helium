# define DATADIR for tests that need input files
add_definitions(-DDATADIR="${CMAKE_SOURCE_DIR}/data/")

set(tools
  #  opencl
)


foreach(tool ${tools})
  add_executable(${tool} ${tool}.cpp)
  target_link_libraries(${tool} ${Helium_LIBRARIES})
  install(TARGETS ${tool}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib)
endforeach(tool ${tools})

set(helium_tool_SRCS
  helium.cpp
  header.cpp
  index.cpp
  fold.cpp
  similarity.cpp
)
  

add_executable(helium_tool ${helium_tool_SRCS})
target_link_libraries(helium_tool ${Helium_LIBRARIES} ${OPENBABEL2_LIBRARIES})
set_target_properties(helium_tool PROPERTIES OUTPUT_NAME helium)
