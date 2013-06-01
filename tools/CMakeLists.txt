# define DATADIR for tests that need input files
add_definitions(-DDATADIR="${CMAKE_SOURCE_DIR}/data/")

find_package(JsonRpc REQUIRED)

############################################################
#
# The main 'helium' tool
#
############################################################

set(helium_tool_SRCS
  helium.cpp
  header.cpp
  index.cpp
  fold.cpp
  transpose.cpp
  similarity.cpp
  sort.cpp
  substructure.cpp
)

add_executable(helium_tool ${helium_tool_SRCS})
target_link_libraries(helium_tool ${Helium_LIBRARIES} ${OPENBABEL2_LIBRARIES})
set_target_properties(helium_tool PROPERTIES OUTPUT_NAME helium)
install(TARGETS ${tool}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

############################################################
#
# The 'obhelium' conversion tool (depends on OpenBabel)
#
############################################################

find_package(OpenBabel2)
if (OPENBABEL2_FOUND)
  message(STATUS "OpenBabel found: building obhelium tool")
  include_directories(${OPENBABEL2_INCLUDE_DIR})
  add_executable(obhelium obhelium.cpp)
  target_link_libraries(obhelium ${Helium_LIBRARIES} ${OPENBABEL2_LIBRARIES})
  install(TARGETS obhelium
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif()


############################################################
#
# Misc. tools
#
############################################################

set(tools
  #  opencl
  fps2hel
)

foreach(tool ${tools})
  add_executable(${tool} ${tool}.cpp)
  target_link_libraries(${tool} ${Helium_LIBRARIES})
  install(TARGETS ${tool}
                  RUNTIME DESTINATION bin
                  LIBRARY DESTINATION lib
                  ARCHIVE DESTINATION lib)
endforeach(tool ${tools})
