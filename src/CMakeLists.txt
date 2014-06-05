############################################################
#
# The Helium library
#
############################################################

set(HELIUM_HDRS
  bitvec.h
  contract.h
  concepts.h
  concurrent.h
  element.h
  error.h
  molecule.h
  tie.h
  timeout.h
  util.h
  hemol.h
  smiles.h
  smiley.h
  distancematrix.h
  ring.h
  smarts.h
  smirks.h
  smartmol.h
)

set(HELIUM_ALGORITHMS_HDRS
  algorithms/bfs.h
  algorithms/canonical.h
  algorithms/components.h
  algorithms/cycles.h
  algorithms/dfs.h
  algorithms/enumeratepaths.h
  algorithms/enumeratesubgraphs.h
  algorithms/extendedconnectivities.h
  algorithms/gtd.h
  algorithms/invariants.h
  algorithms/isomorphism.h
  algorithms/dijkstra.h
  algorithms/floydwarshall.h
  algorithms/smartcomponents.h
)

set(HELIUM_UTIL_HDRS
  util/fileio.h
  util/functor.h
  util/string.h
  util/typetraits.h
  util/vector.h
)

set(HELIUM_FILEIO_HDRS
  fileio/file.h
  fileio/fingerprints.h
  fileio/fps.h
  fileio/moleculefile.h
  fileio/cml.h
)

set(HELIUM_FINGERPRINTS_HDRS
  fingerprints/fingerprints.h
  fingerprints/similarity.h
)

set(HELIUM_STEREO_HDRS
  stereo/stereounit.h
  stereo/razinger.h
)

set(HELIUM_TOOLKITS_HDRS
  toolkits/openbabel.h
)

set(HELIUM_CHEMIST_HDRS
  chemist/molecule.h
)

set(HELIUM_JSON_HDRS
  json/json.h
)

set(HELIUM_DEPICT_HDRS
  depict/painter.h
  depict/svgpainter.h
  depict/depict.h
)

set(HELIUM_ALL_HDRS
  ${HELIUM_HDRS}
  ${HELIUM_ALGORITHMS_HDRS}
  ${HELIUM_UTIL_HDRS}
  ${HELIUM_FILEIO_HDRS}
  ${HELIUM_FINGERPRINTS_HDRS}
  ${HELIUM_STEREO_HDRS}
  ${HELIUM_TOOLKITS_HDRS}
  ${HELIUM_CHEMIST_HDRS}
  ${HELIUM_JSON_HDRS}
  ${HELIUM_DEPICT_HDRS}
)

# Set up the include directory and symlink/copy all headers to it.
# Note: Windows lacks symlinks and so the files must be copied.
make_directory("${CMAKE_BINARY_DIR}/include/Helium")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/util")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/stereo")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/fileio")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/fingerprints")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/algorithms")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/toolkits")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/chemist")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/json")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/depict")
if(UNIX)
  set(header_cmd "create_symlink")
else()
  set(header_cmd "copy_if_different")
endif()
foreach(headerFile ${HELIUM_ALL_HDRS})
  set(from "${CMAKE_CURRENT_SOURCE_DIR}/${headerFile}")
  set(to "${CMAKE_BINARY_DIR}/include/Helium/${headerFile}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E ${header_cmd} "${from}" "${to}")
endforeach()

set(HELIUM_SRCS
  hemol.cpp
  fileio/file.cpp
  smarts.cpp
  smirks.cpp
  jsoncpp.cpp
  diagram.cpp
  depict/svgpainter.cpp
  depict/depict.cpp
)

# Helium library
add_library(helium STATIC ${HELIUM_SRCS})
target_link_libraries(helium ${Boost_LIBRARIES})
install(TARGETS helium
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)


# install target for header files
set(INCLUDE_INSTALL_DIR  "include" CACHE STRING "Install dir for include files")
install(FILES ${HELIUM_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium")
install(FILES ${HELIUM_ALGORITHMS_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/algorithms")
install(FILES ${HELIUM_UTIL_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/util")
install(FILES ${HELIUM_FILEIO_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/fileio")
install(FILES ${HELIUM_FINGERPRINTS_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/fingerprints")
install(FILES ${HELIUM_STEREO_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/stereo")
install(FILES ${HELIUM_TOOLKITS_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/toolkits")
install(FILES ${HELIUM_JSON_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/json")
install(FILES ${HELIUM_CHEMIST_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/chemist")
install(FILES ${HELIUM_CHEMIST_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/depict")
install(FILES ${CMAKE_BINARY_DIR}/include/Helium/config.h DESTINATION "${INCLUDE_INSTALL_DIR}/Helium")
