############################################################
#
# The Helium library
#
############################################################

set(Helium_HDRS
  bitvec.h
  contract.h
  concurrent.h
  molecule.h
  substructure.h
  tie.h
  timeout.h
  util.h
  hemol.h
  smiles.h
  distancematrix.h
  ring.h
  # algorithms
  algorithms/bfs.h
  algorithms/canonical.h
  algorithms/components.h
  algorithms/cycles.h
  algorithms/dfs.h
  algorithms/enumeratepaths.h
  algorithms/enumeratesubgraphs.h
  algorithms/extendedconnectivities.h
  algorithms/invariants.h
  algorithms/isomorphism.h
  algorithms/dijkstra.h
  algorithms/floydwarshall.h
  # util
  util/fileio.h
  util/functor.h
  util/string.h
  util/typetraits.h
  util/vector.h
  # fileio
  fileio/file.h
  fileio/fingerprints.h
  fileio/fps.h
  fileio/molecules.h
  # fingerprints
  fingerprints/fingerprints.h
  fingerprints/similarity.h
  # stereo
  stereo/stereounit.h
  stereo/razinger.h
)

# Set up the include directory and symlink/copy all headers to it.
# Note: Windows lacks symlinks and so the files must be copied.
make_directory("${CMAKE_BINARY_DIR}/include/Helium")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/util")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/stereo")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/fileio")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/fingerprints")
make_directory("${CMAKE_BINARY_DIR}/include/Helium/algorithms")
if(UNIX)
  set(header_cmd "create_symlink")
else()
  set(header_cmd "copy_if_different")
endif()
foreach(headerFile ${Helium_HDRS})
  set(from "${CMAKE_CURRENT_SOURCE_DIR}/${headerFile}")
  set(to "${CMAKE_BINARY_DIR}/include/Helium/${headerFile}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E ${header_cmd} "${from}" "${to}")
endforeach()

set(Helium_SRCS
  hemol.cpp
  fileio/file.cpp
  thirdparty/jsoncpp/jsoncpp.cpp
)

# Helium library
add_library(helium STATIC ${Helium_SRCS})
target_link_libraries(helium ${Boost_LIBRARIES})


# install target for header files
set(INCLUDE_INSTALL_DIR  "include" CACHE STRING "Install dir for include files")
install(FILES ${Helium_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium")
install(FILES ${Helium_util_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/util")
install(FILES ${Helium_fileio_HDRS} DESTINATION "${INCLUDE_INSTALL_DIR}/Helium/fileio")
