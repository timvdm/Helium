find_package(Boost 1.48 COMPONENTS system iostreams program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(heliumservice
  restfulservice.cpp
  mongoose/mongoose.c)

target_link_libraries(heliumservice
  helium
  ${CMAKE_DL_LIBS}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY} )

target_link_libraries(heliumservice )
