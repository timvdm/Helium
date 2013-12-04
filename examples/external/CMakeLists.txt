set(APP example)
set(APP_SRCS
  main.cpp
)

############################################################

cmake_minimum_required(VERSION 2.6)
project(${APP})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(Helium REQUIRED)
include_directories(${HELIUM_INCLUDE_DIR})

add_executable(${APP} ${APP_SRCS})
target_link_libraries(${APP} ${HELIUM_LIBRARIES})
