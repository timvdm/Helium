# - Find jsonrpc
# Once done this will define
#
#  JSONRPC_INCLUDE_DIRS   - where to find jsonrpc.h, etc.
#  JSONRPC_LIBRARIES      - List of libraries when using jsonrpc.
#  JSONRPC_FOUND          - True if jsonrpc found.
#
# An includer may set JSONRPC_ROOT to a jsonrpc installation root to tell
# this module where to look.

#=============================================================================
# Copyright 2001-2011 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

set(_JSONRPC_SEARCHES)

# Search JSONRPC_ROOT first if it is set.
if(JSONRPC_ROOT)
  set(_JSONRPC_SEARCH_ROOT PATHS ${JSONRPC_ROOT} NO_DEFAULT_PATH)
  list(APPEND _JSONRPC_SEARCHES _JSONRPC_SEARCH_ROOT)
endif()

foreach(search ${_JSONRPC_SEARCHES})
  find_path(JSONRPC_INCLUDE_DIR NAMES jsonrpc/rpc.h ${${search}} PATH_SUFFIXES include)
  find_library(JSON_LIBRARY  NAMES json ${${search}} PATH_SUFFIXES lib)
  find_library(JSONRPC_LIBRARY  NAMES jsonrpc ${${search}} PATH_SUFFIXES lib)
  find_library(MONGOOSE_LIBRARY  NAMES mongoose ${${search}} PATH_SUFFIXES lib)
endforeach()

mark_as_advanced(JSONRPC_LIBRARY JSONRPC_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set JSONRPC_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(JSONRPC REQUIRED_VARS JSONRPC_LIBRARY JSONRPC_INCLUDE_DIR)

if(JSONRPC_FOUND)
    set(JSONRPC_INCLUDE_DIRS ${JSONRPC_INCLUDE_DIR})
    set(JSONRPC_LIBRARIES ${JSONRPC_LIBRARY})
    list(APPEND JSONRPC_LIBRARIES ${JSON_LIBRARY})
    list(APPEND JSONRPC_LIBRARIES ${MONGOOSE_LIBRARY})
endif()

