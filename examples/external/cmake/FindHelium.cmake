# - Try to find Helium lib
# Once done this will define
#
#  HELIUM_FOUND - system has helium lib
#  HELIUM_INCLUDE_DIR - the helium include directory
#  HELIUM_LIBRARIES - the helium library

# Copyright (c) 2013 Tim Vandermeersch <tim.vandermeersch@gamil.com>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if(HELIUM_INCLUDE_DIR AND HELIUM_LIBRARIES)
  # in cache already
  set(HELIUM_FOUND TRUE)
else()
  find_path(HELIUM_INCLUDE_DIR NAMES Helium/hemol.h PATHS /usr/include /usr/local/include)
  find_library(HELIUM_LIBRARIES NAMES helium Helium)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(HELIUM DEFAULT_MSG HELIUM_LIBRARIES
    HELIUM_INCLUDE_DIR)
  mark_as_advanced(HELIUM_INCLUDE_DIR HELIUM_LIBRARIES)
endif()

