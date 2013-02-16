set(CTEST_SOURCE_DIRECTORY "/home/tim/Cheminformatics/Helium")
set(CTEST_BINARY_DIRECTORY "/home/tim/Cheminformatics/Helium/build-cdash")

set(CTEST_UPDATE_COMMAND "/usr/bin/git")
set(CTEST_CVS_UPDATE_OPTIONS "pull")
#set(CTEST_CVS_CHECKOUT "${CTEST_CVS_COMMAND} clone git://github.com/cryos/avogadro.git \"${CTEST_SOURCE_DIRECTORY}\"")

# which ctest command to use for running the dashboard
set(CTEST_COMMAND "/usr/bin/ctest -D Nightly")

# what cmake command to use for configuring this dashboard
set(CTEST_CMAKE_COMMAND "/usr/bin/cmake")

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# should ctest wipe the binary tree before running
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
set(CTEST_INITIAL_CACHE "CMAKE_CXX_COMPILER:STRING=g++
ENABLE_CPP11:BOOL=True
BUILDNAME:STRING=g++ on Ubuntu linux (C++11)
SITE:STRING=MolDB.net.timvdm
CVSCOMMAND:FILEPATH=/usr/bin/git
MAKECOMMAND:STRING=make -j5
")

# set any extra envionment varibles here
set(CTEST_ENVIRONMENT)
