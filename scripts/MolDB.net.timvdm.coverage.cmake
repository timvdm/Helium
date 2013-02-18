cmake_minimum_required(VERSION 2.8)

set(CTEST_SITE "MolDB.net.timvdm")
set(CTEST_BUILD_NAME "g++ on Ubuntu linux (Coverage)")

set(CTEST_TEST_TIMEOUT 300)

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}")

set(source_dir "Helium")
set(build_dir "Helium-build-coverage")
set(CTEST_DASHBOARD_ROOT "$ENV{HOME}/Cheminformatics/cdash")
set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${source_dir}")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${build_dir}")

# Use launchers to get better compiler erros/warnings
set(CTEST_USE_LAUNCHERS 1)

set(CTEST_BUILD_FLAGS "-j9")
set(CTEST_COVERAGE_COMMAND "/usr/bin/gcov")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

set(common_flags "-fdiagnostics-show-option -Wall -Wextra -Wpointer-arith -Winvalid-pch -Wcast-align -Wdisabled-optimization -Wwrite-strings -fstack-protector-all -D_FORTIFY_SOURCE=2 -Wconversion -Wno-error=sign-conversion -Wno-error=conversion -Werror=return-type")

set(cov_options "-fprofile-arcs -ftest-coverage")

# Write initial cache.
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
CMAKE_CXX_FLAGS:STRING=${common_flags} ${cov_options} -Woverloaded-virtual -Wstrict-null-sentinel -pipe
CMAKE_C_FLAGS:STRING=${common_flags} ${cov_options} -pipe
CMAKE_EXE_LINKER_FLAGS:STRING=${cov_options}
CMAKE_SHARED_LINKER_FLAGS:STRING=${cov_options}
CMAKE_MODULE_LINKER_FLAGS:STRING=${cov_options}
CMAKE_BUILD_TYPE:STRING=Debug
CTEST_USE_LAUNCHERS:BOOL=ON
")

set(CTEST_UPDATE_COMMAND "git")

set(CTEST_CUSTOM_COVERAGE_EXCLUDE
  ${CTEST_CUSTOM_COVERAGE_EXCLUDE} # keep current exclude expressions
  "/thirdparty/"
  "/test/"
)

#ctest_start(Experimental)
ctest_start(Nightly)
ctest_update()
ctest_configure()
ctest_submit(PARTS Update Configure Notes)
ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
ctest_build(APPEND)
ctest_submit(PARTS Build)
ctest_test(APPEND)
ctest_submit(PARTS Test)

# Perform coverage on this build, and submit that too
ctest_coverage(APPEND)
ctest_submit(PARTS Coverage)
