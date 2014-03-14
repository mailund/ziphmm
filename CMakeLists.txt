cmake_minimum_required(VERSION 2.8)

project(zipHMM CXX C)

enable_testing()

add_subdirectory(zipHMM)

# Packager
set(CPACK_PACKAGE_NAME "zipHMM")
set(CPACK_PACKAGE_VENDOR "BiRC - Bioinformatics Research Center, Aarhus University, Denmark")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(CPACK_SYSTEM_NAME x86_64)
  set(CPACK_TOPLEVEL_TAG x86_64)
else (CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(CPACK_SYSTEM_NAME x86_32)
  set(CPACK_TOPLEVEL_TAG x86_32)
endif (CMAKE_SIZEOF_VOID_P MATCHES "8")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "zipHMM: ")
set(CPACK_RESOURCE_FILE_LICENSE "COPYING")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
  "\\\\./Makefile$"
  "zipHMM/Makefile$"
  "install_manifest\\\\.txt$"
  "CMakeFiles"
  "CTestTestfile\\\\.cmake$"
  "cmake_install\\\\.cmake$"
  "_CPack_Packages"
  "CPackSourceConfig\\\\.cmake$"
  "CPackConfig\\\\.cmake$"
  "CMakeCache\\\\.txt$"
  ""
  "\\\\.svn"
  "\\\\.gz$"
  "\\\\.zip$"
  "\\\\.pyc$"
  "\\\\.a$"
  "\\\\.so$"
  "\\\\.o"
  "~$"
  "rZipHMM\\\\.Rcheck"
  "testtest\\\\.hmm$"
  "tmp.seq$"
  "example_out$"
  "Testing"
  "bin"
  "web"

  ${CPACK_SOURCE_IGNORE_FILES}
  )

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
