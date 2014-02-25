CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
include_directories(${zipHMM_SOURCE_DIR})

## pthreads
find_package(Threads REQUIRED)

## BLAS
IF(UNIX AND NOT APPLE)
  set(ATLAS_ROOT "/com/extra/ATLAS/3.9.84")
  include_directories("${ATLAS_ROOT}/include")
  link_directories("${ATLAS_ROOT}/lib")
ELSE(UNIX AND NOT APPLE)
  find_package(BLAS REQUIRED)
  message("-- BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")
  message("-- BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
  link_libraries(${BLAS_LIBRARIES})
ENDIF(UNIX AND NOT APPLE)

## COMPILER FLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -O3 ${CMAKE_THREAD_LIBS_INIT}")
message("-- CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

## C++ library
ADD_LIBRARY(zipHMM 
  SHARED
  forwarder.cpp
  simple_stop_forwarder.cpp
  matrix.cpp 
  seq_io.cpp 
  hmm_io.cpp 
  io_utils.cpp 
  prob_spaces.cpp
  hmm_utils.cpp
  PThreadProcessingDevice.cpp
  performance_description.cpp
  Stage1JobControl.cpp
  posterior_decoding.cpp
  viterbi.cpp
  test_utils.cpp
  calibrate.cpp
  debug.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(zipHMM rt cblas atlas)
endif(UNIX AND NOT APPLE)
set_target_properties(zipHMM PROPERTIES OUTPUT_NAME libzipHMM PREFIX "" SUFFIX .so)

## PYTHON MODULE
find_package(PythonLibs)
if(PYTHONLIBS_FOUND)
  find_package(PythonInterp)
  message(STATUS "Python bindings will be compiled.")
  include_directories(${PYTHON_INCLUDE_DIRS})
  link_libraries(${PYTHON_LIBRARIES})


  add_library(pyZipHMM SHARED pyZipHMM.cpp)
  set_target_properties(pyZipHMM PROPERTIES OUTPUT_NAME libpyZipHMM PREFIX "" SUFFIX .so)
  target_link_libraries(pyZipHMM zipHMM ${PYTHON_LIBRARIES})
  if(UNIX AND NOT APPLE)
    target_link_libraries(pyZipHMM rt cblas atlas)
  endif(UNIX AND NOT APPLE)
  
  add_custom_command(TARGET pyZipHMM POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy libpyZipHMM.so pyZipHMM/libpyZipHMM.so)
else(PYTHONLIBS_FOUND)
  message(STATUS "Python bindings will not be compiled.")
endif(PYTHONLIBS_FOUND)


## R MODULE
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(LibR)
if(LIBR_FOUND)
  message(STATUS "R bindings will be compiled.")
  include_directories(${LIBR_INCLUDE_DIRS})
  message(STATUS "R include: ${LIBR_INCLUDE_DIRS}")
  message(STATUS "R library: ${LIBR_LIBRARIES}")

  add_custom_target(rZipHMM ALL COMMAND ${LIBR_EXECUTABLE} CMD SHLIB rZipHMM.cpp libzipHMM.so -o librZipHMM.so)

  add_custom_command(TARGET rZipHMM POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy rZipHMM.r rZipHMM/R/rZipHMM.r)
  add_custom_command(TARGET rZipHMM POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy librZipHMM.so rZipHMM/librZipHMM.so)
  add_custom_command(TARGET rZipHMM POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy librZipHMM.so rZipHMM/src/librZipHMM.so)
  add_custom_command(TARGET rZipHMM POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy libzipHMM.so rZipHMM/libzipHMM.so)
  add_custom_command(TARGET rZipHMM POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy libzipHMM.so rZipHMM/src/libzipHMM.so)

else(LIBR_FOUND)
  message(STATUS "R bindings will not be compiled.")
endif(LIBR_FOUND)

## EXECUTABLES
add_executable(build_forwarder build_forwarder.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(build_forwarder rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(build_forwarder zipHMM)

add_executable(forward forward.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(forward rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(forward zipHMM)

add_executable(calibrate calibrate_main.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(calibrate rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(calibrate zipHMM)

add_executable(generate_seq generate_seq.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(generate_seq rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(generate_seq zipHMM)

add_executable(generate_hmm generate_hmm.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(generate_hmm rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(generate_hmm zipHMM)

add_executable(generate_biased_hmm generate_biased_hmm.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(generate_biased_hmm rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(generate_biased_hmm zipHMM)

add_executable(generate_01_seq generate_01_seq.cpp)

add_executable(cpp_example cpp_example.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(cpp_example rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(cpp_example zipHMM)

## TESTS
add_executable(test_device_io test_device_io.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(test_device_io rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(test_device_io zipHMM)
add_test(
  NAME test_device_io 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_device_io
  )

add_executable(test_forwarder test_forwarder.cpp)
if(UNIX AND NOT APPLE)
  target_link_libraries(test_forwarder rt cblas atlas)
endif(UNIX AND NOT APPLE)
target_link_libraries(test_forwarder zipHMM)
add_test(
  NAME test_forwarder 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_forwarder
  )

if(PYTHONLIBS_FOUND)
  add_test(pyZipHMM "python" "../zipHMM/python_test.py")
endif(PYTHONLIBS_FOUND)

if(LIBR_FOUND)
  add_test(rZipHMM "${LIBR_EXECUTABLE}" "--slave -f ../zipHMM/test_rZipHMM.r")
endif(LIBR_FOUND)

## INSTALLATION
file(GLOB zipHMM_hpp *.hpp)

install(TARGETS zipHMM LIBRARY DESTINATION lib)
install(FILES ${zipHMM_hpp} DESTINATION include/zipHMM)

install(TARGETS build_forwarder RUNTIME DESTINATION bin)
install(TARGETS forward RUNTIME DESTINATION bin)
install(TARGETS calibrate RUNTIME DESTINATION bin)
install(TARGETS generate_hmm RUNTIME DESTINATION bin)
install(TARGETS generate_seq RUNTIME DESTINATION bin)

if(PYTHONLIBS_FOUND)
  install(CODE "message(STATUS \"installing Python module\")")
  if(PYTHON_PREFIX)
    install(CODE "execute_process(WORKING_DIRECTORY zipHMM COMMAND ${PYTHON_EXECUTABLE} \"setup.py\" \"install\" \"--prefix=${PYTHON_PREFIX}\")")
  else(PYTHON_PREFIX)
    install(CODE "execute_process(WORKING_DIRECTORY zipHMM COMMAND ${PYTHON_EXECUTABLE} \"setup.py\" \"install\")")
  endif(PYTHON_PREFIX)
endif(PYTHONLIBS_FOUND)

if(LIBR_FOUND)
  install(CODE "message(STATUS \"installing R package\")")
  if(R_PREFIX)
    install(CODE "execute_process(WORKING_DIRECTORY zipHMM COMMAND ${LIBR_EXECUTABLE} CMD INSTALL -l ${R_PREFIX} \"rZipHMM\")")
  else(R_PREFIX)
    install(CODE "execute_process(WORKING_DIRECTORY zipHMM COMMAND ${LIBR_EXECUTABLE} CMD INSTALL \"rZipHMM\")")
  endif(R_PREFIX)
endif(LIBR_FOUND)