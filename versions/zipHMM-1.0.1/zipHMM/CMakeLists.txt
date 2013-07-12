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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -O3 -g ${CMAKE_THREAD_LIBS_INIT}") # remove -g !!
message("-- CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

## C++ library
file(GLOB zipHMM_hpp *.hpp)
file(GLOB zipHMM_cpp *.cpp)
set(zipHMM_lib_content ${zipHMM_cpp} ${zipHMM_hpp})

ADD_LIBRARY(zipHMM ${zipHMM_lib_content})
if(UNIX AND NOT APPLE)
  target_link_libraries(zipHMM rt cblas atlas)
endif(UNIX AND NOT APPLE)

## PYTHON MODULE
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(pyZipHMM
  SHARED 
  pyZipHMM.cpp 
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
  debug.cpp
  )
set_target_properties(pyZipHMM PROPERTIES OUTPUT_NAME libpyZipHMM PREFIX "" SUFFIX .so)
target_link_libraries(pyZipHMM ${PYTHON_LIBRARIES})
if(UNIX AND NOT APPLE)
  target_link_libraries(pyZipHMM rt cblas atlas)
endif(UNIX AND NOT APPLE)

message("-- PYTHON_SITE_PACKAGES: ${PYTHON_SITE_PACKAGES}")

## EXECUTABLES
add_executable(build_forwarder build_forwarder.cpp 
  forwarder.cpp
  hmm_utils.cpp
  PThreadProcessingDevice.cpp
  performance_description.cpp
  Stage1JobControl.cpp
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
  prob_spaces.cpp
  debug.cpp
  )
if(UNIX AND NOT APPLE)
  target_link_libraries(build_forwarder cblas atlas)
endif(UNIX AND NOT APPLE)

add_executable(forward forward.cpp 
  forwarder.cpp
  hmm_utils.cpp
  PThreadProcessingDevice.cpp
  performance_description.cpp
  Stage1JobControl.cpp
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
  prob_spaces.cpp
  debug.cpp
  )
if(UNIX AND NOT APPLE)
  target_link_libraries(forward cblas atlas)
endif(UNIX AND NOT APPLE)

# add_executable(viterbi viterbi_main.cpp 
#   viterbi.cpp
#   PThreadProcessingDevice.cpp
#   performance_description.cpp
#   matrix.cpp
#   seq_io.cpp
#   hmm_io.cpp
#   hmm_utils.cpp
#   io_utils.cpp
#   prob_spaces.cpp
#   debug.cpp
#   )
# if(UNIX AND NOT APPLE)
#   target_link_libraries(forward cblas atlas)
# endif(UNIX AND NOT APPLE)

add_executable(calibrate calibrate.cpp 
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
  prob_spaces.cpp
  forwarder.cpp
  hmm_utils.cpp
  PThreadProcessingDevice.cpp
  performance_description.cpp
  Stage1JobControl.cpp
  debug.cpp
  )
if(UNIX AND NOT APPLE)
  target_link_libraries(calibrate rt cblas atlas)
endif(UNIX AND NOT APPLE)

add_executable(generate_seq generate_seq.cpp 
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
)

add_executable(generate_hmm generate_hmm.cpp 
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
)

add_executable(generate_biased_hmm generate_biased_hmm.cpp 
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
)

add_executable(generate_01_seq generate_01_seq.cpp)

add_executable(cpp_example cpp_example.cpp
  matrix.cpp
  seq_io.cpp
  hmm_io.cpp
  io_utils.cpp
  prob_spaces.cpp
  forwarder.cpp
  hmm_utils.cpp
  PThreadProcessingDevice.cpp
  performance_description.cpp
  Stage1JobControl.cpp
  posterior_decoding.cpp
  viterbi.cpp
  debug.cpp
)
if(UNIX AND NOT APPLE)
  target_link_libraries(cpp_example cblas atlas)
endif(UNIX AND NOT APPLE)

## TESTS
add_executable(test_device_io test_device_io.cpp 
  performance_description.cpp
  io_utils.cpp
  prob_spaces.cpp
  hmm_utils.cpp
  PThreadProcessingDevice.cpp
  Stage1JobControl.cpp
  test_utils.cpp
  debug.cpp
  )
if(UNIX AND NOT APPLE)
  target_link_libraries(test_device_io cblas atlas)
endif(UNIX AND NOT APPLE)

add_test(
  NAME test_device_io 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_device_io
  )

add_executable(test_forwarder test_forwarder.cpp
  hmm_io.cpp
  seq_io.cpp
  io_utils.cpp 
  forwarder.cpp
  simple_stop_forwarder.cpp
  hmm_utils.cpp
  prob_spaces.cpp
  test_utils.cpp
  PThreadProcessingDevice.cpp
  performance_description.cpp
  Stage1JobControl.cpp
  debug.cpp
  )
if(UNIX AND NOT APPLE)
  target_link_libraries(test_forwarder rt cblas atlas)
endif(UNIX AND NOT APPLE)

add_test(
  NAME test_forwarder 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_forwarder
  )

add_test(pyZipHMM "python" "../zipHMM/python_test.py")


## INSTALLATION
set(zipHMM ${zipHMM_hpp} ${zipHMM_cpp})

install(TARGETS zipHMM ARCHIVE DESTINATION lib)
install(FILES ${zipHMM_hpp} DESTINATION include/zipHMM)

execute_process(COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
install(TARGETS pyZipHMM LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES})
install(FILES pyZipHMM.py DESTINATION ${PYTHON_SITE_PACKAGES})

install(TARGETS build_forwarder RUNTIME DESTINATION bin)
install(TARGETS forward RUNTIME DESTINATION bin)
install(TARGETS calibrate RUNTIME DESTINATION bin)
install(TARGETS generate_hmm RUNTIME DESTINATION bin)
install(TARGETS generate_seq RUNTIME DESTINATION bin)
