if(PCMS_ENABLE_OMEGA_H)
  add_library(test_support test/test_support.cpp)
  target_link_libraries(test_support pcms::core) # for omegah and redev
endif()

# TODO use submodule FetchContent/ExternalProject/ExternalData
set(PCMS_TEST_DATA_DIR
    ""
    CACHE PATH "Path to a local copy of the pcms_coupling_data repo.")
# The data is only needed to run tests. Configuring with BUILD_TESTING=OFF
# (e.g. a package build without tests) shouldn't require it.
if(BUILD_TESTING AND NOT EXISTS "${PCMS_TEST_DATA_DIR}")
  message(
    FATAL_ERROR "PCMS_TEST_DATA_DIR \"${PCMS_TEST_DATA_DIR}\" is not accessible"
  )
endif()

set(VALGRIND_EXECUTABLE
    "none"
    CACHE FILEPATH "path to valgrind executable")
set(VALGRIND_ARGS
    "none"
    CACHE
      STRING
      "specify valgrind options; logging (--log-file=%p_<name>.vg) is enabled by default if VALGRIND_EXECUTABLE is set"
)
if(MPIEXEC_PREFLAGS STREQUAL "")
  set(MPIEXEC_PREFLAGS "none")
endif()

set(PCMS_TEST_OMP_NUM_THREADS
    "1"
    CACHE STRING "OMP_NUM_THREADS for each rank of the MPI tests")

function(set_mpi_test_environment TESTNAME)
  set_property(
    TEST ${TESTNAME}
    APPEND
    PROPERTY ENVIRONMENT "OMP_NUM_THREADS=${PCMS_TEST_OMP_NUM_THREADS}")
endfunction()

# Run an MPI test in a directory of its own so that the ADIOS2 files (*.bp,
# *.sst) and logs of concurrently running tests don't collide, and tell ctest
# how many ranks the test uses so that `ctest -j N` treats N as a core budget.
function(isolate_mpi_test TESTNAME NPROCS)
  set(dir ${CMAKE_CURRENT_BINARY_DIR}/mpi_test_runs/${TESTNAME})
  file(MAKE_DIRECTORY ${dir})
  set_tests_properties(${TESTNAME} PROPERTIES WORKING_DIRECTORY ${dir}
                                              PROCESSORS ${NPROCS})
endfunction()

message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_PREFLAGS: ${MPIEXEC_PREFLAGS}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
message(STATUS "VALGRIND_EXECUTABLE: ${VALGRIND_EXECUTABLE}")
message(STATUS "VALGRIND_ARGS: ${VALGRIND_ARGS}")

function(add_exe NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} pcms::core)
  if(PCMS_ENABLE_OMEGA_H)
    target_link_libraries(${NAME} test_support)
  endif()
  target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  if(PCMS_HAS_ASAN)
    target_compile_options(${NAME} PRIVATE -fsanitize=address
                                           -fno-omit-frame-pointer)
    target_link_libraries(${NAME} asan rt)
  endif()
endfunction(add_exe)

function(increaseTimeoutForValgrind TIME TIME_VAR_OUT)
  set(factor 10)
  if(NOT ${VALGRIND_EXECUTABLE} MATCHES "none")
    math(EXPR TIME_OUT "${TIME} * ${factor}")
    set(${TIME_VAR_OUT}
        ${TIME_OUT}
        PARENT_SCOPE)
  endif()
  if(DEFINED PCMS_TIMEOUT)
    set(${TIME_VAR_OUT}
        ${PCMS_TIMEOUT}
        PARENT_SCOPE)
  endif()
endfunction()

# MPI tests don't run in the build directory, so an executable given as a
# target name is replaced with the full path to the target's file. Anything
# else (a path or generator expression) is used as is.
function(resolve_test_exe EXE_VAR)
  if(TARGET ${${EXE_VAR}})
    set(${EXE_VAR}
        $<TARGET_FILE:${${EXE_VAR}}>
        PARENT_SCOPE)
  endif()
endfunction()

function(mpi_test TESTNAME PROCS EXE)
  resolve_test_exe(EXE)
  set(TEST_COMMAND "")
  list(APPEND TEST_COMMAND ${MPIEXEC_EXECUTABLE})
  if(NOT ${MPIEXEC_PREFLAGS} MATCHES "none")
    list(APPEND TEST_COMMAND ${MPIEXEC_PREFLAGS})
  endif()
  list(APPEND TEST_COMMAND ${MPIEXEC_NUMPROC_FLAG} ${PROCS})
  if(NOT ${VALGRIND_EXECUTABLE} MATCHES "none")
    list(APPEND TEST_COMMAND ${VALGRIND_EXECUTABLE} ${VALGRIND_ARGS})
  endif()
  list(APPEND TEST_COMMAND ${EXE} ${ARGN})
  add_test(NAME ${TESTNAME} COMMAND ${TEST_COMMAND})
  set_mpi_test_environment(${TESTNAME})
  isolate_mpi_test(${TESTNAME} ${PROCS})
endfunction(mpi_test)

function(dual_mpi_test)
  set(oneValueArgs
      TESTNAME
      TIMEOUT
      NAME1
      NAME2
      EXE1
      EXE2
      PROCS1
      PROCS2)
  set(multiValueArgs ARGS1 ARGS2)
  cmake_parse_arguments(DUALTEST "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})
  resolve_test_exe(DUALTEST_EXE1)
  resolve_test_exe(DUALTEST_EXE2)
  increasetimeoutforvalgrind(${DUALTEST_TIMEOUT} DUALTEST_TIMEOUT)
  add_test(
    NAME ${DUALTEST_TESTNAME}
    COMMAND
      ${CMAKE_SOURCE_DIR}/ctest/runMultipleMpiJobs.sh --clean-adios2
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${VALGRIND_EXECUTABLE}
      ${VALGRIND_ARGS} ${DUALTEST_NAME1} ${DUALTEST_PROCS1} ${DUALTEST_EXE1}
      "${DUALTEST_ARGS1}" ${DUALTEST_NAME2} ${DUALTEST_PROCS2} ${DUALTEST_EXE2}
      "${DUALTEST_ARGS2}")
  set_tests_properties(${DUALTEST_TESTNAME} PROPERTIES TIMEOUT
                                                       ${DUALTEST_TIMEOUT})
  set_mpi_test_environment(${DUALTEST_TESTNAME})
  math(EXPR nprocs "${DUALTEST_PROCS1} + ${DUALTEST_PROCS2}")
  isolate_mpi_test(${DUALTEST_TESTNAME} ${nprocs})
endfunction(dual_mpi_test)

function(tri_mpi_test)
  set(oneValueArgs
      TESTNAME
      TIMEOUT
      NAME1
      NAME2
      NAME3
      EXE1
      EXE2
      EXE3
      PROCS1
      PROCS2
      PROCS3)
  set(multiValueArgs ARGS1 ARGS2 ARGS3)
  cmake_parse_arguments(TRITEST "${options}" "${oneValueArgs}"
                        "${multiValueArgs}" ${ARGN})
  resolve_test_exe(TRITEST_EXE1)
  resolve_test_exe(TRITEST_EXE2)
  resolve_test_exe(TRITEST_EXE3)
  increasetimeoutforvalgrind(${TRITEST_TIMEOUT} TRITEST_TIMEOUT)
  add_test(
    NAME ${TRITEST_TESTNAME}
    COMMAND
      ${CMAKE_SOURCE_DIR}/ctest/runMultipleMpiJobs.sh --clean-adios2
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${VALGRIND_EXECUTABLE}
      ${VALGRIND_ARGS} ${TRITEST_NAME1} ${TRITEST_PROCS1} ${TRITEST_EXE1}
      "${TRITEST_ARGS1}" ${TRITEST_NAME2} ${TRITEST_PROCS2} ${TRITEST_EXE2}
      "${TRITEST_ARGS2}" ${TRITEST_NAME3} ${TRITEST_PROCS3} ${TRITEST_EXE3}
      "${TRITEST_ARGS3}")
  set_tests_properties(${TRITEST_TESTNAME} PROPERTIES TIMEOUT
                                                      ${TRITEST_TIMEOUT})
  set_mpi_test_environment(${TRITEST_TESTNAME})
  math(EXPR nprocs
       "${TRITEST_PROCS1} + ${TRITEST_PROCS2} + ${TRITEST_PROCS3}")
  isolate_mpi_test(${TRITEST_TESTNAME} ${nprocs})
endfunction(tri_mpi_test)
