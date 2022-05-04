cmake_minimum_required(VERSION 3.15.0...3.21.0)
function(mpi_test)
  set(oneValueArgs EXE1 EXE2 PROCS1 PROCS2)
  set(multiValueArgs ARGS1 ARGS2)
  cmake_parse_arguments(MPITEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  execute_process(
    RESULTS_VARIABLE result
    OUTPUT_VARIABLE stdouterr
    ERROR_VARIABLE stdouterr
    COMMAND_ECHO STDOUT
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPITEST_PROCS1}
            ${VALGRIND} ${VALGRIND_ARGS}
            ${MPITEST_EXE1} ${MPITEST_ARGS1}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPITEST_PROCS2}
            ${VALGRIND} ${VALGRIND_ARGS}
            ${MPITEST_EXE2} ${MPITEST_ARGS2}
  )
  foreach(res IN LISTS result)
    message(STATUS "res ${res}")
    if(res)
      message(STATUS "${stdouterr}")
      message(FATAL_ERROR "ERROR...")
    endif()
  endforeach()
endfunction(mpi_test)

mpi_test(${PROCS1} ${EXE1} ${ARGS1} ${PROCS2} ${EXE2} ${ARGS2})
