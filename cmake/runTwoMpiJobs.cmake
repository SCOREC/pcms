cmake_minimum_required(VERSION 3.15.0...3.21.0)
macro(mpi_test PROCS1 EXE1 ARGS1 PROCS2 EXE2 ARGS2)
  execute_process(
    RESULTS_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
    COMMAND_ECHO STDOUT
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS1}
            ${VALGRIND} ${VALGRIND_ARGS}
            ${EXE1} ${ARGS1}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS2}
            ${VALGRIND} ${VALGRIND_ARGS}
            ${EXE2} ${ARGS2}
  )
  foreach(res IN LISTS result)
    message(STATUS "res ${res}")
    if(res)
      message(FATAL_ERROR "Error stdout ${stdout} stderr ${stderr}")
    endif()
  endforeach()
endmacro(mpi_test)

mpi_test(${PROCS1} ${EXE1} ${ARGS1} ${PROCS2} ${EXE2} ${ARGS2})
