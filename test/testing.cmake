function(mpi_test TESTNAME PROCS EXE)
  message("adding test ${TESTNAME}")
  add_test(
    NAME ${TESTNAME}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS} ${VALGRIND}
    ${VALGRIND_ARGS} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXE} ${ARGN}
    )
endfunction(mpi_test)

mpi_test(t0 4 test_init ${TEST_DATA_DIR}/testdatas/)
