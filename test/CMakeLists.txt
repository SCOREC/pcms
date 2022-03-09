if(BUILD_TESTING)
  set(WDMCPL_TEST_DATA_DIR "" CACHE PATH
    "Path to a local copy of the wdmapp_coupling_data repo.")
  if(NOT EXISTS ${WDMCPL_TEST_DATA_DIR})
    message(FATAL_ERROR "WDMCPL_TEST_DATA_DIR \"${WDMCPL_TEST_DATA_DIR}\" is not accessible")
  endif()

  message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
  message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")

  function(add_exe NAME)
    add_executable(${NAME} ${NAME}.cpp)
    target_link_libraries(${NAME} wdmcpl)
    if(WDMCPL_HAS_ASAN)
      target_compile_options(${NAME} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
      target_link_libraries(${NAME} asan rt)
    endif()
  endfunction(add_exe)

  function(removeBpFiles TESTNAME)
    add_test(NAME ${TESTNAME} COMMAND ${CMAKE_COMMAND}
      -P ${CMAKE_SOURCE_DIR}/cmake/removeBpFiles.cmake)
  endfunction()

  function(mpi_test TESTNAME PROCS EXE)
    removeBpFiles(${TESTNAME}_cleanup)
    add_test(
      NAME ${TESTNAME}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
    )
  endfunction(mpi_test)

  function(dual_mpi_test TESTNAME PROCS1 EXE1 ARGS1 PROCS2 EXE2 ARGS2)
    removeBpFiles(${TESTNAME}_cleanup)
    add_test(NAME ${TESTNAME}
      COMMAND ${CMAKE_COMMAND}
      -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DPROCS1=${PROCS1} -DEXE1=${EXE1} -DARGS1=${ARGS1}
      -DPROCS2=${PROCS2} -DEXE2=${EXE2} -DARGS2=${ARGS2}
      -P ${CMAKE_SOURCE_DIR}/cmake/runTwoMpiJobs.cmake)
  endfunction(dual_mpi_test)

  add_exe(test_ohClassPtn)
  set(testMesh ${WDMCPL_TEST_DATA_DIR}/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh/)
  dual_mpi_test(test_ohClassPtnSendRecv 1 ./test_ohClassPtn 0 ${testMesh}
                                        2 ./test_ohClassPtn 1 ${testMesh})
endif()
