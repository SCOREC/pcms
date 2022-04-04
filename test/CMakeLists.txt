# TODO use submodule FetchContent/ExternalProject/ExternalData
set(WDMCPL_TEST_DATA_DIR "" CACHE PATH
        "Path to a local copy of the wdmapp_coupling_data repo.")
if (NOT EXISTS ${WDMCPL_TEST_DATA_DIR})
    message(FATAL_ERROR "WDMCPL_TEST_DATA_DIR \"${WDMCPL_TEST_DATA_DIR}\" is not accessible")
endif ()

message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")

add_library(test_support test_support.cpp)
target_link_libraries(test_support wdmcpl) #for omegah and redev

function(add_exe NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} wdmcpl test_support)
  target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
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

function(dual_mpi_test)
  set(oneValueArgs NAME EXE1 EXE2 PROCS1 PROCS2 )
  set(multiValueArgs ARGS1 ARGS2)
  cmake_parse_arguments(DUALTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  removeBpFiles(${DUALTEST_NAME}_cleanup)
  add_test(NAME ${DUALTEST_NAME}
    COMMAND ${CMAKE_COMMAND}
    -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DMPITEST_PROCS1=${DUALTEST_PROCS1} -DMPITEST_EXE1=${DUALTEST_EXE1} "-DMPITEST_ARGS1=${DUALTEST_ARGS1}"
    -DMPITEST_PROCS2=${DUALTEST_PROCS2} -DMPITEST_EXE2=${DUALTEST_EXE2} "-DMPITEST_ARGS2=${DUALTEST_ARGS2}"
    -P ${CMAKE_SOURCE_DIR}/cmake/runTwoMpiJobs.cmake)
endfunction(dual_mpi_test)

add_exe(test_ohClassPtn)
set(testMesh ${WDMCPL_TEST_DATA_DIR}/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh/)
message(STATUS "appArgs ${appArgs}")
set(rendezvous 1)
set(notRendezvous 0)
dual_mpi_test(NAME test_ohClassPtnSendRecv
  EXE1 ./test_ohClassPtn PROCS1 1 ARGS1 ${notRendezvous} ${testMesh}
  EXE2 ./test_ohClassPtn PROCS2 2 ARGS2 ${rendezvous} ${testMesh})

add_exe(test_ohClassPtn_appRibPtn)
set(testMesh ${WDMCPL_TEST_DATA_DIR}/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh/)
set(testMesh2p ${WDMCPL_TEST_DATA_DIR}/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/2p.osh)
message(STATUS "appArgs ${appArgs}")
set(rendezvous 1)
set(notRendezvous 0)
dual_mpi_test(NAME test_ohClassPtn_appRibPtn
  EXE1 ./test_ohClassPtn_appRibPtn PROCS1 2 ARGS1 ${rendezvous} ${testMesh}
  EXE2 ./test_ohClassPtn_appRibPtn PROCS2 2 ARGS2 ${notRendezvous} ${testMesh2p})

add_exe(test_meshRcbRdv)
mpi_test(test_meshRcbRdv_1p 1 ./test_meshRcbRdv
        ${WDMCPL_TEST_DATA_DIR}/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh/)


# unit tests
find_package(Catch2 REQUIRED)
add_executable(unit_tests unit_test_main.cpp test_coordinate_transform.cpp test_coordinate.cpp)
target_link_libraries(unit_tests PUBLIC Catch2::Catch2 wdmcpl::wdmcpl)
include(Catch)
catch_discover_tests(unit_tests)
