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
            -P ${CMAKE_SOURCE_DIR}/ctest/removeBpFiles.cmake)
endfunction()

function(mpi_test TESTNAME PROCS EXE)
    removeBpFiles(${TESTNAME}_cleanup)
    add_test(
            NAME ${TESTNAME}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
    )
endfunction(mpi_test)

function(dual_mpi_test)
  set(oneValueArgs TESTNAME TIMEOUT NAME1 NAME2 EXE1 EXE2 PROCS1 PROCS2)
  set(multiValueArgs ARGS1 ARGS2)
  cmake_parse_arguments(DUALTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  removeBpFiles(${DUALTEST_TESTNAME}_cleanup)
  add_test(NAME ${DUALTEST_TESTNAME}
    COMMAND ${CMAKE_SOURCE_DIR}/ctest/runMultipleMpiJobs.sh
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG}
    ${DUALTEST_NAME1} ${DUALTEST_PROCS1} ${DUALTEST_EXE1} "${DUALTEST_ARGS1}"
    ${DUALTEST_NAME2} ${DUALTEST_PROCS2} ${DUALTEST_EXE2} "${DUALTEST_ARGS2}")
  set_tests_properties(${DUALTEST_TESTNAME} PROPERTIES TIMEOUT ${DUALTEST_TIMEOUT})
endfunction(dual_mpi_test)

function(tri_mpi_test)
  set(oneValueArgs TESTNAME TIMEOUT NAME1 NAME2 NAME3 EXE1 EXE2 EXE3 PROCS1 PROCS2 PROCS3)
  set(multiValueArgs ARGS1 ARGS2 ARGS3)
  cmake_parse_arguments(TRITEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  removeBpFiles(${TRITEST_TESTNAME}_cleanup)
  add_test(NAME ${TRITEST_TESTNAME}
    COMMAND ${CMAKE_SOURCE_DIR}/ctest/runMultipleMpiJobs.sh
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG}
    ${TRITEST_NAME1} ${TRITEST_PROCS1} ${TRITEST_EXE1} "${TRITEST_ARGS1}"
    ${TRITEST_NAME2} ${TRITEST_PROCS2} ${TRITEST_EXE2} "${TRITEST_ARGS2}"
    ${TRITEST_NAME3} ${TRITEST_PROCS3} ${TRITEST_EXE3} "${TRITEST_ARGS3}")
  set_tests_properties(${TRITEST_TESTNAME} PROPERTIES TIMEOUT ${TRITEST_TIMEOUT})
endfunction(tri_mpi_test)

#helper variables for readability
set(rendezvous 1)
set(notRendezvous 0)

add_exe(test_ohClassPtn)
set(cyclone1p ${WDMCPL_TEST_DATA_DIR}/cyclone/23elements/mesh.osh/)
dual_mpi_test(TESTNAME test_ohClassPtnSendRecv
  TIMEOUT 4
  NAME1 rdv EXE1 ./test_ohClassPtn PROCS1 2 ARGS1 ${rendezvous} ${cyclone1p}
  NAME2 app EXE2 ./test_ohClassPtn PROCS2 1 ARGS2 ${notRendezvous} ${cyclone1p})

add_exe(test_ohClassPtn_appRibPtn)
set(cyclone2p ${WDMCPL_TEST_DATA_DIR}/cyclone/23elements/2p.osh/)
dual_mpi_test(TESTNAME test_ohClassPtn_appRibPtn
  TIMEOUT 8
  NAME1 rdv EXE1 ./test_ohClassPtn_appRibPtn PROCS1 2 ARGS1 ${rendezvous} ${cyclone1p}
  NAME2 app EXE2 ./test_ohClassPtn_appRibPtn PROCS2 2 ARGS2 ${notRendezvous} ${cyclone2p})

set(d3d1p ${WDMCPL_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc.osh/)
set(d3d4p_cpn ${WDMCPL_TEST_DATA_DIR}/d3d/4p.cpn)
set(d3d16p ${WDMCPL_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc_p16.osh/)
add_exe(test_ohOverlap)
dual_mpi_test(TESTNAME test_ohOverlap_d3d_20p
  TIMEOUT 10
  NAME1 rdv EXE1 ./test_ohOverlap PROCS1 4 ARGS1 ${rendezvous} ${d3d1p} ${d3d4p_cpn}
  NAME2 app EXE2 ./test_ohOverlap PROCS2 16 ARGS2 ${notRendezvous} ${d3d16p} ignored)

set(d3d8p ${WDMCPL_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc_p8.osh/)
add_exe(test_twoClientOverlap)
tri_mpi_test(TESTNAME test_twoClientOverlap_d3d_28p
  TIMEOUT 12
  NAME1 rdv EXE1 ./test_twoClientOverlap PROCS1 4 ARGS1 -1 ${d3d1p} ${d3d4p_cpn}
  NAME2 client0 EXE2 ./test_twoClientOverlap PROCS2 16 ARGS2 0 ${d3d16p} ignored
  NAME3 client1 EXE3 ./test_twoClientOverlap PROCS3 8 ARGS3 1 ${d3d8p} ignored)

add_executable(proxy_coupling XGCProxy.cpp test_proxy_coupling.cpp)
target_link_libraries(proxy_coupling PUBLIC wdmcpl::wdmcpl test_support)
tri_mpi_test(TESTNAME test_proxy_coupling
        TIMEOUT 12
        NAME1 rdv EXE1 ./proxy_coupling PROCS1 4 ARGS1 -1 ${d3d1p} ${d3d4p_cpn}
        NAME2 client0 EXE2 ./proxy_coupling PROCS2 16 ARGS2 0 ${d3d16p} ignored
        NAME3 client1 EXE3 ./proxy_coupling PROCS3 8 ARGS3 1 ${d3d8p} ignored)

# unit tests
find_package(Catch2 2 REQUIRED)
message(STATUS "Found Catch2: ${Catch2_DIR} (found version ${Catch2_VERSION})")
add_executable(unit_tests unit_test_main.cpp test_coordinate_transform.cpp test_coordinate.cpp test_point_search.cpp test_bounding_box.cpp test_uniform_grid.cpp test_field_transfer.cpp)
target_link_libraries(unit_tests PUBLIC Catch2::Catch2 wdmcpl::wdmcpl)
target_include_directories(unit_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(Catch)
catch_discover_tests(unit_tests)
