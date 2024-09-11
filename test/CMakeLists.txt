# TODO use submodule FetchContent/ExternalProject/ExternalData
set(PCMS_TEST_DATA_DIR "" CACHE PATH
        "Path to a local copy of the pcms_coupling_data repo.")
if (NOT EXISTS ${PCMS_TEST_DATA_DIR})
    message(FATAL_ERROR "PCMS_TEST_DATA_DIR \"${PCMS_TEST_DATA_DIR}\" is not accessible")
endif ()

include(ProcessorCount)
ProcessorCount(HOST_NPROC)
message("TESTING ON SYSTEM WITH ${HOST_NPROC} CPU Cores")

message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
set(VALGRIND_EXECUTABLE "none" CACHE FILEPATH "path to valgrind executable")
set(VALGRIND_ARGS "none" CACHE STRING "specify valgrind options; logging (--log-file=%p_<name>.vg) is enabled by default if VALGRIND_EXECUTABLE is set")
message(STATUS "VALGRIND_EXECUTABLE: ${VALGRIND_EXECUTABLE}")
message(STATUS "VALGRIND_ARGS: ${VALGRIND_ARGS}")

if (PCMS_ENABLE_OMEGA_H)
    add_library(test_support test_support.cpp)
    target_link_libraries(test_support pcms::core ) #for omegah and redev
endif ()

function(add_exe NAME)
    add_executable(${NAME} ${NAME}.cpp)
    target_link_libraries(${NAME} pcms::core)
    if (PCMS_ENABLE_OMEGA_H)
        target_link_libraries(${NAME} test_support)
    endif ()
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    if (PCMS_HAS_ASAN)
        target_compile_options(${NAME} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
        target_link_libraries(${NAME} asan rt)
    endif ()
endfunction(add_exe)

function(removeBpFiles TESTNAME)
    add_test(NAME ${TESTNAME} COMMAND ${CMAKE_COMMAND}
            -P ${CMAKE_SOURCE_DIR}/ctest/removeBpFiles.cmake)
endfunction()

function(increaseTimeoutForValgrind TIME TIME_VAR_OUT)
    set(factor 10)
    if (NOT ${VALGRIND_EXECUTABLE} MATCHES "none")
        math(EXPR TIME_OUT "${TIME} * ${factor}")
        set(${TIME_VAR_OUT} ${TIME_OUT} PARENT_SCOPE)
    endif ()
endfunction()

function(mpi_test TESTNAME PROCS EXE)
    removeBpFiles(${TESTNAME}_cleanup)
    if (${VALGRIND_EXECUTABLE} MATCHES "none")
        add_test(NAME ${TESTNAME}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN})
    else ()
        add_test(NAME ${TESTNAME}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS}
                ${VALGRIND_EXECUTABLE} ${VALGRIND_ARGS} ${EXE} ${ARGN})
    endif ()
endfunction(mpi_test)

function(dual_mpi_test)
    set(oneValueArgs TESTNAME TIMEOUT NAME1 NAME2 EXE1 EXE2 PROCS1 PROCS2)
    set(multiValueArgs ARGS1 ARGS2)
    cmake_parse_arguments(DUALTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    removeBpFiles(${DUALTEST_TESTNAME}_cleanup)
    increaseTimeoutForValgrind(${DUALTEST_TIMEOUT} DUALTEST_TIMEOUT)
    add_test(NAME ${DUALTEST_TESTNAME}
            COMMAND ${CMAKE_SOURCE_DIR}/ctest/runMultipleMpiJobs.sh
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG}
            ${VALGRIND_EXECUTABLE} ${VALGRIND_ARGS}
            ${DUALTEST_NAME1} ${DUALTEST_PROCS1} ${DUALTEST_EXE1} "${DUALTEST_ARGS1}"
            ${DUALTEST_NAME2} ${DUALTEST_PROCS2} ${DUALTEST_EXE2} "${DUALTEST_ARGS2}")
    set_tests_properties(${DUALTEST_TESTNAME} PROPERTIES TIMEOUT ${DUALTEST_TIMEOUT})
endfunction(dual_mpi_test)

function(tri_mpi_test)
    set(oneValueArgs TESTNAME TIMEOUT NAME1 NAME2 NAME3 EXE1 EXE2 EXE3 PROCS1 PROCS2 PROCS3)
    set(multiValueArgs ARGS1 ARGS2 ARGS3)
    cmake_parse_arguments(TRITEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    removeBpFiles(${TRITEST_TESTNAME}_cleanup)
    increaseTimeoutForValgrind(${TRITEST_TIMEOUT} TRITEST_TIMEOUT)
    add_test(NAME ${TRITEST_TESTNAME}
            COMMAND ${CMAKE_SOURCE_DIR}/ctest/runMultipleMpiJobs.sh
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG}
            ${VALGRIND_EXECUTABLE} ${VALGRIND_ARGS}
            ${TRITEST_NAME1} ${TRITEST_PROCS1} ${TRITEST_EXE1} "${TRITEST_ARGS1}"
            ${TRITEST_NAME2} ${TRITEST_PROCS2} ${TRITEST_EXE2} "${TRITEST_ARGS2}"
            ${TRITEST_NAME3} ${TRITEST_PROCS3} ${TRITEST_EXE3} "${TRITEST_ARGS3}")
    set_tests_properties(${TRITEST_TESTNAME} PROPERTIES TIMEOUT ${TRITEST_TIMEOUT})
endfunction(tri_mpi_test)

#helper variables for readability
set(rendezvous 1)
set(notRendezvous 0)

if (PCMS_ENABLE_OMEGA_H)
    add_exe(test_ohClassPtn)
    set(cyclone1p ${PCMS_TEST_DATA_DIR}/cyclone/23elements/mesh.osh/)
    dual_mpi_test(TESTNAME test_ohClassPtnSendRecv
            TIMEOUT 4
            NAME1 rdv EXE1 ./test_ohClassPtn PROCS1 2 ARGS1 ${rendezvous} ${cyclone1p}
            NAME2 app EXE2 ./test_ohClassPtn PROCS2 1 ARGS2 ${notRendezvous} ${cyclone1p})

    add_exe(test_ohClassPtn_appRibPtn)
    set(cyclone2p ${PCMS_TEST_DATA_DIR}/cyclone/23elements/2p.osh/)
    dual_mpi_test(TESTNAME test_ohClassPtn_appRibPtn
            TIMEOUT 8
            NAME1 rdv EXE1 ./test_ohClassPtn_appRibPtn PROCS1 2 ARGS1 ${rendezvous} ${cyclone1p}
            NAME2 app EXE2 ./test_ohClassPtn_appRibPtn PROCS2 2 ARGS2 ${notRendezvous} ${cyclone2p})

    set(d3d1p ${PCMS_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc.osh/)
    set(d3d1p_cpn ${PCMS_TEST_DATA_DIR}/d3d/1p.cpn)
    set(d3d2p_cpn ${PCMS_TEST_DATA_DIR}/d3d/2p.cpn)
    set(d3d4p_cpn ${PCMS_TEST_DATA_DIR}/d3d/4p.cpn)
    set(d3d8p_cpn ${PCMS_TEST_DATA_DIR}/d3d/8p.cpn)
    set(d3d16p_cpn ${PCMS_TEST_DATA_DIR}/d3d/16p.cpn)
    set(d3d16p ${PCMS_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc_p16.osh/)
    add_exe(test_ohOverlap)
    if(HOST_NPROC GREATER 10)
        dual_mpi_test(TESTNAME test_ohOverlap_d3d_20p
                TIMEOUT 10
                NAME1 rdv EXE1 ./test_ohOverlap PROCS1 4 ARGS1 ${rendezvous} ${d3d1p} ${d3d4p_cpn}
                NAME2 app EXE2 ./test_ohOverlap PROCS2 16 ARGS2 ${notRendezvous} ${d3d16p} ignored)
    endif()

    set(d3d8p ${PCMS_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc_p8.osh/)
    add_exe(test_twoClientOverlap)
    if(HOST_NPROC GREATER 14)
        tri_mpi_test(TESTNAME test_twoClientOverlap_d3d_28p
                TIMEOUT 12
                NAME1 rdv EXE1 ./test_twoClientOverlap PROCS1 4 ARGS1 -1 ${d3d1p} ${d3d4p_cpn}
                NAME2 client0 EXE2 ./test_twoClientOverlap PROCS2 16 ARGS2 0 ${d3d16p} ignored
                NAME3 client1 EXE3 ./test_twoClientOverlap PROCS3 8 ARGS3 1 ${d3d8p} ignored)
    endif()
    add_executable(proxy_coupling test_proxy_coupling.cpp)
    target_link_libraries(proxy_coupling PUBLIC pcms::core test_support)
    tri_mpi_test(TESTNAME test_proxy_coupling_4p
            TIMEOUT 10
            NAME1 rdv EXE1 ./proxy_coupling PROCS1 2 ARGS1 -1 ${d3d1p} ${d3d2p_cpn}
            NAME2 client0 EXE2 ./proxy_coupling PROCS2 1 ARGS2 0 ${d3d1p} ignored
            NAME3 client1 EXE3 ./proxy_coupling PROCS3 1 ARGS3 1 ${d3d1p} ignored)
    if(HOST_NPROC GREATER 13)
        tri_mpi_test(TESTNAME test_proxy_coupling_26p
                TIMEOUT 12
                NAME1 rdv EXE1 ./proxy_coupling PROCS1 2 ARGS1 -1 ${d3d1p} ${d3d2p_cpn}
                NAME2 client0 EXE2 ./proxy_coupling PROCS2 16 ARGS2 0 ${d3d16p} ignored
                NAME3 client1 EXE3 ./proxy_coupling PROCS3 8 ARGS3 1 ${d3d8p} ignored)
    endif()
    if (HOST_NPROC GREATER 14)
        tri_mpi_test(TESTNAME test_proxy_coupling_28p
                TIMEOUT 12
                NAME1 rdv EXE1 ./proxy_coupling PROCS1 4 ARGS1 -1 ${d3d1p} ${d3d4p_cpn}
                NAME2 client0 EXE2 ./proxy_coupling PROCS2 16 ARGS2 0 ${d3d16p} ignored
                NAME3 client1 EXE3 ./proxy_coupling PROCS3 8 ARGS3 1 ${d3d8p} ignored)
    endif ()
    if(HOST_NPROC GREATER 16)
        tri_mpi_test(TESTNAME test_proxy_coupling_32p
                TIMEOUT 12
                NAME1 rdv EXE1 ./proxy_coupling PROCS1 8 ARGS1 -1 ${d3d1p} ${d3d8p_cpn}
                NAME2 client0 EXE2 ./proxy_coupling PROCS2 16 ARGS2 0 ${d3d16p} ignored
                NAME3 client1 EXE3 ./proxy_coupling PROCS3 8 ARGS3 1 ${d3d8p} ignored)
    endif()
    if(HOST_NPROC GREATER 20)
        tri_mpi_test(TESTNAME test_proxy_coupling_40p
                TIMEOUT 16
                NAME1 rdv EXE1 ./proxy_coupling PROCS1 16 ARGS1 -1 ${d3d1p} ${d3d16p_cpn}
                NAME2 client0 EXE2 ./proxy_coupling PROCS2 16 ARGS2 0 ${d3d16p} ignored
                NAME3 client1 EXE3 ./proxy_coupling PROCS3 8 ARGS3 1 ${d3d8p} ignored)
    endif()
endif ()
# unit tests
find_package(Catch2 3)
if(Catch2_FOUND)
  message(STATUS "Found Catch2: ${Catch2_DIR} (found version ${Catch2_VERSION})")
  set(PCMS_UNIT_TEST_SOURCES
          unit_test_main.cpp
          test_coordinate_transform.cpp
          test_coordinate.cpp
          test_bounding_box.cpp)
  if (PCMS_ENABLE_XGC)
      list(APPEND PCMS_UNIT_TEST_SOURCES
              test_xgc_reverse_classifcation.cpp
              test_xgc_field_adapter.cpp)
  endif ()

  if (PCMS_ENABLE_OMEGA_H)
      add_executable(field_transfer_example field_transfer_example.cpp)
      target_link_libraries(field_transfer_example PUBLIC pcms::core)
      list(APPEND PCMS_UNIT_TEST_SOURCES
              test_field_transfer.cpp
              test_uniform_grid.cpp
              test_omega_h_copy.cpp
              test_point_search.cpp
              )
  endif ()
  add_executable(unit_tests ${PCMS_UNIT_TEST_SOURCES})
  target_link_libraries(unit_tests PUBLIC Catch2::Catch2 pcms::core)
  target_include_directories(unit_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  include(Catch)
  catch_discover_tests(unit_tests)
else()
message(WARNING "Catch2 not found. Disabling Unit Tests")
endif()

if (PCMS_ENABLE_C)
    find_package(Kokkos REQUIRED)
    if (PCMS_ENABLE_XGC)
        add_executable(test_proxy_couple_xgc_c_interface test_proxy_coupling_xgc_client.c)
        target_link_libraries(test_proxy_couple_xgc_c_interface PUBLIC pcms::capi MPI::MPI_C Kokkos::kokkos)
        if (PCMS_ENABLE_OMEGA_H)
            add_executable(test_proxy_couple_xgc_cpp_interface test_proxy_coupling_xgc_server.cpp)
            target_link_libraries(test_proxy_couple_xgc_cpp_interface PUBLIC pcms::core MPI::MPI_C test_support)
            dual_mpi_test(TESTNAME xgc_proxy_to_xgc
                    TIMEOUT 10
                    NAME1 cinterface EXE1 $<TARGET_FILE:test_proxy_couple_xgc_c_interface> PROCS1 2 ARGS1 ${PCMS_TEST_DATA_DIR}/d3d/meshRclassification.txt
                    NAME2 cppinterface EXE2 $<TARGET_FILE:test_proxy_couple_xgc_cpp_interface> PROCS2 1 ARGS2 ${d3d1p} ${d3d1p_cpn} 0)
            dual_mpi_test(TESTNAME xgc_proxy_to_omega
                    TIMEOUT 10
                    NAME1 cinterface EXE1 $<TARGET_FILE:test_proxy_couple_xgc_c_interface> PROCS1 4 ARGS1 ${PCMS_TEST_DATA_DIR}/d3d/meshRclassification.txt
                    NAME2 cppinterface EXE2 $<TARGET_FILE:test_proxy_couple_xgc_cpp_interface> PROCS2 2 ARGS2 ${d3d1p} ${d3d2p_cpn} 1)
        endif ()
    endif ()
endif ()
if(PCMS_ENABLE_Fortran)
    add_executable(test_proxy_coupling_xgc_client_fortran test_proxy_coupling_xgc_client_fortran.f90)
    target_link_libraries(test_proxy_coupling_xgc_client_fortran PUBLIC pcms::fortranapi MPI::MPI_Fortran)
    set_target_properties(test_proxy_coupling_xgc_client_fortran PROPERTIES LINKER_LANGUAGE CXX)
    dual_mpi_test(TESTNAME xgc_fortran_proxy_to_omega
            TIMEOUT 10
            NAME1 finterface EXE1 $<TARGET_FILE:test_proxy_coupling_xgc_client_fortran> PROCS1 2 ARGS1 ${PCMS_TEST_DATA_DIR}/d3d/meshRclassification.txt
            NAME2 cppinterface EXE2 $<TARGET_FILE:test_proxy_couple_xgc_cpp_interface> PROCS2 2 ARGS2 ${d3d1p} ${d3d2p_cpn} 1)
endif()

add_executable(xgc_n0_server xgc_n0_coupling_server.cpp)
target_link_libraries(xgc_n0_server PUBLIC pcms::core test_support)
