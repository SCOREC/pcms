include(ProcessorCount)
ProcessorCount(HOST_NPROC)
message("TESTING ON SYSTEM WITH ${HOST_NPROC} CPU Cores")

# helper variables for readability
set(rendezvous 1)
set(notRendezvous 0)

if(PCMS_ENABLE_OMEGA_H)
  add_exe(test_ohClassPtn)
  set(cyclone1p ${PCMS_TEST_DATA_DIR}/cyclone/23elements/mesh.osh/)
  if(HOST_NPROC GREATER_EQUAL 3)
    dual_mpi_test(
      TESTNAME
      test_ohClassPtnSendRecv
      TIMEOUT
      4
      NAME1
      rdv
      EXE1
      ./test_ohClassPtn
      PROCS1
      2
      ARGS1
      ${rendezvous}
      ${cyclone1p}
      NAME2
      app
      EXE2
      ./test_ohClassPtn
      PROCS2
      1
      ARGS2
      ${notRendezvous}
      ${cyclone1p})
  endif()
  add_exe(test_ohClassPtn_appRibPtn)
  set(cyclone2p ${PCMS_TEST_DATA_DIR}/cyclone/23elements/2p.osh/)
  if(HOST_NPROC GREATER_EQUAL 4)
    dual_mpi_test(
      TESTNAME
      test_ohClassPtn_appRibPtn
      TIMEOUT
      8
      NAME1
      rdv
      EXE1
      ./test_ohClassPtn_appRibPtn
      PROCS1
      2
      ARGS1
      ${rendezvous}
      ${cyclone1p}
      NAME2
      app
      EXE2
      ./test_ohClassPtn_appRibPtn
      PROCS2
      2
      ARGS2
      ${notRendezvous}
      ${cyclone2p})
  endif()
  set(d3d1p ${PCMS_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc.osh/)
  set(d3d1p_cpn ${PCMS_TEST_DATA_DIR}/d3d/1p.cpn)
  set(d3d2p_cpn ${PCMS_TEST_DATA_DIR}/d3d/2p.cpn)
  set(d3d4p_cpn ${PCMS_TEST_DATA_DIR}/d3d/4p.cpn)
  set(d3d8p_cpn ${PCMS_TEST_DATA_DIR}/d3d/8p.cpn)
  set(d3d16p_cpn ${PCMS_TEST_DATA_DIR}/d3d/16p.cpn)
  set(d3d16p ${PCMS_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc_p16.osh/)
  add_exe(test_ohOverlap)
  if(HOST_NPROC GREATER_EQUAL 20)
    dual_mpi_test(
      TESTNAME
      test_ohOverlap_d3d_20p
      TIMEOUT
      20
      NAME1
      rdv
      EXE1
      ./test_ohOverlap
      PROCS1
      4
      ARGS1
      ${rendezvous}
      ${d3d1p}
      ${d3d4p_cpn}
      NAME2
      app
      EXE2
      ./test_ohOverlap
      PROCS2
      16
      ARGS2
      ${notRendezvous}
      ${d3d16p}
      ignored)
  endif()

  set(d3d8p ${PCMS_TEST_DATA_DIR}/d3d/d3d-full_9k_sfc_p8.osh/)
  add_exe(test_twoClientOverlap)
  if(HOST_NPROC GREATER_EQUAL 28)
    tri_mpi_test(
      TESTNAME
      test_twoClientOverlap_d3d_28p
      TIMEOUT
      20
      NAME1
      rdv
      EXE1
      ./test_twoClientOverlap
      PROCS1
      4
      ARGS1
      -1
      ${d3d1p}
      ${d3d4p_cpn}
      NAME2
      client0
      EXE2
      ./test_twoClientOverlap
      PROCS2
      16
      ARGS2
      0
      ${d3d16p}
      ignored
      NAME3
      client1
      EXE3
      ./test_twoClientOverlap
      PROCS3
      8
      ARGS3
      1
      ${d3d8p}
      ignored)
  endif()
  add_executable(proxy_coupling test_proxy_coupling.cpp)
  target_link_libraries(proxy_coupling PUBLIC pcms::core test_support)
  if(HOST_NPROC GREATER_EQUAL 4)
    tri_mpi_test(
      TESTNAME
      test_proxy_coupling_4p
      TIMEOUT
      10
      NAME1
      rdv
      EXE1
      ./proxy_coupling
      PROCS1
      2
      ARGS1
      -1
      ${d3d1p}
      ${d3d2p_cpn}
      NAME2
      client0
      EXE2
      ./proxy_coupling
      PROCS2
      1
      ARGS2
      0
      ${d3d1p}
      ignored
      NAME3
      client1
      EXE3
      ./proxy_coupling
      PROCS3
      1
      ARGS3
      1
      ${d3d1p}
      ignored)
  endif()
  if(HOST_NPROC GREATER_EQUAL 26)
    tri_mpi_test(
      TESTNAME
      test_proxy_coupling_26p
      TIMEOUT
      20
      NAME1
      rdv
      EXE1
      ./proxy_coupling
      PROCS1
      2
      ARGS1
      -1
      ${d3d1p}
      ${d3d2p_cpn}
      NAME2
      client0
      EXE2
      ./proxy_coupling
      PROCS2
      16
      ARGS2
      0
      ${d3d16p}
      ignored
      NAME3
      client1
      EXE3
      ./proxy_coupling
      PROCS3
      8
      ARGS3
      1
      ${d3d8p}
      ignored)
  endif()
  if(HOST_NPROC GREATER_EQUAL 28)
    tri_mpi_test(
      TESTNAME
      test_proxy_coupling_28p
      TIMEOUT
      20
      NAME1
      rdv
      EXE1
      ./proxy_coupling
      PROCS1
      4
      ARGS1
      -1
      ${d3d1p}
      ${d3d4p_cpn}
      NAME2
      client0
      EXE2
      ./proxy_coupling
      PROCS2
      16
      ARGS2
      0
      ${d3d16p}
      ignored
      NAME3
      client1
      EXE3
      ./proxy_coupling
      PROCS3
      8
      ARGS3
      1
      ${d3d8p}
      ignored)
  endif()
  if(HOST_NPROC GREATER_EQUAL 32)
    tri_mpi_test(
      TESTNAME
      test_proxy_coupling_32p
      TIMEOUT
      20
      NAME1
      rdv
      EXE1
      ./proxy_coupling
      PROCS1
      8
      ARGS1
      -1
      ${d3d1p}
      ${d3d8p_cpn}
      NAME2
      client0
      EXE2
      ./proxy_coupling
      PROCS2
      16
      ARGS2
      0
      ${d3d16p}
      ignored
      NAME3
      client1
      EXE3
      ./proxy_coupling
      PROCS3
      8
      ARGS3
      1
      ${d3d8p}
      ignored)
  endif()
  if(HOST_NPROC GREATER_EQUAL 40)
    tri_mpi_test(
      TESTNAME
      test_proxy_coupling_40p
      TIMEOUT
      16
      NAME1
      rdv
      EXE1
      ./proxy_coupling
      PROCS1
      16
      ARGS1
      -1
      ${d3d1p}
      ${d3d16p_cpn}
      NAME2
      client0
      EXE2
      ./proxy_coupling
      PROCS2
      16
      ARGS2
      0
      ${d3d16p}
      ignored
      NAME3
      client1
      EXE3
      ./proxy_coupling
      PROCS3
      8
      ARGS3
      1
      ${d3d8p}
      ignored)
  endif()
endif()
# unit tests
find_package(Catch2 3)
if(Catch2_FOUND)
  message(
    STATUS "Found Catch2: ${Catch2_DIR} (found version ${Catch2_VERSION})")
  set(PCMS_UNIT_TEST_SOURCES unit_test_main.cpp test_coordinate_transform.cpp
                             test_coordinate.cpp test_bounding_box.cpp)
  if(PCMS_ENABLE_XGC)
    list(APPEND PCMS_UNIT_TEST_SOURCES test_xgc_reverse_classification.cpp
         test_xgc_field_adapter.cpp)
  endif()

  if(PCMS_ENABLE_OMEGA_H)
    add_executable(field_transfer_example field_transfer_example.cpp)
    target_link_libraries(field_transfer_example PUBLIC pcms::core
                                                        pcms_interpolator)
    list(
      APPEND
      PCMS_UNIT_TEST_SOURCES
      test_field_transfer.cpp
      test_uniform_grid.cpp
      test_omega_h_copy.cpp
      test_point_search.cpp
      test_mls_basis.cpp
      test_rbf_interp.cpp
      test_normalisation.cpp
      test_spline_interpolator.cpp
      test_svd_serial.cpp)
  endif()
  add_executable(unit_tests ${PCMS_UNIT_TEST_SOURCES})
  target_link_libraries(unit_tests PUBLIC Catch2::Catch2 pcms::core
                                          pcms_interpolator)
  target_include_directories(unit_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  include(Catch)
  catch_discover_tests(unit_tests)
else()
  message(WARNING "Catch2 not found. Disabling Unit Tests")
endif()

if(PCMS_ENABLE_C)
  find_package(Kokkos REQUIRED)
  if(PCMS_ENABLE_XGC)
    add_executable(test_proxy_couple_xgc_c_interface
                   test_proxy_coupling_xgc_client.c)
    target_link_libraries(test_proxy_couple_xgc_c_interface
                          PUBLIC pcms::capi MPI::MPI_C Kokkos::kokkos)
    if(PCMS_ENABLE_OMEGA_H)
      add_executable(test_proxy_couple_xgc_cpp_interface
                     test_proxy_coupling_xgc_server.cpp)
      target_link_libraries(test_proxy_couple_xgc_cpp_interface
                            PUBLIC pcms::core MPI::MPI_C test_support)
      if(HOST_NPROC GREATER_EQUAL 3)
        dual_mpi_test(
          TESTNAME
          xgc_proxy_to_xgc
          TIMEOUT
          10
          NAME1
          cinterface
          EXE1
          $<TARGET_FILE:test_proxy_couple_xgc_c_interface>
          PROCS1
          2
          ARGS1
          ${PCMS_TEST_DATA_DIR}/d3d/meshRclassification.txt
          NAME2
          cppinterface
          EXE2
          $<TARGET_FILE:test_proxy_couple_xgc_cpp_interface>
          PROCS2
          1
          ARGS2
          ${d3d1p}
          ${d3d1p_cpn}
          0)
      endif()
      if(HOST_NPROC GREATER_EQUAL 6)
        dual_mpi_test(
          TESTNAME
          xgc_proxy_to_omega
          TIMEOUT
          10
          NAME1
          cinterface
          EXE1
          $<TARGET_FILE:test_proxy_couple_xgc_c_interface>
          PROCS1
          4
          ARGS1
          ${PCMS_TEST_DATA_DIR}/d3d/meshRclassification.txt
          NAME2
          cppinterface
          EXE2
          $<TARGET_FILE:test_proxy_couple_xgc_cpp_interface>
          PROCS2
          2
          ARGS2
          ${d3d1p}
          ${d3d2p_cpn}
          1)
      endif()
    endif()
  endif()
endif()
if(PCMS_ENABLE_Fortran)
  add_executable(test_proxy_coupling_xgc_client_fortran
                 test_proxy_coupling_xgc_client_fortran.f90)
  target_link_libraries(test_proxy_coupling_xgc_client_fortran
                        PUBLIC pcms::fortranapi MPI::MPI_Fortran)
  set_target_properties(test_proxy_coupling_xgc_client_fortran
                        PROPERTIES LINKER_LANGUAGE Fortran)
  if(HOST_NPROC GREATER_EQUAL 4)
    dual_mpi_test(
      TESTNAME
      xgc_fortran_proxy_to_omega
      TIMEOUT
      10
      NAME1
      finterface
      EXE1
      $<TARGET_FILE:test_proxy_coupling_xgc_client_fortran>
      PROCS1
      2
      ARGS1
      ${PCMS_TEST_DATA_DIR}/d3d/meshRclassification.txt
      NAME2
      cppinterface
      EXE2
      $<TARGET_FILE:test_proxy_couple_xgc_cpp_interface>
      PROCS2
      2
      ARGS2
      ${d3d1p}
      ${d3d2p_cpn}
      1)
  endif()
endif()

add_executable(xgc_n0_server xgc_n0_coupling_server.cpp)
target_link_libraries(xgc_n0_server PUBLIC pcms::core test_support)
