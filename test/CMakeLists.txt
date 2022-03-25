# TODO use submodule FetchContent/ExternalProject/ExternalData
set(WDMCPL_TEST_DATA_DIR "" CACHE PATH
        "Path to a local copy of the wdmapp_coupling_data repo.")
if (NOT EXISTS ${WDMCPL_TEST_DATA_DIR})
    message(FATAL_ERROR "WDMCPL_TEST_DATA_DIR \"${WDMCPL_TEST_DATA_DIR}\" is not accessible")
endif ()

message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")

function(add_exe NAME)
    add_executable(${NAME} ${NAME}.cpp)
    target_link_libraries(${NAME} wdmcpl)
    if (WDMCPL_HAS_ASAN)
        target_compile_options(${NAME} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
        target_link_libraries(${NAME} asan rt)
    endif ()
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

add_exe(test_meshRcbRdv)
mpi_test(test_meshRcbRdv_1p 1 ./test_meshRcbRdv
        ${WDMCPL_TEST_DATA_DIR}/Cyclone_ITG/Cyclone_ITG_deltaf_23mesh/mesh.osh/)


# unit tests
find_package(Catch2 REQUIRED)
add_executable(unit_tests unit_test_main.cpp test_coordinate_transform.cpp test_coordinate.cpp)
target_link_libraries(unit_tests PUBLIC Catch2::Catch2 wdmcpl::wdmcpl)
include(Catch)
catch_discover_tests(unit_tests)