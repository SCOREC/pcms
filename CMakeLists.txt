cmake_minimum_required(VERSION 3.19)
project(
  pcms
  VERSION 0.3.0
  LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# use _ROOT variables for configurations
cmake_policy(SET CMP0074 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

option(PCMS_ENABLE_ASAN "enable address sanitizer" OFF)
set(PCMS_HAS_ASAN OFF)
if(PCMS_ENABLE_ASAN AND CMAKE_COMPILER_IS_GNUCXX MATCHES 1)
  set(PCMS_HAS_ASAN ON)
endif()
option(PCMS_ENABLE_CLANG_TIDY "enable clang tidy target" OFF)

# based on Professional CMake: A Practical Guide
# this approach works around some issues with CMakes internal handling using co-compilation
# this approach affords a separate clang-tidy run step in the CI
if(PROJECT_IS_TOP_LEVEL)
  if(PCMS_ENABLE_CLANG_TIDY)
    set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
  
    find_program(CLANG_TIDY_EXECUTABLE clang-tidy REQUIRED)
    find_program(RUN_CLANG_TIDY_EXECUTABLE run-clang-tidy REQUIRED)
  
    set(CLANG_TIDY_CONFIG_FILE ${CMAKE_CURRENT_LIST_DIR}/.clang-tidy)
    set(CLANG_TIDY_SOURCE_FILTER [=[^(?!.*(\.c|\.i|\.f90)$).*$]=])
  
    message(STATUS "CLANG_TIDY_EXECUTABLE=${CLANG_TIDY_EXECUTABLE}")
    message(STATUS "RUN_CLANG_TIDY_EXECUTABLE= ${RUN_CLANG_TIDY_EXECUTABLE}")
    message(STATUS "CLANG_TIDY_CONFIG_FILE=${CLANG_TIDY_CONFIG_FILE}")
    message(STATUS "CLANG_TIDY_SOURCE_FILTER=${CLANG_TIDY_SOURCE_FILTER}")
  endif()
endif()

option(PCMS_ENABLE_SERVER "enable the coupling server implementation" ON)
option(PCMS_ENABLE_CLIENT "enable the coupling client implementation" ON)

option(PCMS_ENABLE_XGC "enable xgc field adapter" ON)
option(PCMS_ENABLE_OMEGA_H "enable Omega_h field adapter" OFF)
option(PCMS_ENABLE_C "Enable pcms C api" ON)
option(PCMS_ENABLE_Python "Enable pcms Python api" OFF)

option(PCMS_ENABLE_PRINT "PCMS print statements enabled" ON)

option(PETSC_LINK_STATIC "Use pkg-config --static results for PETSc"
       ${_pcms_link_petsc_static_default})

# find package before fortran enabled, so we don't require the adios2 fortran
# interfaces this is important because adios2 build with clang/gfortran is
# broken
find_package(redev 4.3.0 REQUIRED)

if (PCMS_ENABLE_Python)
  set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()
if(PCMS_ENABLE_C)
  enable_language(C)
  option(PCMS_ENABLE_Fortran "Enable pcms fortran api" ON)
  if(PCMS_ENABLE_Fortran)
    enable_language(Fortran)
  endif()
endif()
if(PCMS_ENABLE_SERVER)
  set(PCMS_ENABLE_OMEGA_H
      ON
      CACHE BOOL "enable Omega_h field adapter" FORCE)
endif()

set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)

# Create the run-clang-tidy target here (after find_package(MPI)) so we can
# forward any MPI include paths that are missing from compile_commands.json.
# When CMAKE_CXX_COMPILER is an MPI wrapper (e.g. mpicxx), FindMPI may leave
# MPI_CXX_INCLUDE_DIRS empty because the wrapper handles includes implicitly,
# meaning compile_commands.json has no explicit -I flags for mpi.h.
if(PROJECT_IS_TOP_LEVEL)
  if(PCMS_ENABLE_CLANG_TIDY)
    set(_clang_tidy_extra_args "")
    get_target_property(_mpi_cxx_incs MPI::MPI_CXX INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT _mpi_cxx_incs)
      execute_process(
        COMMAND ${MPI_CXX_COMPILER} -show
        OUTPUT_VARIABLE _mpi_show
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      string(REGEX MATCHALL "-I([^ ]+)" _mpi_inc_flags "${_mpi_show}")
      foreach(_flag IN LISTS _mpi_inc_flags)
        list(APPEND _clang_tidy_extra_args "-extra-arg=${_flag}")
      endforeach()
      message(STATUS "clang-tidy MPI extra args: ${_clang_tidy_extra_args}")
    endif()
  
    add_custom_target(run-clang-tidy
      VERBATIM COMMAND ${RUN_CLANG_TIDY_EXECUTABLE}
      -clang-tidy-binary ${CLANG_TIDY_EXECUTABLE}
      -p ${CMAKE_BINARY_DIR}
      -quiet
      -config-file ${CLANG_TIDY_CONFIG_FILE}
      -source-filter ${CLANG_TIDY_SOURCE_FILTER}
      ${_clang_tidy_extra_args}
    )
  endif()
endif()

message(STATUS "Found redev: ${redev_DIR} (found version ${redev_VERSION})")
if(PCMS_ENABLE_OMEGA_H)
  find_package(Omega_h REQUIRED VERSION 10)
  message(
    STATUS "Found Omega_h: ${Omega_h_DIR} (found version ${Omega_h_VERSION})")
  if(NOT Omega_h_USE_MPI)
    message(FATAL_ERROR "Omega_h must be built with MPI enabled.")
  endif()
endif()
# adios2 adds C and Fortran depending on how it was built
find_package(ADIOS2 CONFIG 2.10.2 REQUIRED)
find_package(Kokkos CONFIG 4.5 REQUIRED)

option(PCMS_ENABLE_MESHFIELDS "Enable MeshFields support" ON)
option(PCMS_ENABLE_PETSC "Enable PETSc support" OFF)

if(PCMS_ENABLE_MESHFIELDS)
  find_package(meshfields REQUIRED)
  message(STATUS "Found MeshFields: ${meshfields_DIR} (found version ${meshfields_VERSION})")
endif()

if(PCMS_ENABLE_PETSC)
  find_package(PETSc REQUIRED)
endif()

add_subdirectory(src)

include(CTest)
include(test/testing.cmake)
# define 'smoke tests' to test the install
add_custom_target(test_install DEPENDS check) # maintain test_install target
if(NOT BUILD_TESTING)
  add_custom_target(
    check
    COMMAND ${CMAKE_CTEST_COMMAND} -R smoke_test
    COMMENT "Test installed PUMIPic utilities")
endif()

add_executable(smoke_test test/test_ohClassPtn.cpp)
target_link_libraries(smoke_test pcms::core test_support)
dual_mpi_test(
  TESTNAME
  smoke_test
  TIMEOUT
  4
  NAME1
  rdv
  EXE1
  ./smoke_test
  PROCS1
  2
  ARGS1
  1
  ${PCMS_TEST_DATA_DIR}/cyclone/23elements/mesh.osh/
  NAME2
  app
  EXE2
  ./smoke_test
  PROCS2
  1
  ARGS2
  0
  ${PCMS_TEST_DATA_DIR}/cyclone/23elements/mesh.osh/)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
add_subdirectory(tools)
