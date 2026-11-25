cmake_minimum_required(VERSION 3.26)
project(
  pcms
  VERSION 0.4.0
  LANGUAGES CXX)

include(GNUInstallDirs)

# Report the accelerator inventory of the CI runner (diagnostic aid).
execute_process(
  COMMAND sh -c "nvidia-smi -L 2>&1; nvidia-smi --query-gpu=name,uuid,memory.total,memory.free,driver_version --format=csv,noheader 2>&1; hostname; id"
  OUTPUT_VARIABLE PCMS_ACCEL_INVENTORY
  ERROR_VARIABLE PCMS_ACCEL_INVENTORY_ERR)
message(STATUS "=== accelerator inventory ===")
message(STATUS "${PCMS_ACCEL_INVENTORY}")
message(STATUS "=== end accelerator inventory ===")
# The package config files go to <prefix>/lib/cmake/pcms on every platform.
set(PCMS_INSTALL_CMAKEDIR "lib/cmake/pcms"
    CACHE STRING "Install directory of the pcms CMake package files")
include(CMakePackageConfigHelpers)
# Installed libraries find each other and their private dependencies through
# their RPATH, so consumers of a shared PCMS do not need the private
# dependencies on their link line (the package config finds them only for a
# static PCMS).
if(NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()
if(NOT DEFINED CMAKE_INSTALL_RPATH AND NOT APPLE)
  set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()

# use _ROOT variables for configurations
cmake_policy(SET CMP0074 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
# map options renamed in 0.4 (e.g. PCMS_ENABLE_Python) to their new names
include(PcmsRenamedOptions)

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

# Every optional capability is off by default; a build with none of them needs
# only Kokkos and MPI (see the presets in CMakePresets.json). Each flag names
# the backend namespace or feature it builds (pcms::omega_h, pcms.omega_h in
# Python, PCMS_ENABLE_OMEGA_H in pcms/configuration.hpp, ...).
option(PCMS_ENABLE_OMEGA_H "Omega_h backend (needs Omega_h)" OFF)
option(PCMS_ENABLE_OMEGA_H_MESHFIELDS
       "MeshFields Lagrange spaces on Omega_h meshes (needs OMEGA_H, MeshFields)" OFF)
option(PCMS_ENABLE_CONSERVATIVE_PROJECTION
       "Conservative and control-variate projections (needs OMEGA_H_MESHFIELDS, PETSc)" OFF)
option(PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION
       "Polynomial-reconstruction spaces on point clouds (needs KokkosKernels)" OFF)
option(PCMS_ENABLE_XGC "XGC backend (needs OMEGA_H)" OFF)
option(PCMS_ENABLE_MFEM "MFEM backend (needs MFEM)" OFF)
option(PCMS_ENABLE_COUPLER "Coupler (needs redev and ADIOS2)" OFF)
option(PCMS_ENABLE_C
       "C API: the client part with COUPLER and XGC, the transfer part with OMEGA_H and POLYNOMIAL_RECONSTRUCTION"
       OFF)
option(PCMS_ENABLE_FORTRAN "Fortran API (needs C)" OFF)
option(PCMS_ENABLE_PYTHON "Python module (needs pybind11)" OFF)
option(PCMS_ENABLE_EXAMPLES "Build the in-tree examples" OFF)
option(PCMS_ENABLE_PRINT "PCMS print statements enabled" ON)

# Options with no effect any more
pcms_removed_option(PCMS_ENABLE_SERVER
                    "the coupler no longer has a separate server part")
pcms_removed_option(PCMS_ENABLE_CLIENT
                    "the coupler no longer has a separate client part")
pcms_removed_option(PCMS_ENABLE_TRANSIENT "the timestepper is always built")

# A capability whose prerequisite capability is off is an error, not a silent
# change of the configuration.
function(pcms_require capability prerequisite)
  if(PCMS_ENABLE_${capability} AND NOT PCMS_ENABLE_${prerequisite})
    message(FATAL_ERROR "PCMS_ENABLE_${capability} needs PCMS_ENABLE_${prerequisite}: "
                        "configure with -DPCMS_ENABLE_${prerequisite}=ON.")
  endif()
endfunction()
pcms_require(OMEGA_H_MESHFIELDS OMEGA_H)
pcms_require(CONSERVATIVE_PROJECTION OMEGA_H_MESHFIELDS)
pcms_require(XGC OMEGA_H)
pcms_require(FORTRAN C)
# PCMS_ENABLE_C builds each part of the C (and Fortran) API whose capabilities
# are on: the client wraps the coupler with XGC fields, the transfer part
# interpolates on Omega_h meshes by polynomial reconstruction.
set(PCMS_BUILD_CAPI_CLIENT OFF)
set(PCMS_BUILD_CAPI_TRANSFER OFF)
if(PCMS_ENABLE_C)
  if(PCMS_ENABLE_COUPLER AND PCMS_ENABLE_XGC)
    set(PCMS_BUILD_CAPI_CLIENT ON)
  endif()
  if(PCMS_ENABLE_OMEGA_H AND PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION)
    set(PCMS_BUILD_CAPI_TRANSFER ON)
  endif()
  if(NOT PCMS_BUILD_CAPI_CLIENT AND NOT PCMS_BUILD_CAPI_TRANSFER)
    message(FATAL_ERROR "PCMS_ENABLE_C builds the C client API with "
                        "PCMS_ENABLE_COUPLER and PCMS_ENABLE_XGC and the C "
                        "transfer API with PCMS_ENABLE_OMEGA_H and "
                        "PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION: turn on at "
                        "least one of the two pairs.")
  endif()
endif()

# Find the dependency of a capability (the package and the arguments of
# find_package, without REQUIRED); fail naming the capability if it is missing.
macro(pcms_find_package capability package)
  find_package(${package} ${ARGN})
  if(NOT ${package}_FOUND)
    message(FATAL_ERROR "PCMS_ENABLE_${capability} needs ${package}, which was "
                        "not found: add its prefix to CMAKE_PREFIX_PATH or "
                        "set ${package}_DIR.")
  endif()
endmacro()

# Per-thread on-stack ring-buffer capacities for the adjacency-based
# intersection/adjacency BFS (queue_visited.hpp). Defaults are sized for 3D
# tetrahedral meshes on the serial/OpenMP backends. Lower them for device
# (GPU) builds, where each search thread's stack must hold these arrays.
set(PCMS_INTERSECTION_QUEUE_SIZE 1024 CACHE STRING
    "Max BFS queue size per target element in the intersection search")
set(PCMS_INTERSECTION_TRACK_SIZE 2048 CACHE STRING
    "Max visited-set size per target element in the intersection search")
set(PCMS_INTERSECTION_ABS_TOL 1e-18 CACHE STRING
    "Absolute measure tolerance for accepting a clipped intersection")
set(PCMS_INTERSECTION_REL_TOL 1e-12 CACHE STRING
    "Relative measure tolerance for accepting a clipped intersection")
set(PCMS_PROJECTION_KSP_RTOL 1e-14 CACHE STRING
    "Default relative tolerance for the projection mass-matrix solve")
set(PCMS_PROJECTION_KSP_ATOL 1e-50 CACHE STRING
    "Default absolute tolerance for the projection mass-matrix solve")

option(PCMS_PETSC_LINK_STATIC "Use pkg-config --static results for PETSc"
       ${_pcms_link_petsc_static_default})

# find package before fortran enabled, so we don't require the adios2 fortran
# interfaces this is important because adios2 build with clang/gfortran is
# broken
if(PCMS_ENABLE_COUPLER)
  pcms_find_package(COUPLER redev 4.3.0)
endif()

if (PCMS_ENABLE_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()

if(PCMS_ENABLE_C)
  enable_language(C)
endif()
if(PCMS_ENABLE_FORTRAN)
  enable_language(Fortran)
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

if(PCMS_ENABLE_COUPLER)
  message(STATUS "Found redev: ${redev_DIR} (found version ${redev_VERSION})")
  # adios2 adds C and Fortran depending on how it was built
  pcms_find_package(COUPLER ADIOS2 CONFIG 2.10.2)
endif()
find_package(Kokkos CONFIG 4.5 REQUIRED)
if(PCMS_ENABLE_OMEGA_H)
  pcms_find_package(OMEGA_H Omega_h)
  message(
    STATUS "Found Omega_h: ${Omega_h_DIR} (found version ${Omega_h_VERSION})")
  if(NOT Omega_h_USE_MPI)
    message(FATAL_ERROR "Omega_h must be built with MPI enabled.")
  endif()
  if(PCMS_ENABLE_PYTHON AND Omega_h_VERSION VERSION_LESS 11.4)
    message(FATAL_ERROR "PCMS Python with PCMS_ENABLE_OMEGA_H requires "
                        "Omega_h 11.4 or newer.")
  endif()
endif()
if(PCMS_ENABLE_OMEGA_H_MESHFIELDS)
  pcms_find_package(OMEGA_H_MESHFIELDS meshfields)
  message(STATUS "Found MeshFields: ${meshfields_DIR} (found version ${meshfields_VERSION})")
endif()
if(PCMS_ENABLE_CONSERVATIVE_PROJECTION)
  pcms_find_package(CONSERVATIVE_PROJECTION PETSc)
endif()
if(PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION)
  # the MLS solve uses KokkosBatched
  pcms_find_package(POLYNOMIAL_RECONSTRUCTION KokkosKernels)
endif()
if(PCMS_ENABLE_MFEM)
  pcms_find_package(MFEM MFEM)
  message(STATUS "Found MFEM: ${MFEM_DIR} (found version ${MFEM_VERSION})")
  if(NOT MFEM_USE_MPI)
    message(FATAL_ERROR "MFEM must be built with MPI enabled.")
  endif()
endif()

include(CTest)
# Every installed header is in a file set of its library; with testing on,
# the target all_verify_interface_header_sets compiles each one on its own
# (the test public_headers_self_contained).
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTING)
  set(CMAKE_VERIFY_INTERFACE_HEADER_SETS ON)
endif()
add_subdirectory(src)

include(test/testing.cmake)
# define 'smoke tests' to test the install
add_custom_target(test_install DEPENDS check) # maintain test_install target
if(NOT BUILD_TESTING)
  # include(CTest) only enables testing when BUILD_TESTING is on. Enable it so
  # smoke_test is still registered and 'check' can run it.
  enable_testing()
  add_custom_target(
    check
    COMMAND ${CMAKE_CTEST_COMMAND} -R smoke_test
    COMMENT "Test installed PUMIPic utilities")
endif()

# the smoke test exercises the coupler and needs the Omega_h test support
if(TARGET test_support)
  add_executable(smoke_test test/test_ohClassPtn.cpp)
  target_link_libraries(smoke_test pcms::coupler test_support)
  # Uses a mesh vendored in test/data so the install check works without
  # PCMS_TEST_DATA_DIR (e.g. with BUILD_TESTING=OFF).
  dual_mpi_test(
    TESTNAME
    smoke_test
    TIMEOUT
    4
    NAME1
    rdv
    EXE1
    smoke_test
    PROCS1
    2
    ARGS1
    1
    ${CMAKE_CURRENT_SOURCE_DIR}/test/data/cyclone23/mesh.osh/
    NAME2
    app
    EXE2
    smoke_test
    PROCS2
    1
    ARGS2
    0
    ${CMAKE_CURRENT_SOURCE_DIR}/test/data/cyclone23/mesh.osh/)
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
add_subdirectory(tools)
if(PCMS_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
