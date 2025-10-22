cmake_minimum_required(VERSION 3.19)
project(
  pcms
  VERSION 0.2.0
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

option(PCMS_ENABLE_SERVER "enable the coupling server implementation" ON)
option(PCMS_ENABLE_CLIENT "enable the coupling client implementation" ON)

option(PCMS_ENABLE_XGC "enable xgc field adapter" ON)
option(PCMS_ENABLE_OMEGA_H "enable Omega_h field adapter" OFF)
option(PCMS_ENABLE_C "Enable pcms C api" ON)

option(PCMS_PRINT_ENABLED "PCMS print statements enabled" ON)
find_package(spdlog QUIET)

# find package before fortran enabled, so we don't require the adios2 fortran
# interfaces this is important because adios2 build with clang/gfortran is
# broken
find_package(redev 4.3.0 REQUIRED)

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
find_package(ADIOS2 CONFIG 2.5 REQUIRED)
find_package(Kokkos CONFIG 4.5 REQUIRED)

option(PCMS_ENABLE_MESHFIELDS "Enable MeshFields support" OFF)
option(PCMS_ENABLE_PETSC "Enable PETSc support" OFF)

if (PCMS_ENABLE_MESHFIELDS)
	find_package(meshfields REQUIRED)
	message(STATUS "Found MeshFields: ${meshfields_DIR} (found version ${meshfields_VERSION})")
endif()

if (PCMS_ENABLE_PETSC)
	if (NOT PETSC_DIR OR NOT PETSC_ARCH)
		message(FATAL_ERROR "PETSC_DIR and PETSC_ARCH must be provided when PETSc is not CMake-enabled")
	endif()
	set(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include" "${PETSC_DIR}/${PETSC_ARCH}/include")
	set(PETSC_LIB_DIR "${PETSC_DIR}/${PETSC_ARCH}/lib")
	find_library(PETSC_LIBRARY NAMES petsc HINTS ${PETSC_LIB_DIR} REQUIRED)

	if (NOT PETSC_LIBRARY)
		message(FATAL_ERROR "Couldnot find libpetsc in ${PETSC_LIB_DIR}")
	endif()

	set(PCMS_PETSC_INCLUDES ${PETSC_INCLUDE_DIRS} CACHE INTERNAL "")
	set(PCMS_PETSC_LIBRARIES ${PETSC_LIBRARY} CACHE INTERNAL "")

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
