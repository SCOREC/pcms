cmake_minimum_required(VERSION 3.19)
project(wdm_coupling VERSION 0.0.1 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(WDMCPL_ENABLE_ASAN "enable address sanitizer" OFF)
set(WDMCPL_HAS_ASAN OFF)
if(WDMCPL_ENABLE_ASAN AND CMAKE_COMPILER_IS_GNUCXX MATCHES 1)
  set(WDMCPL_HAS_ASAN ON)
endif()

set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)
find_package(redev REQUIRED)
message(STATUS "Found redev: ${redev_DIR} (found version ${redev_VERSION})")
find_package(Omega_h REQUIRED VERSION 10)
message(STATUS "Found Omega_h: ${Omega_h_DIR} (found version ${Omega_h_VERSION})")

## use pkgconfig since the fftw autoconf install produces
## broken cmake config files
## https://github.com/FFTW/fftw3/issues/130
find_package(PkgConfig REQUIRED)
pkg_check_modules(fftw REQUIRED IMPORTED_TARGET fftw3>=3.3)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH) 
  set(WDMCPL_TEST_DATA_DIR "" CACHE PATH
    "Path to a local copy of the wdmapp_coupling_data repo.")
  if(NOT EXISTS ${WDMCPL_TEST_DATA_DIR})
    message(FATAL_ERROR "WDMCPL_TEST_DATA_DIR \"${WDMCPL_TEST_DATA_DIR}\" is not accessible")
  endif()
endif()

