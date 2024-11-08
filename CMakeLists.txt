cmake_minimum_required(VERSION 3.19)
project(pcms VERSION 0.0.5 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# use _ROOT variables for configurations
cmake_policy(SET CMP0074 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

option(PCMS_ENABLE_ASAN "enable address sanitizer" OFF)
set(PCMS_HAS_ASAN OFF)
if (PCMS_ENABLE_ASAN AND CMAKE_COMPILER_IS_GNUCXX MATCHES 1)
  set(PCMS_HAS_ASAN ON)
endif ()

option(PCMS_ENABLE_SERVER "enable the coupling server implementation" ON)
option(PCMS_ENABLE_CLIENT "enable the coupling client implementation" ON)

option(PCMS_ENABLE_XGC "enable xgc field adapter" ON)
option(PCMS_ENABLE_OMEGA_H "enable Omega_h field adapter" OFF)
option(PCMS_ENABLE_C "Enable pcms C api" ON)

# find package before fortran enabled, so we don't require the adios2 fortran interfaces
# this is important because adios2 build with clang/gfortran is broken
find_package(redev 4.3.0 REQUIRED)

if (PCMS_ENABLE_C)
  enable_language(C)
  option(PCMS_ENABLE_Fortran "Enable pcms fortran api" ON)
  if(PCMS_ENABLE_Fortran)
    enable_language(Fortran)
  endif()
endif ()
if (PCMS_ENABLE_SERVER)
  set(PCMS_ENABLE_OMEGA_H ON CACHE BOOL "enable Omega_h field adapter" FORCE)
endif ()


set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)
message(STATUS "Found redev: ${redev_DIR} (found version ${redev_VERSION})")
if (PCMS_ENABLE_OMEGA_H)
  find_package(Omega_h REQUIRED VERSION 10)
  message(STATUS "Found Omega_h: ${Omega_h_DIR} (found version ${Omega_h_VERSION})")
  if(NOT Omega_h_USE_MPI)
    message(FATAL_ERROR "Omega_h must be built with MPI enabled.")
  endif()
endif()

## use pkgconfig since the fftw autoconf install produces
## broken cmake config files
## https://github.com/FFTW/fftw3/issues/130
find_package(PkgConfig REQUIRED)
pkg_check_modules(fftw REQUIRED IMPORTED_TARGET fftw3>=3.3)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
add_subdirectory(tools)
