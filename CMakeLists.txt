cmake_minimum_required(VERSION 3.13)
project(wdm_coupling VERSION 0.0.1 LANGUAGES CXX)

set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)

#adios2 adds C and Fortran depending on how it was built
find_package(ADIOS2 2.5 REQUIRED)
find_package(Kokkos 3.0 REQUIRED)

## use pkgconfig since the fftw autoconf install produces
## broken cmake config files
## https://github.com/FFTW/fftw3/issues/130
find_package(PkgConfig REQUIRED)
pkg_check_modules(fftw REQUIRED IMPORTED_TARGET fftw3>=3.3)

add_subdirectory(src)