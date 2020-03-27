cmake_minimum_required(VERSION 3.12)
project(wdm_coupling VERSION 0.0.1 LANGUAGES CXX)
#adios2 adds C and Fortran depending on how it was built
find_package(ADIOS2 REQUIRED)
find_package(Kokkos REQUIRED)
add_subdirectory(src)
