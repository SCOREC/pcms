set(SOURCES coupling.cc)
set(HEADERS coupling.h)

add_library(coupler ${SOURCES})
target_link_libraries(coupler adios2::adios2 Kokkos::kokkos PkgConfig::fftw)
add_executable(cpl cpl.cc)
target_link_libraries(cpl coupler)
