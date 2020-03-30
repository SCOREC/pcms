set(SOURCES coupling.cc)
set(HEADERS)

add_executable(cpl ${SOURCES})
target_link_libraries(cpl adios2::adios2 Kokkos::kokkos PkgConfig::fftw)
