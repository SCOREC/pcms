set(SOURCES coupling.cc)
set(HEADERS coupling.h)

add_library(coupler ${SOURCES})
set_target_properties(coupler PROPERTIES PUBLIC_HEADER ${HEADERS})
target_link_libraries(coupler adios2::adios2 Kokkos::kokkos PkgConfig::fftw)

add_executable(cpl cpl.cc)
target_link_libraries(cpl coupler)

install(TARGETS cpl coupler)
