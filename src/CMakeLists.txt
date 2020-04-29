set(SOURCES 
  coupling.cc
  commpart1.cc
  importpart3mesh.cc
  BoundExchange.cc
  commpart1.cc
  fourierdataproc.cc
  initdatastruc.cc 
  interpo.cc
  testutilities.cc
)

set(HEADERS
    coupling.h 
    couplingConstants.h
    commpart1.h
    dataprocess.h
    importpart3mesh.h
    commwrapping.h
    parameters.h
    testutilities.h
)

add_library(coupler ${SOURCES})
set_target_properties(coupler PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_link_libraries(coupler adios2::adios2 Kokkos::kokkos PkgConfig::fftw)
target_include_directories(coupler PUBLIC ./)
