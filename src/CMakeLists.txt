set(SOURCES 
  coupling.cc
  commpart1.cc
  importpart3mesh.cc
  BoundaryDescr3D.cc
  commpart1.cc
  fourierdataproc.cc
  dataprocess.cc 
  interpo.cc
  inittestenv.cc
)

set(HEADERS
    coupling.h 
    couplingConstants.h
    couplingTypes.h
    commpart1.h
    dataprocess.h
    importpart3mesh.h
    BoundaryDescr3D.h
    commwrapping.h
    testutilities.h
)

add_library(coupler ${SOURCES})
set_target_properties(coupler PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_link_libraries(coupler adios2::adios2 Kokkos::kokkos PkgConfig::fftw)
target_include_directories(coupler PUBLIC ./)
