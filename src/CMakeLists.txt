set(SOURCES 
  adios2Routines.cc
  commpart1.cc
  importpart3mesh.cc
  BoundaryDescr3D.cc
  commpart1.cc
  fourierdataproc.cc
  dataprocess.cc 
  interpo.cc
  boundbuffer.cc
)

set(HEADERS
    adios2Routines.h
    couplingConstants.h
    couplingTypes.h
    commpart1.h
    dataprocess.h
    importpart3mesh.h
    BoundaryDescr3D.h
    testutilities.h
)

add_library(coupler ${SOURCES})
set_target_properties(coupler PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_link_libraries(coupler adios2::adios2 Kokkos::kokkos PkgConfig::fftw)
target_include_directories(coupler PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS coupler)
