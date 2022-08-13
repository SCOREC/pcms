#include <wdmcpl/transfer_field.h>
#include <wdmcpl/omega_h_field.h>
#include <catch2/catch.hpp>
#include <Omega_h_mesh.hpp>


TEST_CASE( "field copy", "[field transfer]" ) {
  Omega_h::Library lib;
  Omega_h::Mesh mesh(&lib);

  //MPI_Comm mpi_comm = lib.world()->get_impl();
 wdmcpl::OmegaHField f1("source", mesh);
 wdmcpl::OmegaHField f2("target", mesh);
 wdmcpl::copy_field(f1,f2);
}




