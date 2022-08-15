#include <wdmcpl/transfer_field.h>
#include <wdmcpl/omega_h_field.h>
#include <catch2/catch.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Kokkos_Core.hpp>

TEST_CASE("field copy", "[field transfer]")
{
  Omega_h::Library lib;
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  wdmcpl::OmegaHField<wdmcpl::LO> f1("source", mesh);
  Omega_h::Write<int> data(mesh.nents(0));
  Omega_h::parallel_for(
    data.size(), OMEGA_H_LAMBDA(int i) { data[i] = i; });
  mesh.add_tag<wdmcpl::LO>(0, "source", 1, data);
  wdmcpl::OmegaHField<wdmcpl::LO> f2("target", mesh);
  wdmcpl::copy_field(f1, f2);
  auto target_array = mesh.get_array<int>(0, "target");
  REQUIRE(target_array.size() == mesh.nents(0));
  //int result = 0;
  //Kokkos::parallel_reduce(target_array.size(),KOKKOS_LAMBDA(int i, int& lsum){lsum+=target_array[i];},result);
  // FIXME use parallel reduction rather than copy to host (currently linker errors when using kokkos)
  auto host_target_array = Omega_h::HostRead(target_array);
  auto result = std::accumulate(&host_target_array[0],&host_target_array[host_target_array.size()],0);
  auto n = target_array.size()-1;
  REQUIRE(result == n*(n+1)/2);
}

TEST_CASE("field lagrange projection", "[.][field transfer]")
{
  Omega_h::Library lib;
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  wdmcpl::OmegaHField<wdmcpl::LO> f1("source", mesh);
  wdmcpl::OmegaHField<wdmcpl::LO> f2("target", mesh);
  wdmcpl::project_field(f1, f2, wdmcpl::Lagrange<1>{});
}
