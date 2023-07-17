#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <wdmcpl/transfer_field.h>
#include <wdmcpl/omega_h_field.h>
#include <Kokkos_Core.hpp>

TEST_CASE("copy omega_h_field data")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 100, 100, 0, false);
  const auto nverts = mesh.nents(0);
  Omega_h::Write<int> ids(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) { ids[i] = i; });
  mesh.add_tag<int>(0,"test_ids",1,Omega_h::Read(ids));
  const bool tag_already_exists = GENERATE(true,false);
  if(tag_already_exists) {
    Omega_h::Write<int> zeros(nverts, 0);
    mesh.add_tag<int>(0,"copied",1,zeros);
  }
  SECTION("No filter") {
    wdmcpl::OmegaHField<int,double> original("test_ids",mesh);
    wdmcpl::OmegaHField<int,double> copied("copied",mesh);
    REQUIRE(original.Size() == copied.Size());
    wdmcpl::copy_field(original,copied);
    auto copied_array = wdmcpl::get_nodal_data(copied);
    REQUIRE(copied_array.size() == copied.Size());
    REQUIRE(copied_array.size() == nverts);
    int sum=0;
    Kokkos::parallel_reduce(nverts, KOKKOS_LAMBDA(int i, int &local_sum) {
                                    local_sum += (ids[i] == copied_array[i]);
                                    },sum);
    REQUIRE(sum == nverts);
  }
  SECTION("trivial positive mask"){
    Omega_h::Write<Omega_h::I8> mask(nverts,1);
    wdmcpl::OmegaHField<int,double> original("test_ids",mesh, mask);
    wdmcpl::OmegaHField<int,double> copied("copied",mesh,mask);
    REQUIRE(original.Size() == copied.Size());
    wdmcpl::copy_field(original,copied);
    auto copied_array = wdmcpl::get_nodal_data(copied);
    REQUIRE(copied_array.size() == copied.Size());
    REQUIRE(copied_array.size() == nverts);
    int sum=0;
    Kokkos::parallel_reduce(nverts, KOKKOS_LAMBDA(int i, int &local_sum) {
      local_sum += (ids[i] == copied_array[i]);
    },sum);
    REQUIRE(sum == nverts);
  }
  SECTION("every other mask"){
    Omega_h::Write<Omega_h::I8> mask(nverts,0);
    Omega_h::parallel_for(
      nverts, OMEGA_H_LAMBDA(int i) { mask[i] = i%2; });
    wdmcpl::OmegaHField<int,double> original("test_ids",mesh,mask);
    wdmcpl::OmegaHField<int,double> copied("copied",mesh,mask);
    REQUIRE(original.Size() == copied.Size());
    wdmcpl::copy_field(original,copied);
    auto copied_array = wdmcpl::get_nodal_data(copied);
    auto original_array = wdmcpl::get_nodal_data(original);
    REQUIRE(copied_array.size() == copied.Size());
    REQUIRE(original_array.size() == original.Size());
    int sum=0;
    Kokkos::parallel_reduce(original_array.size(), KOKKOS_LAMBDA(int i, int &local_sum) {
      local_sum += (original_array[i] == copied_array[i]);
    },sum);
    REQUIRE(sum == original_array.size());
  }
}