#include <wdmcpl/transfer_field.h>
#include <wdmcpl/omega_h_field.h>
#include <catch2/catch_test_macros.hpp>
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
  int result = 0;
  Kokkos::parallel_reduce(
    target_array.size(),
    KOKKOS_LAMBDA(int i, int& lsum) { lsum += target_array[i]; }, result);
  auto n = target_array.size() - 1;
  REQUIRE(result == n * (n + 1) / 2);
}
static int sum_array(const Omega_h::Read<int>& target_array)
{
  int result = 0;
  Kokkos::parallel_reduce(
    target_array.size(),
    KOKKOS_LAMBDA(int i, int& lsum) { lsum += target_array[i]; }, result);
  return result;
}

TEST_CASE("field interpolation (identical fields)", "[field transfer]")
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

  SECTION("Nearest Neighbor")
  {
    wdmcpl::interpolate_field(f1, f2, wdmcpl::NearestNeighbor{});
    auto target_array = mesh.get_array<int>(0, "target");
    REQUIRE(target_array.size() == mesh.nents(0));
    int result = sum_array(target_array);
    auto n = target_array.size() - 1;
    REQUIRE(result == n * (n + 1) / 2);
  }
  SECTION("Lagrange<1>")
  {
    wdmcpl::interpolate_field(f1, f2, wdmcpl::Lagrange<1>{});
    auto target_array = mesh.get_array<int>(0, "target");
    REQUIRE(target_array.size() == mesh.nents(0));
    int result = sum_array(target_array);
    auto n = target_array.size() - 1;
    REQUIRE(result == n * (n + 1) / 2);
  }
}
