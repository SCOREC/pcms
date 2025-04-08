#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <pcms/transfer_field2.h>
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <Kokkos_Core.hpp>

TEST_CASE("copy omega_h_field2 data")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 100, 100, 0, false);
  const auto nverts = mesh.nents(0);
  Omega_h::Write<double> ids(nverts);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) { ids[i] = i; });
  mesh.add_tag<double>(0, "test_ids", 1, Omega_h::Read(ids));
  const bool tag_already_exists = GENERATE(true, false);
  if (tag_already_exists) {
    Omega_h::Write<double> zeros(nverts, 0);
    mesh.add_tag<double>(0, "copied", 1, zeros);
  }
  auto layout =
    pcms::OmegaHFieldLayout(mesh, pcms::OmegaHFieldLayoutLocation::PieceWise, 2);
  pcms::OmegaHField2 original("test_ids", pcms::CoordinateSystem::Cartesian, layout, mesh);
  pcms::OmegaHField2 copied("copied", pcms::CoordinateSystem::Cartesian, layout, mesh);
  pcms::copy_field2(original, copied);
  auto copied_array = copied.GetDOFHolderData().GetValues();
  REQUIRE(copied_array.size() == nverts);
  int sum = 0;
  Kokkos::parallel_reduce(
    nverts,
    KOKKOS_LAMBDA(int i, int& local_sum) {
      local_sum += std::abs(ids[i] - copied_array[i]) < 1e-12;
    },
    sum);
  REQUIRE(sum == nverts);
}
