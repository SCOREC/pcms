#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include "pcms/transfer/copy.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/assert.h"
#include <Kokkos_Core.hpp>

using pcms::Real;

void test_copy(Omega_h::CommPtr world, int dim, int order, int num_components)
{
  int nx = 100;
  int ny = dim > 1 ? 100 : 0;
  int nz = dim > 2 ? 100 : 0;

  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, nx, ny, nz, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, order, num_components, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto layout = factory.GetLayout();
  int ndata = layout->GetNumOwnedDofHolder() * num_components;
  Omega_h::HostWrite<Real> ids(ndata);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, ndata),
    [=](int i) { ids[i] = i; });

  auto original = factory.CreateField<Real>(pcms::FieldMetadata{});
  original.SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(
      ids.data(), layout->GetNumOwnedDofHolder(), num_components));

  auto copied = factory.CreateField<Real>(pcms::FieldMetadata{});
  pcms::Copy<Real> copy(factory, factory);
  copy.Apply(original, copied);
  auto copied_array = pcms::FlattenToRank1View(copied.GetDOFHolderDataHost());

  REQUIRE(copied_array.size() == ndata);
  int sum = 0;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, ndata),
    [=](int i, int& local_sum) {
      local_sum += std::abs(ids[i] - copied_array[i]) < 1e-12;
    },
    sum);
  REQUIRE(sum == ndata);
}

TEST_CASE("copy omega_h_field2 data")
{
  auto lib = Omega_h::Library{};
  test_copy(lib.world(), 2, 0, 1);
  test_copy(lib.world(), 2, 1, 1);
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 100,
                                 100, 0, false);
  REQUIRE_THROWS_AS(pcms::LagrangeFunctionSpace::FromMesh(
                      mesh, 2, 1, pcms::CoordinateSystem::Cartesian, "global",
                      pcms::LagrangeFunctionSpace::Backend::OmegaH),
                    pcms::pcms_error);
}
