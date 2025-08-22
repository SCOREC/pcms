#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <pcms/transfer_field2.h>
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <Kokkos_Core.hpp>

void test_copy(Omega_h::CommPtr world, int dim, std::array<int, 4> nodes_per_dim,
               int num_components)
{
  int nx = 100;
  int ny = dim > 1 ? 100 : 0;
  int nz = dim > 2 ? 100 : 0;

  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, nx, ny, nz, false);
  pcms::OmegaHFieldLayout layout(mesh, nodes_per_dim, num_components, pcms::CoordinateSystem::Cartesian);
  int ndata = layout.GetNumOwnedDofHolder() * num_components;
  Omega_h::Write<double> ids(ndata);
  Omega_h::parallel_for(
    ndata, OMEGA_H_LAMBDA(int i) { ids[i] = i; });

  pcms::OmegaHField2 original("",  layout, mesh);
  pcms::Rank1View<const double, pcms::HostMemorySpace> array_view{
    std::data(ids), std::size(ids)};
  pcms::FieldDataView<const double, pcms::HostMemorySpace> field_data_view{
    array_view, original.GetCoordinateSystem()};
  original.SetDOFHolderData(field_data_view);

  pcms::OmegaHField2 copied("", layout, mesh);
  pcms::copy_field2(original, copied);
  auto copied_array = copied.GetDOFHolderData().GetValues();

  REQUIRE(copied_array.size() == ndata);
  int sum = 0;
  Kokkos::parallel_reduce(
    ndata,
    KOKKOS_LAMBDA(int i, int& local_sum) {
      local_sum += std::abs(ids[i] - copied_array[i]) < 1e-12;
    },
    sum);
  REQUIRE(sum == ndata);
}

TEST_CASE("copy omega_h_field2 data")
{
  auto lib = Omega_h::Library{};
  test_copy(lib.world(), 2, {1, 0, 0, 0}, 1);
  test_copy(lib.world(), 2, {1, 1, 0, 0}, 1);
}
