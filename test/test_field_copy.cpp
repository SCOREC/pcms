#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <pcms/transfer_field2.h>
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/create_field.h"
#include <Kokkos_Core.hpp>

void test_copy(Omega_h::CommPtr world, int dim, int order, int num_components)
{
  int nx = 100;
  int ny = dim > 1 ? 100 : 0;
  int nz = dim > 2 ? 100 : 0;

  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, nx, ny, nz, false);
  auto layout = pcms::CreateLagrangeLayout(mesh, order, num_components,
                                           pcms::CoordinateSystem::Cartesian);
  int ndata = layout->GetNumOwnedDofHolder() * num_components;
  Omega_h::Write<double> ids(ndata);
  Omega_h::parallel_for(ndata, OMEGA_H_LAMBDA(int i) { ids[i] = i; });

  auto original = layout->CreateField();
  original->SetDOFHolderData(pcms::make_const_array_view(ids));

  auto copied = layout->CreateField();
  pcms::copy_field2(*original, *copied);
  auto copied_array = copied->GetDOFHolderData();

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
  test_copy(lib.world(), 2, 1, 1);
  test_copy(lib.world(), 2, 2, 1);
}
