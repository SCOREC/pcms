//
// Created by hasanm4 on 1/17/25.
//

#include <catch2/catch_test_macros.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/interpolator/interpolation_base.h>

#include <vector>
#include <iostream>

void create_sinx_cosy_data(
  Omega_h::Mesh& mesh,
  pcms::ScalarArrayView<double, pcms::HostMemorySpace>& sinxcosy)
{
  // get the bounding box of the mesh
  Omega_h::BBox<2> bb = Omega_h::get_bounding_box<2>(&mesh);
  double dx = bb.max[0] - bb.min[0];
  double dy = bb.max[1] - bb.min[1];

  auto coords = mesh.coords();
  auto nnodes = mesh.nverts();
  Omega_h::Write<Omega_h::Real> sinxcosytag(nnodes);

  auto assignSinCos = OMEGA_H_LAMBDA(int node)
  {
    auto x = (coords[2 * node + 0] - bb.min[0]) / dx * 2.0 * M_PI;
    auto y = (coords[2 * node + 1] - bb.min[1]) / dy * 2.0 * M_PI;
    sinxcosytag[node] = 2 + sin(x) * cos(y);
  };
  Omega_h::parallel_for(nnodes, assignSinCos, "assignSinCos");
  Omega_h::HostWrite<Omega_h::Real> sinxcosyhost(sinxcosytag);

  OMEGA_H_CHECK_PRINTF(
    nnodes = sinxcosy.size(),
    "Given ScalarArrayView is not properly sized %d vs %zu\n", nnodes,
    sinxcosy.size());

  for (int i = 0; i < nnodes; ++i) {
    sinxcosy[i] = sinxcosyhost[i];
  }
}

TEST_CASE("Test MLSInterpolationHandler: Single Mesh")
{
  auto lib = Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  printf("[INFO] Mesh created with %d vertices and %d faces\n", source_mesh.nverts(),
         source_mesh.nfaces());

  auto mls_single = MLSInterpolationHandler(source_mesh, 0.09, true);

  Omega_h::HostWrite<double> source_data_host_write(source_mesh.nverts());
  Omega_h::HostWrite<double> interpolated_data_hwrite(source_mesh.nfaces());
  Omega_h::HostWrite<double> target_data_expected(source_mesh.nfaces());

  pcms::ScalarArrayView<double, pcms::HostMemorySpace> sourceArrayView(
    source_data_host_write.data(), source_data_host_write.size());
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> interpolatedArrayView(
    interpolated_data_hwrite.data(), interpolated_data_hwrite.size());
  pcms::ScalarArrayView<double, pcms::HostMemorySpace>
    sinxcosyarrayview_target_expected(target_data_expected.data(),
                                      target_data_expected.size());

  create_sinx_cosy_data(source_mesh, sourceArrayView);

  mls_single.eval(sourceArrayView, interpolatedArrayView);
}