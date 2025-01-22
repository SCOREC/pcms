//
// Created by hasanm4 on 1/17/25.
//

#include <catch2/catch_test_macros.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/interpolator/interpolation_base.h>

#include <vector>
#include <iostream>

bool isClose(Omega_h::HostWrite<Omega_h::Real>& array1,
             Omega_h::HostWrite<Omega_h::Real>& array2,
             double percent_diff = 0.1);

/**
 * @brief Create sin(x)cos(y) + 2 at each vertex of the mesh
 */
void createSinxCosyAtVertex(Omega_h::Mesh& mesh,
                            Omega_h::Write<Omega_h::Real>& sinxcosy)
{
  OMEGA_H_CHECK_PRINTF(sinxcosy.size() == mesh.nverts(),
                       "Given ScalarArrayView is not properly sized %d vs %d\n",
                       mesh.nverts(), sinxcosy.size());

  // get the bounding box of the mesh
  Omega_h::BBox<2> bb = Omega_h::get_bounding_box<2>(&mesh);
  double dx = bb.max[0] - bb.min[0];
  double dy = bb.max[1] - bb.min[1];

  const auto& coords = mesh.coords();

  auto assignSinCos = OMEGA_H_LAMBDA(int node)
  {
    auto x = (coords[2 * node + 0] - bb.min[0]) / dx * 2.0 * M_PI;
    auto y = (coords[2 * node + 1] - bb.min[1]) / dy * 2.0 * M_PI;
    sinxcosy[node] = 2 + Kokkos::sin(x) * Kokkos::cos(y);
  };
  Omega_h::parallel_for(mesh.nverts(), assignSinCos, "assignSinCos");
}

void node2CentroidInterpolation(Omega_h::Mesh& mesh,
                                Omega_h::Write<Omega_h::Real>& node_vals,
                                Omega_h::Write<Omega_h::Real>& centroid_vals)
{
  OMEGA_H_CHECK_PRINTF(node_vals.size() == mesh.nverts(),
                       "Given ScalarArrayView is not properly sized %d vs %d\n",
                       mesh.nverts(), node_vals.size());
  OMEGA_H_CHECK_PRINTF(centroid_vals.size() == mesh.nfaces(),
                       "Given ScalarArrayView is not properly sized %d vs %d\n",
                       mesh.nfaces(), centroid_vals.size());

  const auto& face2node = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  auto coords = mesh.coords();
  Omega_h::LO nfaces = mesh.nfaces();

  auto averageSinCos = OMEGA_H_LAMBDA(Omega_h::LO face)
  {
    auto faceNodes = Omega_h::gather_verts<3>(face2node, face);
    Omega_h::Real sum = node_vals[faceNodes[0]] + node_vals[faceNodes[1]] +
                        node_vals[faceNodes[2]];
    centroid_vals[face] = sum / 3.0;
  };
  Omega_h::parallel_for(nfaces, averageSinCos, "averageSinCos");
}

TEST_CASE("Test MLSInterpolationHandler: Single Mesh")
{
  auto lib = Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  printf("[INFO] Mesh created with %d vertices and %d faces\n",
         source_mesh.nverts(), source_mesh.nfaces());

  auto mls_single = MLSInterpolationHandler(source_mesh, 0.12, 12, 3, true);

  Omega_h::Write<Omega_h::Real> sinxcosy_node(source_mesh.nverts(),
                                              "sinxcosy_node");
  createSinxCosyAtVertex(source_mesh, sinxcosy_node);
  Omega_h::Write<Omega_h::Real> sinxcosy_centroid(source_mesh.nfaces(),
                                                  "sinxcosy_centroid");
  node2CentroidInterpolation(source_mesh, sinxcosy_node, sinxcosy_centroid);

  Omega_h::HostWrite<double> source_data_host_write(sinxcosy_centroid);
  Omega_h::HostWrite<double> interpolated_data_hwrite(source_mesh.nverts());
  Omega_h::HostWrite<double> exact_values_at_nodes(sinxcosy_node);

  pcms::ScalarArrayView<double, pcms::HostMemorySpace> sourceArrayView(
    source_data_host_write.data(), source_data_host_write.size());
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> interpolatedArrayView(
    interpolated_data_hwrite.data(), interpolated_data_hwrite.size());

  OMEGA_H_CHECK_PRINTF(sourceArrayView.size() == mls_single.getSourceSize(),
                       "Source size mismatch: %zu vs %zu\n",
                       sourceArrayView.size(), mls_single.getSourceSize());
  OMEGA_H_CHECK_PRINTF(
    interpolatedArrayView.size() == mls_single.getTargetSize(),
    "Target size mismatch: %zu vs %zu\n", interpolatedArrayView.size(),
    mls_single.getTargetSize());

  mls_single.eval(sourceArrayView, interpolatedArrayView);

  REQUIRE(isClose(exact_values_at_nodes, interpolated_data_hwrite, 10.0) ==
          true);
}

bool isClose(Omega_h::HostWrite<Omega_h::Real>& array1,
             Omega_h::HostWrite<Omega_h::Real>& array2, double percent_diff)
{
  if (array1.size() != array2.size()) {
    fprintf(stderr, "[ERROR] Arrays are not of the same size: %d vs %d\n",
            array1.size(), array2.size());
    return false;
  }

  double eps = percent_diff / 100.0;
  for (int i = 0; i < array1.size(); i++) {
    if (std::abs(array1[i] - array2[i]) > eps * std::abs(array2[i])) {
      fprintf(stderr, "[ERROR] Arrays differ at index %d: %.16f vs %.16f\n", i,
              array1[i], array2[i]);
      return false;
    }
  }

  return true;
}