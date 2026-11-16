//
// Created by hasanm4 on 1/17/25.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/transfer/interpolation_base.h>
#include <pcms/utility/mesh_geometry.h>
#include <pcms/utility/print.h>
#include "field_test_utils.h"

#include <vector>
#include <iostream>

void translate_mesh(Omega_h::Mesh* mesh, Omega_h::Vector<2> translation_vector)
{
  auto coords = mesh->coords();
  auto nverts = mesh->nverts();
  auto out = Omega_h::Write<Omega_h::Real>(coords.size());

  auto f = OMEGA_H_LAMBDA(Omega_h::LO i)
  {
    auto coord = Omega_h::get_vector<2>(coords, i);
    coord = coord + translation_vector;
    Omega_h::set_vector<2>(out, i, coord);
  };

  Omega_h::parallel_for(nverts, f);
  mesh->set_coords(Omega_h::Reals(out));
}

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

TEST_CASE("Test MLSMeshInterpolation")
{
  pcms::printInfo("[INFO] Starting MLS Interpolation Test...\n");
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 20, 20, 0, false);
  pcms::printInfo("[INFO] Mesh created with %d vertices and %d faces\n",
                  source_mesh.nverts(), source_mesh.nfaces());

  Omega_h::Write<Omega_h::Real> source_sinxcosy_node(source_mesh.nverts(),
                                                     "source_sinxcosy_node");
  createSinxCosyAtVertex(source_mesh, source_sinxcosy_node);

  SECTION("Single Mesh")
  {
    pcms::printInfo("\n-------------------- Single Mesh Interpolation Test "
                    "Started --------------------\n");
    pcms::printInfo("Mesh based search...\n");
    auto mls_single =
      pcms::MLSMeshInterpolation(source_mesh, 0.12, 15, 3, true, 0.0, 5.0);

    auto source_points_vec = pcms::test::CopyOmegaHRealsToVector(
      pcms::get_entity_centroids(source_mesh, Omega_h::FACE));
    auto source_points_view = pcms::make_array_view(source_points_vec);

    auto target_points_vec =
      pcms::test::CopyOmegaHRealsToVector(source_mesh.coords());
    auto target_points_view = pcms::make_array_view(target_points_vec);
    REQUIRE(source_mesh.dim() == 2);
    pcms::printInfo("Point cloud based search...\n");
    auto point_mls = pcms::MLSPointCloudInterpolation(
      source_points_view, target_points_view, 2, 0.12, 15, 3, true, 0.0, 5.0);

    Omega_h::Write<Omega_h::Real> sinxcosy_centroid(source_mesh.nfaces(),
                                                    "sinxcosy_centroid");
    node2CentroidInterpolation(source_mesh, source_sinxcosy_node,
                               sinxcosy_centroid);

    Omega_h::HostWrite<double> source_data_host_write(sinxcosy_centroid);
    Omega_h::HostWrite<double> interpolated_data_hwrite(source_mesh.nverts());
    Omega_h::HostWrite<double> point_cloud_interpolated_data_hwrite(
      source_mesh.nverts());
    Omega_h::HostWrite<double> exact_values_at_nodes(source_sinxcosy_node);

    auto sourceArrayView = pcms::make_array_view(source_data_host_write);
    auto interpolatedArrayView =
      pcms::make_array_view(interpolated_data_hwrite);
    auto point_cloud_interpolatedArrayView =
      pcms::make_array_view(point_cloud_interpolated_data_hwrite);

    OMEGA_H_CHECK_PRINTF(sourceArrayView.size() == mls_single.getSourceSize(),
                         "Source size mismatch: %zu vs %zu\n",
                         sourceArrayView.size(), mls_single.getSourceSize());
    OMEGA_H_CHECK_PRINTF(
      interpolatedArrayView.size() == mls_single.getTargetSize(),
      "Target size mismatch: %zu vs %zu\n", interpolatedArrayView.size(),
      mls_single.getTargetSize());

    pcms::printInfo("Evaluating Mesh based MLS Interpolation...\n");
    mls_single.eval(sourceArrayView, interpolatedArrayView);
    pcms::printInfo("Evaluating Point Cloud Based MLS Interpolation...\n");
    point_mls.eval(sourceArrayView, point_cloud_interpolatedArrayView);

    // write the meshes in vtk format with the interpolated values as tag to
    // visualize the interpolated values
    source_mesh.add_tag<double>(
      Omega_h::VERT, "mesh_interpolated_sinxcosy", 1,
      Omega_h::Write<double>(interpolated_data_hwrite));
    source_mesh.add_tag<double>(
      Omega_h::VERT, "point_cloud_interpolated_sinxcosy", 1,
      Omega_h::Write(point_cloud_interpolated_data_hwrite));
    source_mesh.add_tag<double>(Omega_h::VERT, "exact_sinxcosy", 1,
                                Omega_h::Write(exact_values_at_nodes));
    source_mesh.add_tag<double>(Omega_h::FACE, "centroid_sinxcosy", 1,
                                sinxcosy_centroid);
    Omega_h::vtk::write_parallel("source_mesh_with_tags_to_debug.vtk",
                                 &source_mesh);

    REQUIRE(isClose(exact_values_at_nodes, interpolated_data_hwrite, 10.0) ==
            true);
    pcms::printInfo(
      "[****] Single Mesh Interpolation Test Passed with %.2f%% tolerance!\n",
      10.0);

    REQUIRE(isClose(exact_values_at_nodes, point_cloud_interpolated_data_hwrite,
                    10.0) == true);
    //////////////////////////////////// *************** Test Supports are same
    ///*****************************//
    auto mesh_based_supports = mls_single.getSupports();
    auto point_cloud_based_supports = point_mls.getSupports();
    pcms::test::CheckSupportResultsEquivalent(point_cloud_based_supports,
                                              mesh_based_supports);

    // Check if the point cloud interpolation is same as the MLS interpolation
    pcms::printDebugInfo("Interpolated data size: %d\n",
                         point_cloud_interpolated_data_hwrite.size());
    REQUIRE(point_cloud_interpolated_data_hwrite.size() ==
            interpolated_data_hwrite.size());

    for (int i = 0; i < interpolated_data_hwrite.size(); i++) {
      pcms::printDebugInfo("Interpolated data: %d, %.16f, %.16f\n", i,
                           interpolated_data_hwrite[i],
                           point_cloud_interpolated_data_hwrite[i]);

      REQUIRE_THAT(
        point_cloud_interpolated_data_hwrite[i],
        Catch::Matchers::WithinAbs(interpolated_data_hwrite[i], 1e-4));
    }
  }

  SECTION("Double Mesh")
  {
    pcms::printInfo("\n-------------------- Double Mesh Interpolation Test "
                    "Started --------------------\n");
    auto target_mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 17, 17, 0, false);
    pcms::printInfo(
      "[INFO] Target Mesh created with %d vertices and %d faces\n",
      target_mesh.nverts(), target_mesh.nfaces());

    // TODO: This is a way around.
    // https://github.com/SCOREC/pcms/pull/148#discussion_r1926204199
    // translate_mesh(&target_mesh, Omega_h::Vector<2>{(1.0 - 0.999) / 2.0,
    //                                                (1.0 - 0.999) / 2.0});

    auto mls_double = pcms::MLSMeshInterpolation(source_mesh, target_mesh, 0.12,
                                                 15, 3, true, 0.0, 5.0);

    Omega_h::HostWrite<double> source_data_host_write(source_sinxcosy_node);
    Omega_h::HostWrite<double> interpolated_data_hwrite(
      mls_double.getTargetSize());

    auto sourceArrayView = pcms::make_array_view(source_data_host_write);
    auto interpolatedArrayView =
      pcms::make_array_view(interpolated_data_hwrite);

    mls_double.eval(sourceArrayView, interpolatedArrayView);

    Omega_h::Write<Omega_h::Real> exact_target_sinxcosy_node(
      target_mesh.nverts(), "target_sinxcosy_node");
    createSinxCosyAtVertex(target_mesh, exact_target_sinxcosy_node);
    Omega_h::HostWrite<double> exact_target_sinxcosy_node_hwrite(
      exact_target_sinxcosy_node);

    REQUIRE(isClose(interpolated_data_hwrite, exact_target_sinxcosy_node_hwrite,
                    10.0) == true);

    pcms::printInfo(
      "[INFO] Double Mesh Interpolation Test Passed with %.2f%% tolerance!\n",
      10.0);
  }
}

TEST_CASE("MLSPointCloudInterpolation honors provided dimension in eval")
{
  auto source_points = Omega_h::HostWrite<Omega_h::Real>(27 * 3);
  int idx = 0;
  for (int z = 0; z < 3; ++z) {
    for (int y = 0; y < 3; ++y) {
      for (int x = 0; x < 3; ++x) {
        source_points[idx++] = x;
        source_points[idx++] = y;
        source_points[idx++] = z;
      }
    }
  }
  auto target_points = source_points;

  auto source_points_view = pcms::make_array_view(source_points);
  auto target_points_view = pcms::make_array_view(target_points);

  // Degree-1 polynomial that depends on z to catch accidental 2D behavior.
  auto source_values = Omega_h::HostWrite<Omega_h::Real>(27);
  for (int i = 0; i < 27; ++i) {
    const auto x = source_points[3 * i + 0];
    const auto y = source_points[3 * i + 1];
    const auto z = source_points[3 * i + 2];
    source_values[i] = 1.0 + 2.0 * x - 3.0 * y + 5.0 * z;
  }

  auto output_values = Omega_h::HostWrite<Omega_h::Real>(27, "output_values");
  auto source_values_view = pcms::make_array_view(source_values);
  auto output_values_view = pcms::make_array_view(output_values);

  auto mls = pcms::MLSPointCloudInterpolation(
    source_points_view, target_points_view, 3, 2.5, 10, 1, true, 0.0, 5.0);

  REQUIRE_NOTHROW(mls.eval(source_values_view, output_values_view));
  for (int i = 0; i < output_values.size(); ++i) {
    CHECK_THAT(output_values[i],
               Catch::Matchers::WithinAbs(source_values[i], 1e-6));
  }
}

bool isClose(Omega_h::HostWrite<Omega_h::Real>& array1,
             Omega_h::HostWrite<Omega_h::Real>& array2, double percent_diff)
{
  if (array1.size() != array2.size()) {
    pcms::printError("[ERROR] Arrays are not of the same size: %d vs %d\n",
                     array1.size(), array2.size());
    return false;
  }

  double eps = percent_diff / 100.0;
  for (int i = 0; i < array1.size(); i++) {
    if (std::abs(array1[i] - array2[i]) > eps * std::abs(array2[i])) {
      pcms::printError("[ERROR] Arrays differ at index %d: %.16f vs %.16f\n", i,
                       array1[i], array2[i]);
      return false;
    }
  }

  return true;
}
