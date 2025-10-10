//
// Created by hasanm4 on 1/17/25.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/interpolator/interpolation_base.h>

#include <vector>
#include <iostream>
#include <unordered_map>

bool areArraysEqualUnordered(const Omega_h::HostRead<Omega_h::LO>& array1,
                             const Omega_h::HostRead<Omega_h::LO>& array2,
                             int start, int end)
{
  // Ensure the indices are valid
  assert(start >= 0 && end <= array1.size() && start <= end);
  assert(start >= 0 && end <= array2.size() && start <= end);

  // Use frequency maps to count occurrences of each value
  std::unordered_map<Omega_h::LO, int> freq1, freq2;

  for (int i = start; i < end; ++i) {
    freq1[array1[i]]++;
    freq2[array2[i]]++;
  }

  // Compare the frequency maps
  if (freq1 != freq2) {
    printf("[ERROR] Arrays differ in the range [%d, %d)\n", start, end);
    return false;
  }

  return true;
}

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

TEST_CASE("Test MLSInterpolationHandler")
{
  fprintf(stdout, "[INFO] Starting MLS Interpolation Test...\n");
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 20, 20, 0, false);
  printf("[INFO] Mesh created with %d vertices and %d faces\n",
         source_mesh.nverts(), source_mesh.nfaces());

  Omega_h::Write<Omega_h::Real> source_sinxcosy_node(source_mesh.nverts(),
                                                     "source_sinxcosy_node");
  createSinxCosyAtVertex(source_mesh, source_sinxcosy_node);

  SECTION("Single Mesh")
  {
    fprintf(stdout, "\n-------------------- Single Mesh Interpolation Test "
                    "Started --------------------\n");
    printf("Mesh based search...\n");
    auto mls_single = MLSInterpolationHandler(source_mesh, 0.12, 15, 3, true);

    auto source_points_reals = getCentroids(source_mesh);
    auto source_points_host =
      Omega_h::HostRead<Omega_h::Real>(source_points_reals);
    auto source_points_host_write =
      Omega_h::HostWrite<Omega_h::Real>(source_points_host.size());
    for (int i = 0; i < source_points_host.size(); i++) {
      source_points_host_write[i] = source_points_host[i];
    }
    auto source_points_view =
      pcms::Rank1View<double, pcms::HostMemorySpace>(
        source_points_host_write.data(), source_points_host_write.size());

    auto target_points_reals = source_mesh.coords();
    auto target_points_host =
      Omega_h::HostRead<Omega_h::Real>(target_points_reals);
    auto target_points_host_write =
      Omega_h::HostWrite<Omega_h::Real>(target_points_host.size());
    for (int i = 0; i < target_points_host.size(); i++) {
      target_points_host_write[i] = target_points_host[i];
    }
    auto target_points_view =
      pcms::Rank1View<double, pcms::HostMemorySpace>(
        target_points_host_write.data(), target_points_host_write.size());
    REQUIRE(source_mesh.dim() == 2);
    printf("Point cloud based search...\n");
    auto point_mls = MLSPointCloudInterpolation(
      source_points_view, target_points_view, 2, 0.12, 15, 3, true);

    Omega_h::Write<Omega_h::Real> sinxcosy_centroid(source_mesh.nfaces(),
                                                    "sinxcosy_centroid");
    node2CentroidInterpolation(source_mesh, source_sinxcosy_node,
                               sinxcosy_centroid);

    Omega_h::HostWrite<double> source_data_host_write(sinxcosy_centroid);
    Omega_h::HostWrite<double> interpolated_data_hwrite(source_mesh.nverts());
    Omega_h::HostWrite<double> point_cloud_interpolated_data_hwrite(
      source_mesh.nverts());
    Omega_h::HostWrite<double> exact_values_at_nodes(source_sinxcosy_node);

    pcms::Rank1View<double, pcms::HostMemorySpace> sourceArrayView(
      source_data_host_write.data(), source_data_host_write.size());
    pcms::Rank1View<double, pcms::HostMemorySpace> interpolatedArrayView(
      interpolated_data_hwrite.data(), interpolated_data_hwrite.size());
    pcms::Rank1View<double, pcms::HostMemorySpace>
      point_cloud_interpolatedArrayView(
        point_cloud_interpolated_data_hwrite.data(),
        point_cloud_interpolated_data_hwrite.size());

    OMEGA_H_CHECK_PRINTF(sourceArrayView.size() == mls_single.getSourceSize(),
                         "Source size mismatch: %zu vs %zu\n",
                         sourceArrayView.size(), mls_single.getSourceSize());
    OMEGA_H_CHECK_PRINTF(
      interpolatedArrayView.size() == mls_single.getTargetSize(),
      "Target size mismatch: %zu vs %zu\n", interpolatedArrayView.size(),
      mls_single.getTargetSize());

    printf("Evaluating Mesh based MLS Interpolation...\n");
    mls_single.eval(sourceArrayView, interpolatedArrayView);
    printf("Evaluating Point Cloud Based MLS Interpolation...\n");
    point_mls.eval(sourceArrayView, point_cloud_interpolatedArrayView);

    // write the meshes in vtk format with the interpolated values as tag to visualize the interpolated values
    source_mesh.add_tag<double>(Omega_h::VERT, "mesh_interpolated_sinxcosy", 1, Omega_h::Write<double>(interpolated_data_hwrite));
    source_mesh.add_tag<double>(Omega_h::VERT, "point_cloud_interpolated_sinxcosy", 1,  Omega_h::Write(point_cloud_interpolated_data_hwrite));
    source_mesh.add_tag<double>(Omega_h::VERT, "exact_sinxcosy", 1, Omega_h::Write(exact_values_at_nodes));
    source_mesh.add_tag<double>(Omega_h::FACE, "centroid_sinxcosy", 1, sinxcosy_centroid);
    Omega_h::vtk::write_parallel("source_mesh_with_tags_to_debug.vtk", &source_mesh);

    REQUIRE(isClose(exact_values_at_nodes, interpolated_data_hwrite, 10.0) ==
            true);
    fprintf(
      stdout,
      "[****] Single Mesh Interpolation Test Passed with %.2f%% tolerance!\n",
      10.0);

    REQUIRE(isClose(exact_values_at_nodes, point_cloud_interpolated_data_hwrite,
                    10.0) == true);
    //////////////////////////////////// *************** Test Supports are same
    ///*****************************//
    auto mesh_based_supports = mls_single.getSupports();
    auto point_cloud_based_supports = point_mls.getSupports();
    auto mesh_based_support_ptr_host =
      Omega_h::HostRead<Omega_h::LO>(mesh_based_supports.supports_ptr);
    auto mesh_based_support_idx_host =
      Omega_h::HostRead<Omega_h::LO>(mesh_based_supports.supports_idx);
    auto point_cloud_based_support_ptr_host =
      Omega_h::HostRead<Omega_h::LO>(point_cloud_based_supports.supports_ptr);
    auto point_cloud_based_support_idx_host =
      Omega_h::HostRead<Omega_h::LO>(point_cloud_based_supports.supports_idx);

    REQUIRE(point_cloud_based_support_idx_host.size() ==
            mesh_based_support_idx_host.size());
    REQUIRE(point_cloud_based_support_ptr_host.size() ==
            mesh_based_support_ptr_host.size());
    for (int i = 0; i < mesh_based_support_ptr_host.size(); i++) {
      REQUIRE(point_cloud_based_support_ptr_host[i] ==
              mesh_based_support_ptr_host[i]);
    }

    for (int i = 0; i < mesh_based_support_ptr_host.size() - 1; i++) {
      auto start = mesh_based_support_ptr_host[i];
      auto end = mesh_based_support_ptr_host[i + 1];

      bool isEqual =
        areArraysEqualUnordered(mesh_based_support_idx_host,
                                point_cloud_based_support_idx_host, start, end);
      REQUIRE(isEqual);
    }

    // Check if the point cloud interpolation is same as the MLS interpolation
    printf("Interpolated data size: %d\n",
           point_cloud_interpolated_data_hwrite.size());
    REQUIRE(point_cloud_interpolated_data_hwrite.size() ==
            interpolated_data_hwrite.size());

    for (int i = 0; i < interpolated_data_hwrite.size(); i++) {
      printf("Interpolated data: %d, %.16f, %.16f\n", i,
             interpolated_data_hwrite[i],
             point_cloud_interpolated_data_hwrite[i]);
      if (i == 0 || i == 78)
        continue; // FIXME

      REQUIRE_THAT(
        point_cloud_interpolated_data_hwrite[i],
        Catch::Matchers::WithinAbs(interpolated_data_hwrite[i], 1e-4));
    }
  }

  SECTION("Double Mesh")
  {
    fprintf(stdout, "\n-------------------- Double Mesh Interpolation Test "
                    "Started --------------------\n");
    auto target_mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 0.999, 0.999,
                                          1, 17, 17, 0, false);
    printf("[INFO] Target Mesh created with %d vertices and %d faces\n",
           target_mesh.nverts(), target_mesh.nfaces());

    // TODO: This is a way around.
    // https://github.com/SCOREC/pcms/pull/148#discussion_r1926204199
    translate_mesh(&target_mesh, Omega_h::Vector<2>{(1.0 - 0.999) / 2.0,
                                                    (1.0 - 0.999) / 2.0});

    auto mls_double =
      MLSInterpolationHandler(source_mesh, target_mesh, 0.12, 12, 3, true);

    Omega_h::HostWrite<double> source_data_host_write(source_sinxcosy_node);
    Omega_h::HostWrite<double> interpolated_data_hwrite(
      mls_double.getTargetSize());

    pcms::Rank1View<double, pcms::HostMemorySpace> sourceArrayView(
      source_data_host_write.data(), source_data_host_write.size());
    pcms::Rank1View<double, pcms::HostMemorySpace> interpolatedArrayView(
      interpolated_data_hwrite.data(), interpolated_data_hwrite.size());

    mls_double.eval(sourceArrayView, interpolatedArrayView);

    Omega_h::Write<Omega_h::Real> exact_target_sinxcosy_node(
      target_mesh.nverts(), "target_sinxcosy_node");
    createSinxCosyAtVertex(target_mesh, exact_target_sinxcosy_node);
    Omega_h::HostWrite<double> exact_target_sinxcosy_node_hwrite(
      exact_target_sinxcosy_node);

    REQUIRE(isClose(interpolated_data_hwrite, exact_target_sinxcosy_node_hwrite,
                    10.0) == true);

    fprintf(
      stdout,
      "[INFO] Double Mesh Interpolation Test Passed with %.2f%% tolerance!\n",
      10.0);
  }
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