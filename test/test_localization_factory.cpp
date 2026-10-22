#include <catch2/catch_test_macros.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>

#include "pcms/field/layout/omega_h_entity.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/layout/point_cloud.h"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/localization/adj_search.hpp"
#include "pcms/localization/localization_path_selection.h"
#include "pcms/localization/mesh_localization.h"
#include "pcms/localization/mls_support_helpers.h"
#include "pcms/localization/point_cloud_localization.h"
#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/omega_h_array_utils.h"
#include "field_test_utils.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"

namespace
{

std::shared_ptr<pcms::OmegaHEntityLayout> MakeMeshEntityLayout(
  Omega_h::Mesh& mesh, int entity_dim)
{
  return std::make_shared<pcms::OmegaHEntityLayout>(
    mesh, entity_dim, 1, pcms::csys::Cartesian::Deferred());
}

} // namespace

TEST_CASE(
  "LocalizationFactory: point-cloud Build matches BuildPointCloudSupports")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 8, 8, 0, false);

  pcms::MLSOptions options;
  options.radius = 0.2;
  options.min_req_supports = 10;
  options.adapt_radius = true;

  const int dim = mesh.dim();
  auto source_coords = mesh.coords();
  auto target_coords = pcms::test::CopyOmegaHRealsToVector(source_coords);
  auto target_device = pcms::test::CreateDeviceCoordinateView(
    target_coords, pcms::csys::Cartesian::Deferred(), dim);

  auto coords_dev = pcms::ConvertCoordsTo2D(source_coords, mesh.nverts(), dim);

  auto layout = std::make_shared<pcms::PointCloudLayout>(
    dim, coords_dev, pcms::csys::Cartesian::Deferred());
  pcms::PointCloudLocalizationFactory factory(layout, options);

  auto actual = factory.Build(target_device.coordinate_view);
  auto expected = pcms::BuildPointCloudSupports(
    source_coords, source_coords, dim, options.radius, options.min_req_supports,
    options.adapt_radius);

  pcms::test::CheckSupportResultsEquivalent(actual, expected);
}

TEST_CASE("LocalizationFactory: vertex adjacency Build matches two-mesh "
          "searchNeighbors")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 12, 12, 0, false);
  auto target_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  pcms::MLSOptions options;
  options.radius = 0.12;
  options.min_req_supports = 15;
  options.adapt_radius = true;

  pcms::AdjacencyLocalizationFactory factory(source_mesh, Omega_h::VERT,
                                             options);
  auto target_coords =
    pcms::test::CopyOmegaHRealsToVector(target_mesh.coords());
  auto target_device = pcms::test::CreateDeviceCoordinateView(
    target_coords, pcms::csys::Cartesian::Deferred(), target_mesh.dim());

  auto actual = factory.Build(target_device.coordinate_view);

  Omega_h::Real radius_sq = options.radius * options.radius;
  auto expected = pcms::searchNeighbors(
    source_mesh, target_mesh, radius_sq,
    static_cast<Omega_h::LO>(options.min_req_supports),
    static_cast<Omega_h::LO>(3 * options.min_req_supports),
    options.adapt_radius);

  pcms::test::CheckSupportResultsEquivalent(actual, expected);
}

TEST_CASE("Localization path selection: vertex source uses adjacency for "
          "pointwise requests")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 8, 8, 0, false);

  auto source_layout = std::make_shared<pcms::OmegaHLagrangeLayout>(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred());

  REQUIRE(pcms::detail::SelectLocalizationPath(*source_layout) ==
          pcms::detail::LocalizationPath::VertexAdjacencySearch);
}

TEST_CASE(
  "Localization path selection: centroid-to-vertex uses same-mesh adjacency")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 8, 8, 0, false);

  auto source_layout = MakeMeshEntityLayout(mesh, pcms::Face);
  auto target_layout = std::make_shared<pcms::OmegaHLagrangeLayout>(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred());

  REQUIRE(
    pcms::detail::SelectLocalizationPath(*source_layout, target_layout.get()) ==
    pcms::detail::LocalizationPath::CentroidToVertexAdjacencySearch);
}

TEST_CASE("Localization path selection: centroid source uses point-cloud "
          "supports for arbitrary points")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 8, 8, 0, false);

  auto source_layout = MakeMeshEntityLayout(mesh, pcms::Face);

  REQUIRE(pcms::detail::SelectLocalizationPath(*source_layout) ==
          pcms::detail::LocalizationPath::PointCloudSupports);
}

TEST_CASE("LocalizationFactory: centroid-to-vertex Build correct on same-mesh "
          "fast path")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 8, 8, 0, false);

  pcms::MLSOptions options;
  options.radius = 0.25;
  options.min_req_supports = 8;
  options.adapt_radius = true;

  // Source: face centroids on the mesh; target: mesh vertices.
  pcms::AdjacencyLocalizationFactory factory(mesh, Omega_h::FACE, options);
  auto actual = factory.BuildSameMeshCentroidToVertex();

  // Reference: direct call to the single-mesh centroid-to-vertex overload.
  Omega_h::Real radius_sq = options.radius * options.radius;
  auto expected = pcms::searchNeighbors(
    mesh, radius_sq, static_cast<Omega_h::LO>(options.min_req_supports),
    options.adapt_radius);

  pcms::test::CheckSupportResultsEquivalent(actual, expected);
}

TEST_CASE("LocalizationFactory: correct with different meshes")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 8, 8, 0, false);
  auto target_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 6, 6, 0, false);

  pcms::MLSOptions options;
  options.radius = 0.25;
  options.min_req_supports = 8;
  options.adapt_radius = true;

  pcms::AdjacencyLocalizationFactory factory(source_mesh, Omega_h::FACE,
                                             options);

  auto target_coords =
    pcms::test::CopyOmegaHRealsToVector(target_mesh.coords());
  auto target_device = pcms::test::CreateDeviceCoordinateView(
    target_coords, pcms::csys::Cartesian::Deferred(), target_mesh.dim());

  // Should not throw — falls back to point-cloud N^2 path.
  auto result = factory.Build(target_device.coordinate_view);
  auto ptr_host = Omega_h::HostRead<Omega_h::LO>(result.supports_ptr);
  REQUIRE(ptr_host.size() == target_mesh.nverts() + 1);
}
