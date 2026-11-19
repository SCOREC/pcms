#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/backends/point_cloud/mls_evaluator.hpp>
#include <pcms/backends/point_cloud/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/omega_h/entity_layout.hpp>
#include <pcms/utility/entity_types.hpp>
#include <pcms/utility/mesh_geometry.hpp>
#include "field_test_utils.hpp"

#include <cmath>
#include <vector>

namespace pcms
{
using detail::MLSBinding;
using omega_h::OmegaHEntityLayout;
using utility::get_entity_centroids;
} // namespace pcms

using Catch::Matchers::WithinAbs;
using pcms::HostMemorySpace;
using pcms::Real;

namespace
{

/// 2 + sin(x) cos(y) over the unit square scaled to [0, 2 pi]^2
KOKKOS_INLINE_FUNCTION Real SinxCosy(Real x, Real y)
{
  return 2 + Kokkos::sin(2 * M_PI * x) * Kokkos::cos(2 * M_PI * y);
}

/// Checks that `values` are within 10% of SinxCosy at the mesh vertices.
void CheckAtVertices(Omega_h::Mesh& mesh,
                     const Kokkos::View<Real**, HostMemorySpace>& values)
{
  const auto coords = pcms::test::CopyOmegaHRealsToVector(mesh.coords());
  REQUIRE(values.extent(0) == static_cast<size_t>(mesh.nverts()));
  for (size_t i = 0; i < static_cast<size_t>(mesh.nverts()); ++i) {
    const Real exact = SinxCosy(coords[2 * i], coords[2 * i + 1]);
    CAPTURE(i);
    REQUIRE_THAT(values(i, 0), WithinAbs(exact, 0.1 * std::abs(exact)));
  }
}

} // namespace

TEST_CASE("MLS reconstruction from a mesh to vertices")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto source_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 20, 20, 0, false);
  const auto system = pcms::csys::Cartesian::Create(2);
  const pcms::NeighbourhoodOptions supports{.radius = 0.12, .min_supports = 15};
  const auto source_geometry = pcms::test::MakeOmegaHGeometry(source_mesh);
  auto vertices = std::make_shared<pcms::OmegaHEntityLayout>(source_geometry,
                                                             pcms::Vertex, 1);

  SECTION("Single Mesh")
  {
    // element centroids to the vertices of the same mesh, by the mesh
    // adjacency and by the distance search over the centroids as a point cloud
    auto mesh_space = pcms::CreatePolynomialReconstructionSpace(
      pcms::point_cloud::CreateGeometry(
        std::make_shared<pcms::OmegaHEntityLayout>(source_geometry, pcms::Face,
                                                   1),
        supports));
    auto centroids = pcms::test::CopyOmegaHRealsToVector(
      pcms::get_entity_centroids(source_mesh, Omega_h::FACE));
    auto cloud_space = pcms::CreatePolynomialReconstructionSpace(
      pcms::point_cloud::CreateGeometry(
        pcms::CoordinateView<HostMemorySpace>(
          system, pcms::Rank2View<const Real, HostMemorySpace>(
                    centroids.data(), source_mesh.nfaces(), 2)),
        supports),
      {});
    const auto mesh_bound = mesh_space->GetGeometry()->Bind(
      vertices->GetDOFHolderCoordinates(), pcms::DOFHolderOwners(*vertices));
    const auto cloud_bound =
      cloud_space->GetGeometry()->Bind(vertices->GetDOFHolderCoordinates());

    auto mesh_field = mesh_space->CreateFunction<Real>();
    pcms::test::SetField(
      mesh_field, OMEGA_H_LAMBDA(Real x, Real y) { return SinxCosy(x, y); });
    auto cloud_field = cloud_space->CreateFunction<Real>();
    pcms::test::SetField(
      cloud_field, OMEGA_H_LAMBDA(Real x, Real y) { return SinxCosy(x, y); });
    const auto mesh_values = pcms::test::EvaluateToHost(
      *mesh_space->CreatePointEvaluator<Real>(mesh_bound), mesh_field,
      mesh_field.GetData().GetValueBasis());
    const auto cloud_values = pcms::test::EvaluateToHost(
      *cloud_space->CreatePointEvaluator<Real>(cloud_bound), cloud_field,
      cloud_field.GetData().GetValueBasis());

    CheckAtVertices(source_mesh, mesh_values);
    CheckAtVertices(source_mesh, cloud_values);
    // both searches find the same supports, so the values agree
    pcms::test::CheckSupportResultsEquivalent(
      dynamic_cast<const pcms::MLSBinding&>(*cloud_bound).Supports(),
      dynamic_cast<const pcms::MLSBinding&>(*mesh_bound).Supports());
    for (int i = 0; i < source_mesh.nverts(); ++i) {
      REQUIRE_THAT(cloud_values(i, 0), WithinAbs(mesh_values(i, 0), 1e-4));
    }
  }

  SECTION("Double Mesh")
  {
    // the vertices of one mesh to the vertices of another
    auto target_mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 17, 17, 0, false);
    auto space = pcms::CreatePolynomialReconstructionSpace(
      pcms::point_cloud::CreateGeometry(vertices, supports));
    auto field = space->CreateFunction<Real>();
    pcms::test::SetField(
      field, OMEGA_H_LAMBDA(Real x, Real y) { return SinxCosy(x, y); });
    const pcms::OmegaHEntityLayout targets(
      pcms::test::MakeOmegaHGeometry(target_mesh), pcms::Vertex, 1);
    const auto values = pcms::test::EvaluateToHost(
      *space->CreatePointEvaluator<Real>(
        space->GetGeometry()->Bind(targets.GetDOFHolderCoordinates())),
      field, field.GetData().GetValueBasis());
    CheckAtVertices(target_mesh, values);
  }
}

TEST_CASE("MLS reconstruction honours the dimension of the points")
{
  std::vector<Real> points;
  for (int z = 0; z < 3; ++z) {
    for (int y = 0; y < 3; ++y) {
      for (int x = 0; x < 3; ++x) {
        points.insert(points.end(), {Real(x), Real(y), Real(z)});
      }
    }
  }
  const auto system = pcms::csys::Cartesian::Create(3);
  auto space = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(
      pcms::CoordinateView<HostMemorySpace>(
        system,
        pcms::Rank2View<const Real, HostMemorySpace>(points.data(), 27, 3)),
      pcms::NeighbourhoodOptions{.radius = 2.5, .min_supports = 10}),
    pcms::MLSOptions{.degree = 1});
  // a degree-1 polynomial that depends on z to catch accidental 2D behavior
  auto field = space->CreateFunction<Real>();
  pcms::test::SetField(
    field, OMEGA_H_LAMBDA(Real x, Real y, Real z) {
      return 1.0 + 2.0 * x - 3.0 * y + 5.0 * z;
    });
  const auto values = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(
      pcms::test::Bind(*space->GetGeometry(), points)),
    field, field.GetData().GetValueBasis());
  for (size_t i = 0; i < 27; ++i) {
    const Real exact = 1.0 + 2.0 * points[3 * i] - 3.0 * points[3 * i + 1] +
                       5.0 * points[3 * i + 2];
    CHECK_THAT(values(i, 0), WithinAbs(exact, 1e-6));
  }
}
