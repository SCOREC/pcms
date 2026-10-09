// Regression tests for target elements that the mesh intersection misses
// entirely.
//
// Conservatively projecting u == 1 from a source mesh onto a target mesh that
// covers the same domain must give u == 1 on every target element: the value
// on target element e is (1/|e|) * integral over e of the source field, so it
// is exactly the fraction of e that the intersection quadrature found.  A
// target element that comes back 0 was never matched to any source element.
//
// Found on a 5.7M-element tetrahedral mesh coarsened 3x with Omega_h's
// adapter: 2528 of 227381 target elements came back exactly 0, losing 0.66% of
// the volume and 1.77% of the transferred field.  The distribution was
// all-or-nothing -- 224853 elements at coverage 1, 2528 at exactly 0, nothing
// between -- and every failing element's centroid lay on a face of the source
// mesh to within 1e-15 in barycentric coordinates.
//
// Nested uniform boxes reproduce it with no mesh files: the 5^3 target's
// element centroids land on source grid planes when the source is 10^3,
// because 10 is a multiple of 5.
//
// The seed of the intersection search is the source element containing the
// target element's centroid, and it is added to the map unconditionally.  A
// correct seed therefore always yields positive coverage -- it contains an
// interior point of the target -- so zero coverage means the point search
// returned an element that does not contain the centroid, after which the
// breadth-first search over the source dual graph has nowhere to expand.
// See the companion "point search locates points on shared faces" case in
// test_point_search.cpp.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/field/function_space/lagrange.h>
#include <pcms/transfer/mass_matrix_type.hpp>
#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/utility/arrays.h>

#include <cstddef>
#include <vector>

#include "field_test_utils.h"

namespace
{

// Per-target-element fraction of its volume that the intersection quadrature
// accounted for, obtained by projecting u == 1 onto a P0 target space.
std::vector<pcms::Real> CoveragePerTargetElement(Omega_h::Mesh& source_mesh,
                                                 Omega_h::Mesh& target_mesh)
{
  auto source_space = pcms::test::MakeP0Space(source_mesh);
  auto target_space = pcms::test::MakeP0Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source, KOKKOS_LAMBDA(pcms::Real, pcms::Real, pcms::Real) { return 1.0; });

  // Lumped: for a P0 target the mass matrix is diagonal either way, so this
  // is an exact division by the element volume and adds no error of its own.
  pcms::OmegaHConservativeProjection projection(
    *source_space, *target_space, pcms::MassMatrixType::Lumped);
  projection.Apply(source, target);

  const auto values =
    pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
  std::vector<pcms::Real> coverage(values.size());
  for (std::size_t i = 0; i < coverage.size(); ++i) {
    coverage[i] = values[i];
  }
  return coverage;
}

} // namespace

TEST_CASE("mesh intersection covers every target element",
          "[transfer][mesh_intersection][3d][regression]")
{
  Omega_h::Library lib;

  // 10 -> 5 is the smallest nested pair that reproduces it. 12 -> 6, 12 -> 7,
  // 12 -> 9, 12 -> 10 and 12 -> 11 fail too, with 1 to 7 empty elements.
  const int source_divisions = 10;
  const int target_divisions = 5;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitCube(lib, source_divisions);
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitCube(lib, target_divisions);

  const auto coverage = CoveragePerTargetElement(source_mesh, target_mesh);
  REQUIRE(static_cast<Omega_h::LO>(coverage.size()) == target_mesh.nelems());

  const auto sizes_h = Omega_h::HostRead<Omega_h::Real>(
    Omega_h::measure_elements_real(&target_mesh));

  // Report every miss rather than just the first, so a change in how many
  // elements are affected is visible in the failure output.
  std::vector<Omega_h::LO> empty;
  pcms::Real covered_volume = 0.0;
  pcms::Real total_volume = 0.0;
  for (Omega_h::LO e = 0; e < target_mesh.nelems(); ++e) {
    covered_volume += coverage[e] * sizes_h[e];
    total_volume += sizes_h[e];
    if (coverage[e] <= 1e-12) {
      empty.push_back(e);
    }
  }
  CAPTURE(empty.size(), covered_volume, total_volume);
  CHECK(empty.empty());

  for (Omega_h::LO e = 0; e < target_mesh.nelems(); ++e) {
    CAPTURE(e, coverage[e]);
    REQUIRE(coverage[e] == Catch::Approx(1.0).margin(1e-9));
  }
  REQUIRE(covered_volume == Catch::Approx(total_volume).margin(1e-9));
}
