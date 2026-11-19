#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/backends/omega_h/geometry.hpp>
#include <pcms/backends/omega_h/mesh_intersection.hpp>
#include <pcms/utility/assert.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "field_test_utils.hpp"

namespace pcms
{
using omega_h::IntersectionResults;
using omega_h::intersectTargets;
using omega_h::OmegaHGeometry;
using omega_h::simplices_have_degenerate_overlap;
} // namespace pcms

namespace
{
/// Assert the two intersection maps are identical: same per-target counts and
/// the same source elements in the same order.
void require_intersections_identical(const pcms::IntersectionResults& filtered,
                                     const pcms::IntersectionResults& exact)
{
  const auto filtered_offsets = Omega_h::HostRead(filtered.tgt2src_offsets);
  const auto exact_offsets = Omega_h::HostRead(exact.tgt2src_offsets);
  REQUIRE(filtered_offsets.size() == exact_offsets.size());
  for (int i = 0; i < filtered_offsets.size(); ++i) {
    REQUIRE(filtered_offsets[i] == exact_offsets[i]);
  }

  const auto filtered_indices = Omega_h::HostRead(filtered.tgt2src_indices);
  const auto exact_indices = Omega_h::HostRead(exact.tgt2src_indices);
  REQUIRE(filtered_indices.size() == exact_indices.size());
  for (int i = 0; i < filtered_indices.size(); ++i) {
    REQUIRE(filtered_indices[i] == exact_indices[i]);
  }
}

/// Assert the two intersection maps hold the same source elements per target,
/// ignoring order. The order of a row is seed-then-BFS, and the seed is
/// whichever containing element the point search reports first; a centroid
/// lying on a shared face can legitimately be located in either neighbor.
void require_intersections_equivalent(const pcms::IntersectionResults& a,
                                      const pcms::IntersectionResults& b)
{
  const auto a_offsets = Omega_h::HostRead(a.tgt2src_offsets);
  const auto b_offsets = Omega_h::HostRead(b.tgt2src_offsets);
  REQUIRE(a_offsets.size() == b_offsets.size());
  const auto a_indices = Omega_h::HostRead(a.tgt2src_indices);
  const auto b_indices = Omega_h::HostRead(b.tgt2src_indices);
  for (int row = 0; row + 1 < a_offsets.size(); ++row) {
    REQUIRE(a_offsets[row] == b_offsets[row]);
    std::vector<int> a_row(a_indices.data() + a_offsets[row],
                           a_indices.data() + a_offsets[row + 1]);
    std::vector<int> b_row(b_indices.data() + b_offsets[row],
                           b_indices.data() + b_offsets[row + 1]);
    std::sort(a_row.begin(), a_row.end());
    std::sort(b_row.begin(), b_row.end());
    REQUIRE(a_row == b_row);
  }
}
} // namespace

namespace
{
template <int Dim>
void check_outside_centroid_overlap(Omega_h::Library& lib)
{
  const auto build_mesh = [&](int divisions) {
    if constexpr (Dim == 3) {
      return pcms::test::BuildReferenceTet(lib);
    } else {
      return Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0,
                                divisions, divisions, 0, false);
    }
  };
  auto source = build_mesh(2);
  auto target = build_mesh(1);
  // 2D: overlap is the slab [0.9, 1] x [0, 1], area 0.1.
  // 3D: two reference tets shifted along x overlap in a tet of volume 0.1^3/6.
  // Every target centroid is outside the source, but real overlaps remain.
  const auto original = target.coords();
  Omega_h::Write<Omega_h::Real> shifted(original.size());
  Omega_h::parallel_for(
    original.size(), OMEGA_H_LAMBDA(int i) {
      shifted[i] = original[i] + (i % Dim == 0 ? 0.9 : 0.0);
    });
  target.set_coords(shifted);
  const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
    source, pcms::csys::Cartesian::Create(Dim));
  const auto source_coords = Omega_h::HostRead(source.coords());
  const auto target_coords = Omega_h::HostRead(target.coords());
  const auto source_nodes = Omega_h::HostRead(source.ask_elem_verts());
  const auto target_nodes = Omega_h::HostRead(target.ask_elem_verts());

  for (bool prefilter : {false, true}) {
    CAPTURE(Dim, prefilter);
    const auto result = pcms::intersectTargets(*geometry, target, prefilter);
    const auto offsets = Omega_h::HostRead(result.tgt2src_offsets);
    const auto indices = Omega_h::HostRead(result.tgt2src_indices);
    double total_measure = 0;
    for (int t = 0; t < target.nelems(); ++t) {
      r3d::Few<r3d::Vector<Dim>, Dim + 1> tgt;
      double centroid_x = 0;
      for (int v = 0; v <= Dim; ++v) {
        for (int d = 0; d < Dim; ++d) {
          tgt[v][d] = target_coords[target_nodes[t * (Dim + 1) + v] * Dim + d];
        }
        centroid_x += tgt[v][0] / (Dim + 1);
      }
      REQUIRE(centroid_x > 1.0);
      // Exhaustive clipping is independent of centroid seeds and BFS.
      for (int s = 0; s < source.nelems(); ++s) {
        r3d::Few<r3d::Vector<Dim>, Dim + 1> src;
        for (int v = 0; v <= Dim; ++v) {
          for (int d = 0; d < Dim; ++d) {
            src[v][d] =
              source_coords[source_nodes[s * (Dim + 1) + v] * Dim + d];
          }
        }
        r3d::Polytope<Dim> intersection;
        r3d::intersect_simplices(intersection, tgt, src);
        const double measure = std::abs(r3d::measure(intersection));
        if (measure <= 1e-12) {
          continue;
        }
        const auto begin = indices.data() + offsets[t];
        const auto end = indices.data() + offsets[t + 1];
        CAPTURE(t, s, measure);
        REQUIRE(std::count(begin, end, s) == 1);
        total_measure += measure;
      }
    }
    const double expected = Dim == 3 ? 0.001 / 6 : 0.1;
    REQUIRE(total_measure == Catch::Approx(expected).margin(1e-12));
  }
}
} // namespace

TEST_CASE("intersection preserves overlaps with outside target centroids",
          "[intersection]")
{
  Omega_h::Library lib;
  SECTION("2D")
  {
    check_outside_centroid_overlap<2>(lib);
  }
  SECTION("3D")
  {
    check_outside_centroid_overlap<3>(lib);
  }
}

// A pre-filter that never fires would satisfy the equivalence test below
// while saving nothing, so pin down that it actually rejects the case it
// exists for, and that it does not reject a real overlap.
TEST_CASE("intersection prefilter rejects degenerate overlaps",
          "[intersection]")
{
  SECTION("3D, face-adjacent tets are rejected")
  {
    // Two tets sharing the face x=0; separated exactly by that plane.
    const r3d::Few<r3d::Vector<3>, 4> left{
      {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}}};
    const r3d::Few<r3d::Vector<3>, 4> right{
      {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 0}}};
    REQUIRE(pcms::simplices_have_degenerate_overlap<3>(left, right));
  }
  SECTION("3D, overlapping tets are not rejected")
  {
    const r3d::Few<r3d::Vector<3>, 4> a{
      {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    // Same tet shifted slightly; the overlap has positive volume.
    const r3d::Few<r3d::Vector<3>, 4> b{
      {{0.1, 0.1, 0.1}, {1.1, 0.1, 0.1}, {0.1, 1.1, 0.1}, {0.1, 0.1, 1.1}}};
    REQUIRE_FALSE(pcms::simplices_have_degenerate_overlap<3>(a, b));
  }
  SECTION("2D, edge-adjacent triangles are rejected")
  {
    const r3d::Few<r3d::Vector<2>, 3> left{{{0, 0}, {0, 1}, {-1, 0}}};
    const r3d::Few<r3d::Vector<2>, 3> right{{{0, 0}, {0, 1}, {1, 0}}};
    REQUIRE(pcms::simplices_have_degenerate_overlap<2>(left, right));
  }
  SECTION("2D, overlapping triangles are not rejected")
  {
    const r3d::Few<r3d::Vector<2>, 3> a{{{0, 0}, {1, 0}, {0, 1}}};
    const r3d::Few<r3d::Vector<2>, 3> b{{{0.1, 0.1}, {1.1, 0.1}, {0.1, 1.1}}};
    REQUIRE_FALSE(pcms::simplices_have_degenerate_overlap<2>(a, b));
  }
}

// The separating-plane pre-filter must not change which source elements are
// found, only how quickly the degenerate pairs are discarded. The unfiltered
// path is the oracle.
TEST_CASE("intersection prefilter preserves the intersection map",
          "[intersection]")
{
  Omega_h::Library lib;
  auto world = lib.world();

  SECTION("2D, identical meshes")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 8, 8, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 8, 8, 0, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(2));
    require_intersections_identical(
      pcms::intersectTargets(*geometry, target, true),
      pcms::intersectTargets(*geometry, target, false));
  }
  SECTION("2D, nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 12, 12, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 4, 4, 0, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(2));
    require_intersections_identical(
      pcms::intersectTargets(*geometry, target, true),
      pcms::intersectTargets(*geometry, target, false));
  }
  SECTION("2D, non-nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 7, 11, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 5, 3, 0, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(2));
    require_intersections_identical(
      pcms::intersectTargets(*geometry, target, true),
      pcms::intersectTargets(*geometry, target, false));
  }
  SECTION("3D, identical meshes")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 4, 4, 4, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 4, 4, 4, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(3));
    require_intersections_identical(
      pcms::intersectTargets(*geometry, target, true),
      pcms::intersectTargets(*geometry, target, false));
  }
  SECTION("3D, nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 6, 6, 6, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 3, 3, 3, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(3));
    require_intersections_identical(
      pcms::intersectTargets(*geometry, target, true),
      pcms::intersectTargets(*geometry, target, false));
  }
  SECTION("3D, non-nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 5, 3, 4, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 4, 3, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(3));
    require_intersections_identical(
      pcms::intersectTargets(*geometry, target, true),
      pcms::intersectTargets(*geometry, target, false));
  }
}

// The conservative transfer hands the intersection the source space's
// geometry, which locates the target centroids with its own point search. The
// set of source elements per target is a property of the two meshes: at what
// grid resolution the geometry searches must not show in it. Row order can
// (see require_intersections_equivalent), so this compares rows as sets.
TEST_CASE("intersection map is independent of the source search resolution",
          "[intersection]")
{
  Omega_h::Library lib;
  auto world = lib.world();

  SECTION("2D, coarse geometry search grid")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 7, 11, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 5, 3, 0, false);
    const auto coarse = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(2),
      pcms::omega_h::OmegaHGeometryOptions{.search = pcms::SearchOptions{3}});
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(2));
    require_intersections_equivalent(pcms::intersectTargets(*coarse, target),
                                     pcms::intersectTargets(*geometry, target));
  }
  SECTION("3D, coarse geometry search grid")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 5, 3, 4, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 4, 3, false);
    const auto coarse = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(3),
      pcms::omega_h::OmegaHGeometryOptions{.search = pcms::SearchOptions{2}});
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      source, pcms::csys::Cartesian::Create(3));
    require_intersections_equivalent(pcms::intersectTargets(*coarse, target),
                                     pcms::intersectTargets(*geometry, target));
  }
  SECTION("a 2D geometry cannot intersect a 3D target")
  {
    auto flat =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 2, 2, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 2, 2, false);
    const auto geometry = std::make_shared<pcms::OmegaHGeometry>(
      flat, pcms::csys::Cartesian::Create(2));
    REQUIRE_THROWS_AS(pcms::intersectTargets(*geometry, target),
                      pcms::pcms_error);
  }
}
