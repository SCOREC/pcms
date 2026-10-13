#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <pcms/transfer/mesh_intersection.hpp>
#include <pcms/utility/assert.h>

#include <algorithm>
#include <vector>

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
    require_intersections_identical(
      pcms::intersectTargets(source, target, true),
      pcms::intersectTargets(source, target, false));
  }
  SECTION("2D, nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 12, 12, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 4, 4, 0, false);
    require_intersections_identical(
      pcms::intersectTargets(source, target, true),
      pcms::intersectTargets(source, target, false));
  }
  SECTION("2D, non-nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 7, 11, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 5, 3, 0, false);
    require_intersections_identical(
      pcms::intersectTargets(source, target, true),
      pcms::intersectTargets(source, target, false));
  }
  SECTION("3D, identical meshes")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 4, 4, 4, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 4, 4, 4, false);
    require_intersections_identical(
      pcms::intersectTargets(source, target, true),
      pcms::intersectTargets(source, target, false));
  }
  SECTION("3D, nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 6, 6, 6, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 3, 3, 3, false);
    require_intersections_identical(
      pcms::intersectTargets(source, target, true),
      pcms::intersectTargets(source, target, false));
  }
  SECTION("3D, non-nested resolutions")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 5, 3, 4, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 4, 3, false);
    require_intersections_identical(
      pcms::intersectTargets(source, target, true),
      pcms::intersectTargets(source, target, false));
  }
}

// The conservative transfer hands the intersection the source space's own
// point search rather than building a second one. The set of source elements
// per target is a property of the two meshes: which search located the target
// centroids, and at what grid resolution, must not show in it. Row order can
// (see require_intersections_equivalent), so this compares rows as sets.
TEST_CASE("intersection map is independent of the supplied source search",
          "[intersection]")
{
  Omega_h::Library lib;
  auto world = lib.world();

  SECTION("2D, coarse caller-supplied grid")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 7, 11, 0, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 5, 3, 0, false);
    const pcms::GridPointSearchVariant search{
      std::in_place_type<pcms::GridPointSearch2D>, source, 3, 3};
    require_intersections_equivalent(
      pcms::intersectTargets(source, target, search),
      pcms::intersectTargets(source, target));
  }
  SECTION("3D, coarse caller-supplied grid")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 5, 3, 4, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 4, 3, false);
    const pcms::GridPointSearchVariant search{
      std::in_place_type<pcms::GridPointSearch3D>, source, 2, 2, 2};
    require_intersections_equivalent(
      pcms::intersectTargets(source, target, search),
      pcms::intersectTargets(source, target));
  }
  SECTION("a search of the wrong dimension is rejected")
  {
    auto source =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 2, 2, false);
    auto target =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 2, 2, 2, false);
    auto flat =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 2, 2, 0, false);
    const pcms::GridPointSearchVariant search{
      std::in_place_type<pcms::GridPointSearch2D>, flat, 2, 2};
    REQUIRE_THROWS_AS(pcms::intersectTargets(source, target, search),
                      pcms::pcms_error);
  }
}
