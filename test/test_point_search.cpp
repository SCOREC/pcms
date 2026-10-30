#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/localization/distributed_point_search.h>
#include <pcms/localization/point_search.h>
#include <pcms/utility/mesh_geometry.h>
#include <pcms/utility/uniform_grid.h>
#include "field_test_utils.h"
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <vector>

using pcms::AABBox;
using pcms::Uniform2DGrid;

TEST_CASE("global to local")
{
  Omega_h::Matrix<2, 3> coords{{0.0, 0.0}, {1, 0}, {0.5, 1}};
  SECTION("check verts")
  {
    for (int i = 0; i < coords.size(); ++i) {
      auto xi = Omega_h::barycentric_from_global<2, 2>(coords[i], coords);
      Omega_h::Vector<3> hand_xi{0, 0, 0};
      hand_xi[i] = 1;
      printf("[%f,%f,%f] == [%f,%f,%f]\n", xi[0], xi[1], xi[2], hand_xi[0],
             hand_xi[1], hand_xi[2]);
      REQUIRE(Omega_h::are_close(xi, hand_xi));
    }
  }
  SECTION("check point")
  {
    Omega_h::Vector<2> point{0.5, 0.5};
    auto xi = Omega_h::barycentric_from_global<2, 2>(point, coords);
    Omega_h::Vector<3> hand_xi{0.25, 0.25, 0.5};
    printf("[%f,%f,%f] == [%f,%f,%f]\n", xi[0], xi[1], xi[2], hand_xi[0],
           hand_xi[1], hand_xi[2]);
    REQUIRE(are_close(xi, hand_xi));
  }
}
TEST_CASE("Triangle AABB intersection")
{
  using pcms::triangle_intersects_bbox;
  AABBox<2> unit_square{.center = {0, 0}, .half_width = {0.5, 0.5}};
  SECTION("Triangle outside bbox")
  {
    REQUIRE(triangle_intersects_bbox({{1, 1}, {2, 1}, {2, 2}}, unit_square) ==
            false);
  }
  SECTION("BBOX inside Triangle")
  {
    REQUIRE(triangle_intersects_bbox({{-10, -10}, {10, -10}, {10, 10}},
                                     unit_square) == true);
  }
  SECTION("All verts inside BBOX")
  {
    REQUIRE(triangle_intersects_bbox({{-10, -10}, {10, -10}, {10, 10}},
                                     unit_square) == true);
  }
  SECTION("Single vert inside bbox")
  {
    REQUIRE(triangle_intersects_bbox({{0, 0}, {1, 0}, {1, 1}}, unit_square) ==
            true);
    REQUIRE(triangle_intersects_bbox({{1, 1}, {0, 0}, {1, 0}}, unit_square) ==
            true);
    REQUIRE(triangle_intersects_bbox({{1, 0}, {1, 1}, {0, 0}}, unit_square) ==
            true);
  }
  SECTION("edge crossed bbox")
  {
    // Triangle 1 with all vert permutations
    REQUIRE(triangle_intersects_bbox({{-0.6, 0}, {0, 0.6}, {-1, 1}},
                                     unit_square) == true);
    REQUIRE(triangle_intersects_bbox({{-1, 1}, {-0.6, 0}, {0, 0.6}},
                                     unit_square) == true);
    REQUIRE(triangle_intersects_bbox({{0, 0.6}, {-1, 1}, {-0.6, 0}},
                                     unit_square) == true);
    // triangle 2 with all vert permutations
    REQUIRE(triangle_intersects_bbox({{-0.6, 0}, {0.6, 0}, {0, -1}},
                                     unit_square) == true);
    REQUIRE(triangle_intersects_bbox({{0, -1}, {-0.6, 0}, {0.6, 0}},
                                     unit_square) == true);
    REQUIRE(triangle_intersects_bbox({{0.6, 0}, {0, -1}, {-0.6, 0}},
                                     unit_square) == true);
  }
}
TEST_CASE("Triangle BBox Intersection Regression")
{
  using pcms::triangle_intersects_bbox;
  AABBox<2> bbox{.center = {0.10833333, 0.10833333},
                 .half_width = {0.00833333, 0.00833333}};
  REQUIRE(triangle_intersects_bbox({{0, 0}, {0.1, 0}, {0.1, 0.1}}, bbox) ==
          true);
  REQUIRE(triangle_intersects_bbox({{0.1, 0}, {0.2, 0}, {0.2, 0.1}}, bbox) ==
          false);
  REQUIRE(triangle_intersects_bbox({{0.1, 0.1}, {0, 0.1}, {0, 0}}, bbox) ==
          true);
  REQUIRE(triangle_intersects_bbox({{0.1, 0.1}, {0.2, 0.1}, {0.2, 0.2}},
                                   bbox) == true);
  REQUIRE(triangle_intersects_bbox({{0, 0.1}, {0.1, 0.1}, {0.1, 0.2}}, bbox) ==
          true);
  REQUIRE(triangle_intersects_bbox({{0.1, 0.2}, {0, 0.2}, {0, 0.1}}, bbox) ==
          false);
  REQUIRE(triangle_intersects_bbox({{0.2, 0.2}, {0.1, 0.2}, {0.1, 0.1}},
                                   bbox) == true);
  REQUIRE(triangle_intersects_bbox({{0.2, 0.1}, {0.1, 0.1}, {0.1, 0}}, bbox) ==
          true);
}

template <typename T>
bool num_candidates_within_range(const T& intersection_map, pcms::LO min,
                                 pcms::LO max)
{
  using size_type = typename T::size_type;
  using MinMax = Kokkos::MinMax<size_type, Kokkos::HostSpace>;
  using result_type = typename MinMax::value_type;
  result_type result;
  Kokkos::parallel_reduce(
    intersection_map.numRows(),
    KOKKOS_LAMBDA(const int i, result_type& update) {
      auto num_candidates =
        intersection_map.row_map(i + 1) - intersection_map.row_map(i);
      if (num_candidates > update.max_val)
        update.max_val = num_candidates;
      if (num_candidates < update.min_val)
        update.min_val = num_candidates;
    },
    MinMax{result});
  bool within_range = (result.min_val >= min) && (result.max_val <= max);
  if (!within_range)
    std::cerr << result.min_val << ' ' << result.max_val << '\n';
  return within_range;
}
// extern Omega_h::Library omega_h_library;

TEST_CASE("construct intersection map")
{
  // auto world = omega_h_library.world();
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  REQUIRE(mesh.dim() == 2);
  SECTION("grid bbox overlap")
  {
    Kokkos::View<Uniform2DGrid[1]> grid_d("uniform grid");
    auto grid_h = Kokkos::create_mirror_view(grid_d);
    grid_h(0) = Uniform2DGrid{
      .edge_length{1, 1}, .bot_left = {0, 0}, .divisions = {10, 10}};
    Kokkos::deep_copy(grid_d, grid_h);
    auto intersection_map = pcms::detail::construct_intersection_map_2d(
      mesh, grid_d, grid_h(0).GetNumCells());
    // assert(cudaSuccess == cudaDeviceSynchronize());
    REQUIRE(intersection_map.numRows() == 100);
    REQUIRE(num_candidates_within_range(intersection_map, 2, 16));
  }
  SECTION("fine grid")
  {
    Kokkos::View<Uniform2DGrid[1]> grid_d("uniform grid");
    auto grid_h = Kokkos::create_mirror_view(grid_d);
    grid_h(0) = Uniform2DGrid{
      .edge_length{1, 1}, .bot_left = {0, 0}, .divisions = {60, 60}};
    Kokkos::deep_copy(grid_d, grid_h);
    // require number of candidates is >=1 and <=6
    auto intersection_map = pcms::detail::construct_intersection_map_2d(
      mesh, grid_d, grid_h(0).GetNumCells());

    REQUIRE(intersection_map.numRows() == 3600);
    REQUIRE(num_candidates_within_range(intersection_map, 1, 6));
  }
}
namespace
{
/// Assert two candidate maps are identical: same row offsets and the same
/// entries in the same order within every row.
template <typename T>
void require_maps_identical(const T& fast, const T& reference)
{
  REQUIRE(fast.numRows() == reference.numRows());
  REQUIRE(fast.entries.size() == reference.entries.size());

  auto fast_row_map =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, fast.row_map);
  auto ref_row_map =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, reference.row_map);
  auto fast_entries =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, fast.entries);
  auto ref_entries =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, reference.entries);

  for (std::size_t i = 0; i < fast_row_map.extent(0); ++i) {
    REQUIRE(fast_row_map(i) == ref_row_map(i));
  }
  for (std::size_t i = 0; i < fast_entries.extent(0); ++i) {
    REQUIRE(fast_entries(i) == ref_entries(i));
  }
}
} // namespace

// The element-major candidate map must produce exactly what the brute-force
// cell-major build does. Entry order matters as well as membership: the search
// returns the first candidate containing the point, so ordering decides which
// element wins a shared-face tie and which is reported for an out-of-bounds
// query.
TEST_CASE("candidate map matches brute-force reference")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  SECTION("2D")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);
    REQUIRE(mesh.dim() == 2);
    for (pcms::LO divisions : {3, 10, 37}) {
      Kokkos::View<Uniform2DGrid[1]> grid_d("uniform grid");
      auto grid_h = Kokkos::create_mirror_view(grid_d);
      grid_h(0) = Uniform2DGrid{.edge_length{1, 1},
                                .bot_left = {0, 0},
                                .divisions = {divisions, divisions}};
      Kokkos::deep_copy(grid_d, grid_h);
      const auto ncells = grid_h(0).GetNumCells();
      require_maps_identical(
        pcms::detail::construct_intersection_map_2d(mesh, grid_d, ncells),
        pcms::detail::construct_intersection_map_reference_2d(mesh, grid_d,
                                                              ncells));
    }
  }
  SECTION("3D")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 5, 5, 5, false);
    REQUIRE(mesh.dim() == 3);
    for (pcms::LO divisions : {2, 5, 11}) {
      Kokkos::View<pcms::Uniform3DGrid[1]> grid_d("uniform grid");
      auto grid_h = Kokkos::create_mirror_view(grid_d);
      grid_h(0) =
        pcms::Uniform3DGrid{.edge_length{1, 1, 1},
                            .bot_left = {0, 0, 0},
                            .divisions = {divisions, divisions, divisions}};
      Kokkos::deep_copy(grid_d, grid_h);
      const auto ncells = grid_h(0).GetNumCells();
      require_maps_identical(
        pcms::detail::construct_intersection_map_3d(mesh, grid_d, ncells),
        pcms::detail::construct_intersection_map_reference_3d(mesh, grid_d,
                                                              ncells));
    }
  }
}

TEST_CASE("uniform grid search")
{
  using pcms::GridPointSearch2D;
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  auto tolerances =
    GridPointSearch2D::PointSearchTolerances{"point search 2d tolerances"};
  auto tolerances_h = Kokkos::create_mirror_view(tolerances);
  tolerances_h(0) = 0.01;
  tolerances_h(1) = 0.01;
  Kokkos::deep_copy(tolerances, tolerances_h);

  GridPointSearch2D search{mesh, 10, 10, tolerances};

  Kokkos::View<pcms::Real* [2]> points("test_points", 8);
  // Kokkos::View<pcms::Real*[2]> points("test_points", 1);
  auto points_h = Kokkos::create_mirror_view(points);
  points_h(0, 0) = 0;
  points_h(0, 1) = 0;
  points_h(1, 0) = 0.55;
  points_h(1, 1) = 0.54;
  points_h(2, 0) = 100;
  points_h(2, 1) = 100;
  points_h(3, 0) = 1;
  points_h(3, 1) = 1;
  points_h(4, 0) = -1;
  points_h(4, 1) = -1;
  points_h(5, 0) = 1.01;
  points_h(5, 1) = 0.95;
  points_h(6, 0) = 0.05;
  points_h(6, 1) = -0.01;
  points_h(7, 0) = 0.05;
  points_h(7, 1) = 0.02;
  Kokkos::deep_copy(points, points_h);
  auto results = search(points);
  auto results_h = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(results_h, results);
  SECTION("global coordinate within mesh")
  {
    {
      auto [dim, idx, coords] = results_h(0);
      const auto face_idx = search.GetOwningElementId(results_h(0));

      CAPTURE(idx);

      REQUIRE(dim == GridPointSearch2D::Result::Dimensionality::VERTEX);
      REQUIRE(idx == 0);
      REQUIRE(face_idx >= 0);
      REQUIRE(coords[0] == Catch::Approx(1));
      REQUIRE(coords[1] == Catch::Approx(0));
      REQUIRE(coords[2] == Catch::Approx(0));
    }
    {
      auto [dim, idx, coords] = results_h(1);
      const auto face_idx = search.GetOwningElementId(results_h(1));
      REQUIRE(dim == GridPointSearch2D::Result::Dimensionality::EDGE);
      REQUIRE(idx == 156);
      REQUIRE(face_idx >= 0);
      REQUIRE(coords[0] == Catch::Approx(0.5));
      REQUIRE(coords[1] == Catch::Approx(0.1));
      REQUIRE(coords[2] == Catch::Approx(0.4));
    }
    {
      auto [dim, idx, coords] = results_h(7);
      const auto face_idx = search.GetOwningElementId(results_h(7));
      REQUIRE(dim == GridPointSearch2D::Result::Dimensionality::FACE);
      REQUIRE(idx == 0);
      REQUIRE(face_idx >= 0);
    }
  }
  SECTION("Global coordinate outside mesh")
  {
    auto out_of_bounds = results_h(2);
    auto top_right = results_h(3);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::VERTEX);
    REQUIRE(out_of_bounds.IsOutside());
    REQUIRE(out_of_bounds.EntityId() == top_right.element_id);
    REQUIRE(search.GetOwningElementId(out_of_bounds) >= 0);

    out_of_bounds = results_h(4);
    auto bot_left = results_h(0);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::VERTEX);
    REQUIRE(out_of_bounds.IsOutside());
    REQUIRE(out_of_bounds.EntityId() == bot_left.element_id);
    REQUIRE(search.GetOwningElementId(out_of_bounds) >= 0);

    out_of_bounds = results_h(5);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::EDGE);
    REQUIRE(out_of_bounds.IsOutside());
    REQUIRE(out_of_bounds.EntityId() == 219);
    REQUIRE(search.GetOwningElementId(out_of_bounds) >= 0);

    // (0.05, -0.01) is within the edge tolerance of the bottom edge, so it is
    // located on that edge rather than reported outside.
    auto on_edge = results_h(6);
    REQUIRE(on_edge.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::EDGE);
    REQUIRE_FALSE(on_edge.IsOutside());
    REQUIRE(on_edge.element_id == 0);
    REQUIRE(search.GetOwningElementId(on_edge) >= 0);
  }
  SECTION("point on extension of an edge")
  {
    Kokkos::View<pcms::Real* [2]> ext_points("ext_test_points", 1);
    auto ext_h = Kokkos::create_mirror_view(ext_points);
    ext_h(0, 0) = 1.5;
    ext_h(0, 1) = 0.0;
    Kokkos::deep_copy(ext_points, ext_h);
    auto ext_results = search(ext_points);
    auto ext_results_h = Kokkos::create_mirror_view(ext_results);
    Kokkos::deep_copy(ext_results_h, ext_results);

    auto res = ext_results_h(0);
    // Must be out-of-bounds (negative element id)
    REQUIRE(res.element_id < 0);
    // Dimensionality must be VERTEX — the nearest entity is the
    // rightmost bottom vertex (1.0, 0).
    REQUIRE(res.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::VERTEX);
    // Owning element should still resolve to a valid face
    REQUIRE(search.GetOwningElementId(res) >= 0);
  }
}

TEST_CASE("point search locates points on shared faces",
          "[point_search][3d][regression]")
{
  using pcms::GridPointSearch3D;
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  const int source_divisions = 10;
  const int query_divisions = 5;
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, source_divisions,
                       source_divisions, source_divisions, false);
  auto query_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, query_divisions,
                       query_divisions, query_divisions, false);

  // Same grid resolution the intersection search uses for this mesh.
  const auto n = pcms::DivisionsPerAxisForMesh<3>(mesh.nelems());
  GridPointSearch3D search{mesh, n, n, n};

  const auto centroids_h =
    Omega_h::HostRead<Omega_h::Real>(pcms::get_entity_centroids(query_mesh, 3));
  const auto npoints = query_mesh.nelems();

  Kokkos::View<pcms::Real* [3]> points("query points", npoints);
  auto points_h = Kokkos::create_mirror_view(points);
  for (Omega_h::LO i = 0; i < npoints; ++i) {
    for (int d = 0; d < 3; ++d) {
      points_h(i, d) = centroids_h[i * 3 + d];
    }
  }
  Kokkos::deep_copy(points, points_h);

  auto results = search(points);
  auto results_h = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(results_h, results);

  std::vector<Omega_h::LO> not_located;
  for (Omega_h::LO i = 0; i < npoints; ++i) {
    if (results_h(i).element_id < 0) {
      not_located.push_back(i);
    }
  }
  CAPTURE(not_located.size(), npoints);
  REQUIRE(not_located.empty());
}

TEST_CASE("point search encodes outside points unambiguously",
          "[point_search][regression]")
{
  auto lib = Omega_h::Library{};

  SECTION("3D: nearest tetrahedron 0")
  {
    using pcms::GridPointSearch3D;
    using Dimensionality = GridPointSearch3D::Result::Dimensionality;
    auto mesh = pcms::test::BuildReferenceTet(lib);
    const auto n = pcms::DivisionsPerAxisForMesh<3>(mesh.nelems());
    GridPointSearch3D search{mesh, n, n, n};

    Kokkos::View<pcms::Real* [3]> points("points", 2);
    auto points_h = Kokkos::create_mirror_view(points);
    // inside the tet
    points_h(0, 0) = 0.1;
    points_h(0, 1) = 0.1;
    points_h(0, 2) = 0.1;
    // in the tet's bounding box but outside it: the nearest candidate is tet 0
    points_h(1, 0) = 0.5;
    points_h(1, 1) = 0.5;
    points_h(1, 2) = 0.5;
    Kokkos::deep_copy(points, points_h);

    auto results = search(points);
    auto results_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);
    auto owners_h = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, search.GetOwningElementIds(results));

    REQUIRE(results_h(0).dimensionality == Dimensionality::REGION);
    REQUIRE_FALSE(results_h(0).IsOutside());
    REQUIRE(results_h(0).element_id == 0);
    REQUIRE(results_h(0).EntityId() == 0);
    REQUIRE(owners_h(0) == 0);

    REQUIRE(results_h(1).dimensionality == Dimensionality::REGION);
    REQUIRE(results_h(1).IsOutside());
    REQUIRE(results_h(1).EntityId() == 0);
    REQUIRE(search.GetOwningElementId(results_h(1)) == 0);
    REQUIRE(owners_h(1) == 0);

    // The grid over a single tet has no empty cell, so decode the
    // no-candidate result directly.
    Kokkos::View<GridPointSearch3D::Result*> none("none", 1);
    Kokkos::deep_copy(none,
                      GridPointSearch3D::Result{Dimensionality::NO_INTERSECT,
                                                -1, Omega_h::zero_vector<4>()});
    auto none_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, none);
    auto none_owners_h = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, search.GetOwningElementIds(none));
    REQUIRE(none_h(0).IsOutside());
    REQUIRE(none_h(0).EntityId() == -1);
    REQUIRE(search.GetOwningElementId(none_h(0)) == -1);
    REQUIRE(none_owners_h(0) == -1);
  }

  SECTION("2D: nearest vertex 0 and an empty query cell")
  {
    using pcms::GridPointSearch2D;
    using Dimensionality = GridPointSearch2D::Result::Dimensionality;
    // Only the triangle below the diagonal, so the cells above it hold no
    // candidates.
    auto mesh = pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3}));
    const auto n = pcms::DivisionsPerAxisForMesh<2>(mesh.nelems());
    GridPointSearch2D search{mesh, n, n};

    Kokkos::View<pcms::Real* [2]> points("points", 3);
    auto points_h = Kokkos::create_mirror_view(points);
    // beyond vertex 0: the nearest entity is vertex 0
    points_h(0, 0) = -0.05;
    points_h(0, 1) = -0.05;
    // in a grid cell the triangle does not reach
    points_h(1, 0) = 0.6;
    points_h(1, 1) = 0.6;
    // inside the triangle
    points_h(2, 0) = 0.2;
    points_h(2, 1) = 0.2;
    Kokkos::deep_copy(points, points_h);

    auto results = search(points);
    auto results_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);
    auto owners_h = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, search.GetOwningElementIds(results));

    REQUIRE(results_h(0).dimensionality == Dimensionality::VERTEX);
    REQUIRE(results_h(0).IsOutside());
    REQUIRE(results_h(0).EntityId() == 0);
    REQUIRE(search.GetOwningElementId(results_h(0)) == 0);
    REQUIRE(owners_h(0) == 0);

    REQUIRE(results_h(1).dimensionality == Dimensionality::NO_INTERSECT);
    REQUIRE(results_h(1).IsOutside());
    REQUIRE(results_h(1).EntityId() == -1);
    REQUIRE(search.GetOwningElementId(results_h(1)) == -1);
    REQUIRE(owners_h(1) == -1);

    REQUIRE(results_h(2).dimensionality == Dimensionality::FACE);
    REQUIRE_FALSE(results_h(2).IsOutside());
    REQUIRE(results_h(2).element_id == 0);
    REQUIRE(owners_h(2) == 0);
  }
}

// Baseline for the reused local search: on a mesh with no ghost layer every
// element is owned and every element centroid must localize to its own element.
// The distributed search treats a local answer as authoritative only when this
// holds, so it is the contract the distributed layer builds on.
TEST_CASE("ghost-free search localizes every element centroid")
{
  using pcms::GridPointSearch2D;
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 4, 4, 0, false);
  REQUIRE(mesh.nghost_layers() == 0);

  const auto n_elems = mesh.nelems();
  const auto owned = Omega_h::HostRead<Omega_h::I8>(mesh.owned(2));
  REQUIRE(static_cast<Omega_h::LO>(owned.size()) == n_elems);
  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    REQUIRE(owned[i] != 0);
  }

  const auto n = pcms::DivisionsPerAxisForGlobalMesh<2>(mesh);
  GridPointSearch2D search{mesh, n, n};
  const auto centroids =
    Omega_h::HostRead<Omega_h::Real>(pcms::get_entity_centroids(mesh, 2));
  Kokkos::View<pcms::Real* [2]> points("query points", n_elems);
  auto points_h = Kokkos::create_mirror_view(points);
  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    points_h(i, 0) = centroids[i * 2 + 0];
    points_h(i, 1) = centroids[i * 2 + 1];
  }
  Kokkos::deep_copy(points, points_h);

  const auto owners_h = Kokkos::create_mirror_view_and_copy(
    Kokkos::HostSpace{}, search.GetOwningElementIds(search(points)));
  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    REQUIRE(owners_h(i) == i);
  }
}

// The distributed search decides what this rank may answer for by asking
// whether the element the local search reported is one it owns. On a mesh with
// no ghost layer every local element is owned, so every located point must be
// answerable here -- which is the baseline the multi-rank classification test
// builds on.
TEST_CASE("distributed search classifies local answers by ownership")
{
  using pcms::DistributedPointSearch;
  using Status = DistributedPointSearch<2>::LocalStatus;

  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 4, 4, 0, false);
  REQUIRE(mesh.nghost_layers() == 0);

  const auto n = pcms::DivisionsPerAxisForGlobalMesh<2>(mesh);
  pcms::GridPointSearch2D local_search{mesh, n, n};
  const DistributedPointSearch<2> search{mesh, local_search};

  REQUIRE_FALSE(search.HasGhostLayer());
  REQUIRE(search.GetMesh().nelems() == mesh.nelems());

  // Query the centroid of every element plus one point outside the domain.
  const auto n_elems = mesh.nelems();
  const auto centroids =
    Omega_h::HostRead<Omega_h::Real>(pcms::get_entity_centroids(mesh, 2));
  Kokkos::View<pcms::Real* [2]> points("query points", n_elems + 1);
  auto points_h = Kokkos::create_mirror_view(points);
  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    points_h(i, 0) = centroids[i * 2 + 0];
    points_h(i, 1) = centroids[i * 2 + 1];
  }
  // Outside the mesh and outside every grid cell that holds an element.
  points_h(n_elems, 0) = 2.0;
  points_h(n_elems, 1) = 2.0;
  Kokkos::deep_copy(points, points_h);

  const auto classified = search.ClassifyLocally(points);
  REQUIRE(classified.size() == static_cast<size_t>(n_elems + 1));

  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    const auto& c = classified[static_cast<size_t>(i)];
    REQUIRE(c.status == Status::LOCATED_HERE);
    // An element centroid lies strictly inside its element, so the reported
    // element is that element and the barycentric coordinates are positive.
    REQUIRE(c.element_id == i);
    for (int d = 0; d < 3; ++d) {
      REQUIRE(c.barycentric[d] > 0);
    }
  }

  // With no ghosts nothing can be "located elsewhere": either this rank answers
  // or the point is outside what it holds.
  for (const auto& c : classified) {
    REQUIRE(c.status != Status::LOCATED_ELSEWHERE);
  }
  REQUIRE(classified[static_cast<size_t>(n_elems)].status ==
          Status::NOT_FOUND_HERE);
  REQUIRE(classified[static_cast<size_t>(n_elems)].element_id == -1);
}
// A predicate chosen at the map-building call must govern *both* halves: the
// range of cells visited and the test applied to each. Checking the map against
// a brute-force build made with that same predicate is what catches the case
// where one half ignores the argument -- the fast map would then disagree.
//
// Only predicates instantiated in point_search.cpp can be named from here; the
// set is small and internal, which is why the map machinery no longer lives in
// the header. See the instantiations at the end of that file.
TEST_CASE("candidate map accepts a predicate other than the default")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);

  Kokkos::View<Uniform2DGrid[1]> grid_d("uniform grid");
  auto grid_h = Kokkos::create_mirror_view(grid_d);
  grid_h(0) = Uniform2DGrid{
    .edge_length{1, 1}, .bot_left = {0, 0}, .divisions = {10, 10}};
  Kokkos::deep_copy(grid_d, grid_h);
  const auto ncells = grid_h(0).GetNumCells();

  // SimplexVsBbox is the 3D default; in 2D it is the conservative bounding-box
  // test, which is a legitimate thing to want. No other caller asks for it in
  // 2D.
  auto conservative =
    pcms::detail::construct_intersection_map<2, pcms::SimplexVsBbox>(
      mesh, grid_d, ncells);
  require_maps_identical(
    conservative,
    pcms::detail::construct_intersection_map_reference<2, pcms::SimplexVsBbox>(
      mesh, grid_d, ncells));

  // The conservative test accepts everything the refined default accepts and
  // more, wherever an element's box reaches a cell its triangle does not. If
  // the predicate were ignored, the two maps would be identical.
  auto refined =
    pcms::detail::construct_intersection_map_2d(mesh, grid_d, ncells);
  REQUIRE(conservative.entries.size() > refined.entries.size());
}

// The search itself must accept a predicate the same way, chosen once when it
// is built. The predicate must be one instantiated in point_search.cpp.
TEST_CASE("point search accepts a predicate other than its default")
{
  using pcms::GridPointSearch2D;
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);

  auto tolerances =
    GridPointSearch2D::PointSearchTolerances{"point search 2d tolerances"};
  auto tolerances_h = Kokkos::create_mirror_view(tolerances);
  tolerances_h(0) = 0.01;
  tolerances_h(1) = 0.01;
  Kokkos::deep_copy(tolerances, tolerances_h);

  GridPointSearch2D default_search{mesh, 10, 10, tolerances};
  // The conservative 2D predicate instead of this search's default.
  GridPointSearch2D custom_search{mesh, 10, 10, tolerances,
                                  pcms::SimplexVsBbox{}};

  Kokkos::View<pcms::Real* [2]> points("test_points", 4);
  auto points_h = Kokkos::create_mirror_view(points);
  points_h(0, 0) = 0.05;
  points_h(0, 1) = 0.02;
  points_h(1, 0) = 0.55;
  points_h(1, 1) = 0.54;
  points_h(2, 0) = 0.9;
  points_h(2, 1) = 0.9;
  points_h(3, 0) = 1.0;
  points_h(3, 1) = 1.0;
  Kokkos::deep_copy(points, points_h);

  auto default_results = default_search(points);
  auto custom_results = custom_search(points);
  auto default_h = Kokkos::create_mirror_view(default_results);
  auto custom_h = Kokkos::create_mirror_view(custom_results);
  Kokkos::deep_copy(default_h, default_results);
  Kokkos::deep_copy(custom_h, custom_results);

  // The query tests point containment in each candidate, so the predicate only
  // decides which elements are offered. For points inside the mesh both
  // searches must therefore answer identically: the value-constructor path is a
  // candidate-set choice, not a change of semantics.
  for (size_t i = 0; i < 4; ++i) {
    CAPTURE(i);
    REQUIRE(custom_h(i).dimensionality == default_h(i).dimensionality);
    REQUIRE(custom_h(i).element_id == default_h(i).element_id);
  }
}

// The 3D search takes a division count per axis. Unequal counts must locate
// points exactly as equal counts do: elements are binned by each cell's box, so
// a cell box computed for the wrong axis would drop candidates silently.
TEST_CASE("3D point search with unequal divisions per axis",
          "[point_search][3d]")
{
  using pcms::GridPointSearch3D;
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 4, 4, 4, false);
  GridPointSearch3D search{mesh, 3, 5, 7};

  const auto n_elems = mesh.nelems();
  const auto centroids =
    Omega_h::HostRead<Omega_h::Real>(pcms::get_entity_centroids(mesh, 3));
  Kokkos::View<pcms::Real* [3]> points("query points", n_elems);
  auto points_h = Kokkos::create_mirror_view(points);
  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    for (int d = 0; d < 3; ++d) {
      points_h(i, d) = centroids[i * 3 + d];
    }
  }
  Kokkos::deep_copy(points, points_h);

  const auto owners_h = Kokkos::create_mirror_view_and_copy(
    Kokkos::HostSpace{}, search.GetOwningElementIds(search(points)));
  Omega_h::LO n_wrong = 0;
  for (Omega_h::LO i = 0; i < n_elems; ++i) {
    n_wrong += owners_h(i) != i;
  }
  CAPTURE(n_wrong, n_elems);
  REQUIRE(n_wrong == 0);
}

// The grid and overlap layer under the point search is written for any
// dimension, though the search classes themselves are 2D and 3D only. These
// checks pin that down above 3D. Every grid has a different division count on
// each axis, and every probe uses different indices per axis, so an axis mix-up
// cannot hide behind a symmetric case.
namespace
{

/// Cells of `grid` whose id, per-axis indices and box disagree, counted on the
/// device: the same functions the search runs in its kernels.
template <unsigned Dim>
pcms::LO CountCellRoundTripFailuresOnDevice(const pcms::UniformGrid<Dim>& grid)
{
  pcms::LO failures = 0;
  Kokkos::parallel_reduce(
    "cell round trip", grid.GetNumCells(),
    KOKKOS_LAMBDA(pcms::LO cell, pcms::LO & n) {
      const auto axis = grid.AxisIndicesFromCellIndex(cell);
      bool ok = grid.CellIndexFromAxisIndices(axis) == cell;
      Omega_h::Vector<Dim> center;
      const auto box = grid.GetCellBBOX(cell);
      for (unsigned d = 0; d < Dim; ++d) {
        ok = ok && axis[d] >= 0 && axis[d] < grid.divisions[d];
        center[d] = box.center[d];
      }
      ok = ok && grid.ClosestCellID(center) == cell;
      n += ok ? 0 : 1;
    },
    failures);
  return failures;
}

} // namespace

TEST_CASE("grid and overlap layer above 3D", "[point_search][nd]")
{
  using pcms::LO;

  SECTION("4D grid: cell ids, per-axis indices and cell boxes agree")
  {
    // Cell size 1 on every axis, so a cell's box is easy to state.
    const pcms::UniformGrid<4> grid{.edge_length = {2.0, 3.0, 4.0, 5.0},
                                    .bot_left = {0.0, 0.0, 0.0, 0.0},
                                    .divisions = {2, 3, 4, 5}};
    REQUIRE(grid.GetNumCells() == 120);

    // Ids number cells x-fastest: x + Nx * (y + Ny * (z + Nz * w)).
    const Kokkos::Array<LO, 4> last{1, 2, 3, 4};
    REQUIRE(grid.CellIndexFromAxisIndices(last) == 119);
    const Kokkos::Array<LO, 4> probe{1, 0, 2, 3};
    const LO probe_id = 1 + 2 * (0 + 3 * (2 + 4 * 3));
    REQUIRE(grid.CellIndexFromAxisIndices(probe) == probe_id);
    const auto back = grid.AxisIndicesFromCellIndex(probe_id);
    for (int d = 0; d < 4; ++d) {
      REQUIRE(back[d] == probe[d]);
    }
    const auto box = grid.GetCellBBOX(probe_id);
    for (int d = 0; d < 4; ++d) {
      REQUIRE(box.center[d] == Catch::Approx(probe[d] + 0.5));
      REQUIRE(box.half_width[d] == Catch::Approx(0.5));
    }

    REQUIRE(CountCellRoundTripFailuresOnDevice(grid) == 0);
  }

  SECTION("5D grid: every cell round-trips on the device")
  {
    const pcms::UniformGrid<5> grid{.edge_length = {1.0, 2.0, 3.0, 2.0, 1.0},
                                    .bot_left = {0.0, -1.0, 2.0, 0.5, 0.0},
                                    .divisions = {2, 3, 4, 3, 2}};
    REQUIRE(grid.GetNumCells() == 144);
    REQUIRE(CountCellRoundTripFailuresOnDevice(grid) == 0);
  }

  SECTION("4D candidate cells keep each axis separate")
  {
    const pcms::UniformGrid<4> grid{.edge_length = {2.0, 3.0, 4.0, 5.0},
                                    .bot_left = {0.0, 0.0, 0.0, 0.0},
                                    .divisions = {2, 3, 4, 5}};
    // Inside cell (0, 1, 2, 3); padded by one cell, clamped to the grid.
    const AABBox<4> box{.center = {0.5, 1.5, 2.5, 3.5},
                        .half_width = {0.25, 0.25, 0.25, 0.25}};
    const auto range = pcms::candidate_cells<4, pcms::SimplexVsBbox>(box, grid);
    const Kokkos::Array<LO, 4> lo{0, 0, 1, 2};
    const Kokkos::Array<LO, 4> hi{1, 2, 3, 4};
    for (int d = 0; d < 4; ++d) {
      CAPTURE(d);
      REQUIRE(range.lo[d] == lo[d]);
      REQUIRE(range.hi[d] == hi[d]);
    }
    REQUIRE(range.Contains(Kokkos::Array<LO, 4>{0, 1, 2, 3}));
    REQUIRE_FALSE(range.Contains(Kokkos::Array<LO, 4>{0, 1, 2, 1}));
  }

  SECTION("4D simplex-vs-box overlap needs every axis to overlap")
  {
    const AABBox<4> box{.center = {1.0, 1.0, 1.0, 1.0},
                        .half_width = {0.5, 0.5, 0.5, 0.5}};
    // A simplex spanning [0.8, 1.2] on axes 0-2. On axis 3 it either reaches
    // into the box or stops short of it.
    const auto simplex = [](pcms::Real axis3_lo, pcms::Real axis3_hi) {
      Omega_h::Matrix<4, 5> coords{};
      for (int v = 0; v < 5; ++v) {
        for (int d = 0; d < 3; ++d) {
          coords(d, v) = (v % 2 == 0) ? 0.8 : 1.2;
        }
        coords(3, v) = (v == 0) ? axis3_lo : axis3_hi;
      }
      return coords;
    };
    REQUIRE(
      pcms::Overlap<4, pcms::SimplexVsBbox>::apply(simplex(1.4, 2.0), box));
    // Overlapping on axes 0-2 but entirely past the box on axis 3.
    REQUIRE_FALSE(
      pcms::Overlap<4, pcms::SimplexVsBbox>::apply(simplex(1.6, 2.0), box));
  }
}
