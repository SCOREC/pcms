#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/backends/omega_h/point_search.hpp>
#include <pcms/localization/uniform_grid_search.hpp>
#include <pcms/utility/mesh_geometry.hpp>
#include <pcms/utility/uniform_grid.hpp>
#include "field_test_utils.hpp"
#include <pcms/backends/omega_h/geometry.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace pcms
{
using backend::omega_h::GridPointSearch2D;
using backend::omega_h::GridPointSearch3D;
using backend::omega_h::triangle_intersects_bbox;
using localization::UniformGridPointSearch;
using utility::get_entity_centroids;
} // namespace pcms

using pcms::AABBox;
using pcms::Uniform2DGrid;

namespace
{

/// Lifts each result to the lowest element adjacent to its entity, with the
/// result's weights on the entity's vertices and zero on the others, as an
/// evaluator does; a functor, not a lambda, so that it compiles under nvcc in
/// this anonymous namespace.
template <int Dim, typename Result>
struct LiftResults
{
  Omega_h::Few<Omega_h::Adj, Dim> up;
  Omega_h::Few<Omega_h::LOs, Dim + 1> verts_of;
  Kokkos::View<const Result*> results;
  Kokkos::View<Omega_h::LO*> elems;
  Kokkos::View<pcms::Real**, Kokkos::LayoutRight> bary;

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    const auto elem =
      pcms::backend::omega_h::detail::ElementOfResult<Dim>(up, results(i));
    elems(i) = elem;
    Omega_h::LO vertices[Dim + 1];
    pcms::Real weights[Dim + 1];
    const int nv = pcms::backend::omega_h::detail::EntityVertexWeights<Dim>(
      verts_of, results(i), vertices, weights);
    for (int k = 0; k <= Dim; ++k) {
      bary(i, k) = 0.0;
    }
    for (int j = 0; elem >= 0 && j < nv; ++j) {
      for (int k = 0; k <= Dim; ++k) {
        if (verts_of[Dim][elem * (Dim + 1) + k] == vertices[j]) {
          bary(i, k) = weights[j];
        }
      }
    }
  }
};

/// The element and barycentric coordinates ([point][vertex]) of the results,
/// as an evaluator takes them (LiftResults), computed on the device and
/// returned on the host.
template <int Dim, typename Result>
std::pair<Kokkos::View<Omega_h::LO*, Kokkos::HostSpace>,
          Kokkos::View<pcms::Real**, Kokkos::LayoutRight, Kokkos::HostSpace>>
LiftOnHost(Omega_h::Mesh& mesh, Kokkos::View<Result*> results)
{
  const auto n = static_cast<int>(results.extent(0));
  Kokkos::View<Omega_h::LO*> elems("elems", n);
  // LayoutRight on every backend, so the host mirror has the declared type
  Kokkos::View<pcms::Real**, Kokkos::LayoutRight> bary("bary", n, Dim + 1);
  Kokkos::parallel_for(
    "LiftOnHost", n,
    LiftResults<Dim, Result>{
      pcms::backend::omega_h::detail::UpAdjacencies<Dim>(mesh),
      pcms::backend::omega_h::detail::VerticesOfEntities<Dim>(mesh), results,
      elems, bary});
  return {Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, elems),
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, bary)};
}

/// Whether `elem` has all of `verts` among its vertices.
bool ElementHasVertices(Omega_h::Mesh& mesh, Omega_h::LO elem,
                        const std::vector<Omega_h::LO>& verts)
{
  const auto elem_verts = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  for (const auto v : verts) {
    bool found = false;
    for (int k = 0; k <= mesh.dim(); ++k) {
      found = found || elem_verts[elem * (mesh.dim() + 1) + k] == v;
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

/// The point the barycentric coordinates `bary` give in element `elem`.
template <typename BaryView>
std::vector<pcms::Real> PositionInElement(Omega_h::Mesh& mesh, Omega_h::LO elem,
                                          const BaryView& bary, int row)
{
  const int dim = mesh.dim();
  const auto elem_verts = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  const auto coords = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  std::vector<pcms::Real> x(dim, 0.0);
  for (int k = 0; k <= dim; ++k) {
    for (int d = 0; d < dim; ++d) {
      x[d] += bary(row, k) * coords[elem_verts[elem * (dim + 1) + k] * dim + d];
    }
  }
  return x;
}

/// The lowest element adjacent to entity `entity` of dimension `entity_dim`.
Omega_h::LO LowestAdjacentElement(Omega_h::Mesh& mesh, int entity_dim,
                                  Omega_h::LO entity)
{
  const auto up = mesh.ask_up(entity_dim, mesh.dim());
  const auto a2ab = Omega_h::HostRead<Omega_h::LO>(up.a2ab);
  const auto ab2b = Omega_h::HostRead<Omega_h::LO>(up.ab2b);
  Omega_h::LO lowest = ab2b[a2ab[entity]];
  for (Omega_h::LO i = a2ab[entity]; i < a2ab[entity + 1]; ++i) {
    lowest = std::min(lowest, ab2b[i]);
  }
  return lowest;
}

} // namespace

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
    auto intersection_map =
      pcms::backend::omega_h::detail::construct_intersection_map_2d(
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
    auto intersection_map =
      pcms::backend::omega_h::detail::construct_intersection_map_2d(
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
        pcms::backend::omega_h::detail::construct_intersection_map_2d(
          mesh, grid_d, ncells),
        pcms::backend::omega_h::detail::construct_intersection_map_reference_2d(
          mesh, grid_d, ncells));
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
        pcms::backend::omega_h::detail::construct_intersection_map_3d(
          mesh, grid_d, ncells),
        pcms::backend::omega_h::detail::construct_intersection_map_reference_3d(
          mesh, grid_d, ncells));
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
  GridPointSearch2D search{mesh, 10, 10, 0.01};

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
  auto [elems_h, bary_h] = LiftOnHost<2>(mesh, results);
  const auto edge_verts = Omega_h::HostRead<Omega_h::LO>(
    mesh.ask_down(Omega_h::EDGE, Omega_h::VERT).ab2b);
  const auto mesh_coords = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  SECTION("global coordinate within mesh")
  {
    {
      auto [dim, idx, parametric] = results_h(0);
      const auto elem = elems_h(0);

      CAPTURE(idx);

      REQUIRE(dim == GridPointSearch2D::Result::Dimensionality::VERTEX);
      REQUIRE(idx == 0);
      // a vertex has no coordinates of its own
      REQUIRE(parametric[0] == 0);
      REQUIRE(parametric[1] == 0);
      // the lift is the lowest element around the vertex, with all weight on
      // the vertex
      REQUIRE(elem == LowestAdjacentElement(mesh, 0, idx));
      REQUIRE(ElementHasVertices(mesh, elem, {idx}));
      const auto x = PositionInElement(mesh, elem, bary_h, 0);
      REQUIRE(x[0] == Catch::Approx(0).margin(1e-14));
      REQUIRE(x[1] == Catch::Approx(0).margin(1e-14));
    }
    {
      auto [dim, idx, parametric] = results_h(1);
      const auto elem = elems_h(1);
      REQUIRE(dim == GridPointSearch2D::Result::Dimensionality::EDGE);
      REQUIRE(idx == 156);
      // the position along the edge is that of the point's barycentric
      // coordinates on the edge's vertices renormalized, which differs from
      // the projection's by less than the tolerance
      const auto va = edge_verts[idx * 2];
      const auto vb = edge_verts[idx * 2 + 1];
      const double abx = mesh_coords[vb * 2] - mesh_coords[va * 2];
      const double aby = mesh_coords[vb * 2 + 1] - mesh_coords[va * 2 + 1];
      const double t = ((0.55 - mesh_coords[va * 2]) * abx +
                        (0.54 - mesh_coords[va * 2 + 1]) * aby) /
                       (abx * abx + aby * aby);
      REQUIRE(parametric[0] == Catch::Approx(t).margin(0.01));
      REQUIRE(parametric[1] == 0);
      // the lift is the lowest element around the edge, with the weight on the
      // edge's vertices
      REQUIRE(elem == LowestAdjacentElement(mesh, 1, idx));
      REQUIRE(ElementHasVertices(mesh, elem, {va, vb}));
      const auto x = PositionInElement(mesh, elem, bary_h, 1);
      REQUIRE(x[0] == Catch::Approx(mesh_coords[va * 2] + parametric[0] * abx));
      REQUIRE(x[1] ==
              Catch::Approx(mesh_coords[va * 2 + 1] + parametric[0] * aby));
      REQUIRE(x[0] == Catch::Approx(0.55).margin(0.02));
      REQUIRE(x[1] == Catch::Approx(0.54).margin(0.02));
    }
    {
      auto [dim, idx, parametric] = results_h(7);
      const auto elem = elems_h(7);
      REQUIRE(dim == GridPointSearch2D::Result::Dimensionality::FACE);
      REQUIRE(idx == 0);
      // a face is its own element, with the reduced coordinates as weights
      REQUIRE(elem == idx);
      REQUIRE(bary_h(7, 0) == Catch::Approx(1 - parametric[0] - parametric[1]));
      REQUIRE(bary_h(7, 1) == Catch::Approx(parametric[0]));
      REQUIRE(bary_h(7, 2) == Catch::Approx(parametric[1]));
      const auto x = PositionInElement(mesh, elem, bary_h, 7);
      REQUIRE(x[0] == Catch::Approx(0.05));
      REQUIRE(x[1] == Catch::Approx(0.02));
    }
  }
  SECTION("Global coordinate outside mesh")
  {
    auto out_of_bounds = results_h(2);
    auto top_right = results_h(3);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::VERTEX);
    REQUIRE(out_of_bounds.IsOutside());
    REQUIRE(out_of_bounds.EntityId() == top_right.entity);
    REQUIRE(elems_h(2) >= 0);
    REQUIRE(ElementHasVertices(mesh, elems_h(2), {out_of_bounds.EntityId()}));

    out_of_bounds = results_h(4);
    auto bot_left = results_h(0);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::VERTEX);
    REQUIRE(out_of_bounds.IsOutside());
    REQUIRE(out_of_bounds.EntityId() == bot_left.entity);
    REQUIRE(elems_h(4) >= 0);
    REQUIRE(ElementHasVertices(mesh, elems_h(4), {out_of_bounds.EntityId()}));

    out_of_bounds = results_h(5);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::EDGE);
    REQUIRE(out_of_bounds.IsOutside());
    REQUIRE(out_of_bounds.EntityId() == 219);
    REQUIRE(elems_h(5) >= 0);
    REQUIRE(ElementHasVertices(mesh, elems_h(5),
                               {edge_verts[219 * 2], edge_verts[219 * 2 + 1]}));
    // the outside point's position on its nearest edge is its projection
    REQUIRE(out_of_bounds.parametric[0] >= 0);
    REQUIRE(out_of_bounds.parametric[0] <= 1);

    // (0.05, -0.01) is within the edge tolerance of the bottom edge, so it is
    // located on that edge rather than reported outside.
    auto on_edge = results_h(6);
    REQUIRE(on_edge.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::EDGE);
    REQUIRE_FALSE(on_edge.IsOutside());
    REQUIRE(on_edge.entity == 0);
    REQUIRE(elems_h(6) >= 0);
    REQUIRE(
      ElementHasVertices(mesh, elems_h(6), {edge_verts[0], edge_verts[1]}));
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
    auto [ext_elems_h, ext_bary_h] = LiftOnHost<2>(mesh, ext_results);

    auto res = ext_results_h(0);
    // Must be out-of-bounds (negative entity id)
    REQUIRE(res.entity < 0);
    // The outside classification is approximate: the candidate with the
    // largest signed plane distance is put on the facets it is beyond, which
    // here is an edge at the mesh corner rather than the corner vertex.
    REQUIRE(res.dimensionality ==
            GridPointSearch2D::Result::Dimensionality::EDGE);
    // The lift still resolves to a valid face around the edge
    REQUIRE(ext_elems_h(0) >= 0);
    REQUIRE(ext_elems_h(0) == LowestAdjacentElement(mesh, 1, res.EntityId()));
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
  const auto n = pcms::detail::DivisionsPerAxisForMesh<3>(mesh.nelems());
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
    if (results_h(i).entity < 0) {
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

  SECTION("3D: the only tetrahedron")
  {
    using pcms::GridPointSearch3D;
    using Dimensionality = GridPointSearch3D::Result::Dimensionality;
    auto mesh = pcms::test::BuildReferenceTet(lib);
    const auto n = pcms::detail::DivisionsPerAxisForMesh<3>(mesh.nelems());
    GridPointSearch3D search{mesh, n, n, n};

    Kokkos::View<pcms::Real* [3]> points("points", 2);
    auto points_h = Kokkos::create_mirror_view(points);
    // inside the tet
    points_h(0, 0) = 0.1;
    points_h(0, 1) = 0.1;
    points_h(0, 2) = 0.1;
    // in the tet's bounding box but outside it: the nearest entity is a face
    points_h(1, 0) = 0.5;
    points_h(1, 1) = 0.5;
    points_h(1, 2) = 0.5;
    Kokkos::deep_copy(points, points_h);

    auto results = search(points);
    auto results_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);
    const auto owners_h = LiftOnHost<3>(mesh, results).first;

    REQUIRE(results_h(0).dimensionality == Dimensionality::REGION);
    REQUIRE_FALSE(results_h(0).IsOutside());
    REQUIRE(results_h(0).entity == 0);
    REQUIRE(results_h(0).EntityId() == 0);
    // a region is its own element
    REQUIRE(owners_h(0) == 0);

    // the nearest entity is the face opposite vertex 0; face 0 or not, the
    // point lifts to the only tetrahedron around it
    REQUIRE(results_h(1).dimensionality == Dimensionality::FACE);
    REQUIRE(results_h(1).IsOutside());
    REQUIRE(owners_h(1) == 0);

    // The grid over a single tet has no empty cell, so decode the
    // no-candidate result directly.
    Kokkos::View<GridPointSearch3D::Result*> none("none", 1);
    Kokkos::deep_copy(none,
                      GridPointSearch3D::Result{Dimensionality::NO_INTERSECT,
                                                -1, Omega_h::zero_vector<3>()});
    auto none_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, none);
    const auto none_owners_h = LiftOnHost<3>(mesh, none).first;
    REQUIRE(none_h(0).IsOutside());
    REQUIRE(none_h(0).EntityId() == -1);
    REQUIRE(none_owners_h(0) == -1);
  }

  SECTION("2D: nearest vertex 0 and an empty query cell")
  {
    using pcms::GridPointSearch2D;
    using Dimensionality = GridPointSearch2D::Result::Dimensionality;
    // Only the triangle below the diagonal, so the cells above it hold no
    // candidates.
    auto mesh = pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3}));
    const auto n = pcms::detail::DivisionsPerAxisForMesh<2>(mesh.nelems());
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
    const auto owners_h = LiftOnHost<2>(mesh, results).first;

    REQUIRE(results_h(0).dimensionality == Dimensionality::VERTEX);
    REQUIRE(results_h(0).IsOutside());
    REQUIRE(results_h(0).EntityId() == 0);
    REQUIRE(owners_h(0) == 0);

    REQUIRE(results_h(1).dimensionality == Dimensionality::NO_INTERSECT);
    REQUIRE(results_h(1).IsOutside());
    REQUIRE(results_h(1).EntityId() == -1);
    REQUIRE(owners_h(1) == -1);

    REQUIRE(results_h(2).dimensionality == Dimensionality::FACE);
    REQUIRE_FALSE(results_h(2).IsOutside());
    REQUIRE(results_h(2).entity == 0);
    REQUIRE(owners_h(2) == 0);
  }
}

namespace
{
using pcms::GridPointSearch3D;
using Dim3 = GridPointSearch3D::Result::Dimensionality;

/// The id of the entity of dimension `dim` with exactly the vertices `verts`.
Omega_h::LO FindEntity(Omega_h::Mesh& mesh, int dim,
                       std::vector<Omega_h::LO> verts)
{
  std::sort(verts.begin(), verts.end());
  const auto adj = mesh.ask_down(dim, Omega_h::VERT);
  const auto ev = Omega_h::HostRead<Omega_h::LO>(adj.ab2b);
  for (Omega_h::LO e = 0; e < mesh.nents(dim); ++e) {
    std::vector<Omega_h::LO> ids;
    for (int k = 0; k <= dim; ++k) {
      ids.push_back(ev[e * (dim + 1) + k]);
    }
    std::sort(ids.begin(), ids.end());
    if (ids == verts) {
      return e;
    }
  }
  return -1;
}

/// The weights of `weights` (by vertex id) reordered to the vertex order of
/// entity `e` of dimension `dim`, without the first: the reduced coordinates.
std::vector<double> Reduced(
  Omega_h::Mesh& mesh, int dim, Omega_h::LO e,
  const std::vector<std::pair<Omega_h::LO, double>>& w)
{
  const auto adj = mesh.ask_down(dim, Omega_h::VERT);
  const auto ev = Omega_h::HostRead<Omega_h::LO>(adj.ab2b);
  std::vector<double> out;
  for (int k = 1; k <= dim; ++k) {
    double v = 0;
    for (auto [id, weight] : w) {
      if (id == ev[e * (dim + 1) + k]) {
        v = weight;
      }
    }
    out.push_back(v);
  }
  return out;
}

Omega_h::HostRead<Omega_h::Real> HostCoords(Omega_h::Mesh& mesh)
{
  return Omega_h::HostRead<Omega_h::Real>(mesh.coords());
}

/// The point sum_k w_k x_{id_k}.
std::array<double, 3> Combine(
  Omega_h::Mesh& mesh, const std::vector<std::pair<Omega_h::LO, double>>& w)
{
  const auto x = HostCoords(mesh);
  std::array<double, 3> p{0, 0, 0};
  for (auto [id, weight] : w) {
    for (int d = 0; d < 3; ++d) {
      p[d] += weight * x[id * 3 + d];
    }
  }
  return p;
}

Kokkos::View<GridPointSearch3D::Result*, Kokkos::HostSpace> Search3D(
  const GridPointSearch3D& search,
  const std::vector<std::array<double, 3>>& pts)
{
  Kokkos::View<pcms::Real* [3]> points("points", pts.size());
  auto points_h = Kokkos::create_mirror_view(points);
  for (size_t i = 0; i < pts.size(); ++i) {
    for (int d = 0; d < 3; ++d) {
      points_h(i, d) = pts[i][d];
    }
  }
  Kokkos::deep_copy(points, points_h);
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                             search(points));
}
} // namespace

TEST_CASE("3D point search classifies points by entity",
          "[point_search][3d][classification]")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildReferenceTet(lib);
  const Omega_h::LO v0 = 0, v1 = 1, v2 = 2, v3 = 3;
  const auto edge12 = FindEntity(mesh, 1, {v1, v2});
  const auto edge03 = FindEntity(mesh, 1, {v0, v3});
  const auto face123 = FindEntity(mesh, 2, {v1, v2, v3});
  const auto face012 = FindEntity(mesh, 2, {v0, v1, v2});
  REQUIRE(edge12 >= 0);
  REQUIRE(face123 >= 0);

  const double tolerances = 1e-3;
  GridPointSearch3D search{mesh, 4, 4, 4, tolerances};

  const std::vector<std::pair<Omega_h::LO, double>> edge_w{{v1, 0.7},
                                                           {v2, 0.3}};
  const std::vector<std::pair<Omega_h::LO, double>> edge03_w{{v0, 0.4},
                                                             {v3, 0.6}};
  const std::vector<std::pair<Omega_h::LO, double>> face_w{
    {v1, 0.2}, {v2, 0.3}, {v3, 0.5}};
  const std::vector<std::pair<Omega_h::LO, double>> face012_w{
    {v0, 0.5}, {v1, 0.3}, {v2, 0.2}};
  const std::vector<std::pair<Omega_h::LO, double>> region_w{
    {v0, 0.1}, {v1, 0.2}, {v2, 0.3}, {v3, 0.4}};

  const auto off = [](std::array<double, 3> p, std::array<double, 3> dir,
                      double h) {
    const auto n =
      std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    for (int d = 0; d < 3; ++d) {
      p[d] += h * dir[d] / n;
    }
    return p;
  };
  const auto pv3 = Combine(mesh, {{v3, 1.0}});
  const auto pe = Combine(mesh, edge_w);
  const auto pe03 = Combine(mesh, edge03_w);
  const auto pf = Combine(mesh, face_w);
  const auto pf012 = Combine(mesh, face012_w);
  const auto pr = Combine(mesh, region_w);

  SECTION("on the entity")
  {
    const auto r = Search3D(search, {pv3, pe, pe03, pf, pf012, pr});
    REQUIRE(r(0).dimensionality == Dim3::VERTEX);
    REQUIRE(r(0).entity == v3);
    REQUIRE(r(0).parametric[0] == 0.0);
    REQUIRE(r(0).parametric[1] == 0.0);
    REQUIRE(r(0).parametric[2] == 0.0);

    REQUIRE(r(1).dimensionality == Dim3::EDGE);
    REQUIRE(r(1).entity == edge12);
    REQUIRE(r(1).parametric[0] ==
            Catch::Approx(Reduced(mesh, 1, edge12, edge_w)[0]));
    REQUIRE(r(1).parametric[1] == 0.0);
    REQUIRE(r(1).parametric[2] == 0.0);
    REQUIRE(r(2).dimensionality == Dim3::EDGE);
    REQUIRE(r(2).entity == edge03);
    REQUIRE(r(2).parametric[0] ==
            Catch::Approx(Reduced(mesh, 1, edge03, edge03_w)[0]));

    REQUIRE(r(3).dimensionality == Dim3::FACE);
    REQUIRE(r(3).entity == face123);
    const auto exp = Reduced(mesh, 2, face123, face_w);
    REQUIRE(r(3).parametric[0] == Catch::Approx(exp[0]));
    REQUIRE(r(3).parametric[1] == Catch::Approx(exp[1]));
    REQUIRE(r(3).parametric[2] == 0.0);
    REQUIRE(r(4).dimensionality == Dim3::FACE);
    REQUIRE(r(4).entity == face012);
    const auto exp2 = Reduced(mesh, 2, face012, face012_w);
    REQUIRE(r(4).parametric[0] == Catch::Approx(exp2[0]));
    REQUIRE(r(4).parametric[1] == Catch::Approx(exp2[1]));

    // strictly inside: the tetrahedron's barycentrics without the first
    REQUIRE(r(5).dimensionality == Dim3::REGION);
    REQUIRE(r(5).entity == 0);
    const auto exp3 = Reduced(mesh, 3, 0, region_w);
    for (int j = 0; j < 3; ++j) {
      REQUIRE(r(5).parametric[j] == Catch::Approx(exp3[j]));
    }
  }

  SECTION("just off the entity within tolerance classifies the same")
  {
    const auto r =
      Search3D(search, {off(pv3, {1, 1, 1}, 2e-4), off(pe, {0, 0, -1}, 1e-4),
                        off(pf, {1, 1, 1}, 1e-4), off(pf, {-1, -1, -1}, 1e-4)});
    REQUIRE(r(0).dimensionality == Dim3::VERTEX);
    REQUIRE(r(0).entity == v3);
    REQUIRE_FALSE(r(0).IsOutside());

    REQUIRE(r(1).dimensionality == Dim3::EDGE);
    REQUIRE(r(1).entity == edge12);
    REQUIRE(r(1).parametric[0] ==
            Catch::Approx(Reduced(mesh, 1, edge12, edge_w)[0]));
    const auto exp = Reduced(mesh, 2, face123, face_w);
    for (int i : {2, 3}) {
      REQUIRE(r(i).dimensionality == Dim3::FACE);
      REQUIRE(r(i).entity == face123);
      REQUIRE(r(i).parametric[0] == Catch::Approx(exp[0]).margin(1e-4));
      REQUIRE(r(i).parametric[1] == Catch::Approx(exp[1]).margin(1e-4));
    }
  }

  SECTION("outside the mesh: nearest entity")
  {
    GridPointSearch3D tight{mesh, 4, 4, 4};
    const auto r =
      Search3D(tight, {{-0.1, -0.1, -0.1}, {0.6, 0.6, -0.1}, {0.5, 0.5, 0.5}});
    REQUIRE(r(0).dimensionality == Dim3::VERTEX);
    REQUIRE(r(0).IsOutside());
    REQUIRE(r(0).EntityId() == v0);
    REQUIRE(r(0).entity == GridPointSearch3D::Result::OutsideId(v0));

    REQUIRE(r(1).dimensionality == Dim3::EDGE);
    REQUIRE(r(1).IsOutside());
    REQUIRE(r(1).EntityId() == edge12);
    REQUIRE(r(1).parametric[0] ==
            Catch::Approx(Reduced(mesh, 1, edge12, {{v1, 0.5}, {v2, 0.5}})[0]));

    REQUIRE(r(2).dimensionality == Dim3::FACE);
    REQUIRE(r(2).IsOutside());
    REQUIRE(r(2).EntityId() == face123);
    const auto exp =
      Reduced(mesh, 2, face123, {{v1, 1. / 3}, {v2, 1. / 3}, {v3, 1. / 3}});
    REQUIRE(r(2).parametric[0] == Catch::Approx(exp[0]));
    REQUIRE(r(2).parametric[1] == Catch::Approx(exp[1]));
  }

  SECTION("no candidates")
  {
    // two tetrahedra far apart: the cells between them hold no candidate
    Omega_h::Mesh two(&lib);
    Omega_h::build_from_elems_and_coords(
      &two, OMEGA_H_SIMPLEX, 3, Omega_h::LOs({0, 1, 2, 3, 4, 5, 6, 7}),
      Omega_h::Reals({0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                      3, 3, 3, 4, 3, 3, 3, 4, 3, 3, 3, 4}));
    GridPointSearch3D sparse{two, 8, 8, 8, tolerances};
    const auto r = Search3D(sparse, {{2.0, 2.0, 2.0}});
    REQUIRE(r(0).dimensionality == Dim3::NO_INTERSECT);
    REQUIRE(r(0).entity == -1);
    REQUIRE(r(0).IsOutside());
    for (int j = 0; j < 3; ++j) {
      REQUIRE(r(0).parametric[j] == 0.0);
    }
  }
}

TEST_CASE("3D point search matches brute force on interior points",
          "[point_search][3d][regression]")
{
  auto lib = Omega_h::Library{};
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 1, 6, 6, 6, false);
  const auto n = pcms::detail::DivisionsPerAxisForMesh<3>(mesh.nelems());
  GridPointSearch3D search{mesh, n, n, n};

  std::mt19937 gen(12345);
  std::uniform_real_distribution<double> uni(0.001, 0.999);
  std::vector<std::array<double, 3>> pts(300);
  for (auto& p : pts) {
    p = {uni(gen), uni(gen), uni(gen)};
  }
  const auto r = Search3D(search, pts);

  const auto x = HostCoords(mesh);
  const auto ev = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  for (size_t i = 0; i < pts.size(); ++i) {
    Omega_h::Vector<3> pt{pts[i][0], pts[i][1], pts[i][2]};
    Omega_h::LO expected = -1;
    Omega_h::Vector<4> bary;
    for (Omega_h::LO e = 0; e < mesh.nelems() && expected < 0; ++e) {
      Omega_h::Few<Omega_h::Vector<3>, 4> vc;
      for (int k = 0; k < 4; ++k) {
        for (int d = 0; d < 3; ++d) {
          vc[k][d] = x[ev[e * 4 + k] * 3 + d];
        }
      }
      bary = Omega_h::barycentric_from_global<3, 3>(pt, vc);
      if (Omega_h::is_barycentric_inside(bary)) {
        expected = e;
      }
    }
    REQUIRE(expected >= 0);
    CAPTURE(i);
    REQUIRE(r(i).dimensionality == Dim3::REGION);
    REQUIRE(r(i).entity == expected);
    for (int j = 0; j < 3; ++j) {
      REQUIRE(r(i).parametric[j] == Catch::Approx(bary[j + 1]));
    }
  }
}

namespace
{
/// The result of Eq. 13 of the barycentric tolerance notes, evaluated on the
/// host independently of the search: the entity a point lies on, by the zero
/// pattern of its barycentric coordinates in every element of the mesh.
template <int Dim>
class BandReference
{
public:
  struct Answer
  {
    int dim = -1;
    Omega_h::LO id = -1;
    bool outside = false;
    /// outside, with several candidates at the same distance up to rounding
    bool tied = false;
    std::array<double, Dim> xi{};
  };

  BandReference(Omega_h::Mesh& mesh, double tau) : mesh_(mesh), tau_(tau)
  {
    coords_ = HostCoords(mesh);
    const auto ev = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
    for (Omega_h::LO i = 0; i < ev.size(); ++i) {
      elem_verts_.push_back(ev[i]);
    }
    for (int k = 1; k < Dim; ++k) {
      const auto a =
        Omega_h::HostRead<Omega_h::LO>(mesh.ask_down(k, Omega_h::VERT).ab2b);
      for (Omega_h::LO i = 0; i < a.size(); ++i) {
        entity_verts_[k].push_back(a[i]);
      }
    }
  }

  /// lambda and ||grad lambda|| of a point in element `e`, from the inverse
  /// of the matrix of edge vectors
  void Barycentric(Omega_h::LO e, const std::array<double, Dim>& p,
                   std::array<double, Dim + 1>& lambda,
                   std::array<double, Dim + 1>& grad_norm) const
  {
    double a[Dim][2 * Dim];
    for (int r = 0; r < Dim; ++r) {
      for (int c = 0; c < Dim; ++c) {
        a[r][c] = X(e, c + 1)[r] - X(e, 0)[r];
        a[r][Dim + c] = r == c ? 1.0 : 0.0;
      }
    }
    for (int i = 0; i < Dim; ++i) { // Gauss-Jordan with partial pivoting
      int piv = i;
      for (int r = i + 1; r < Dim; ++r) {
        if (std::fabs(a[r][i]) > std::fabs(a[piv][i])) {
          piv = r;
        }
      }
      for (int c = 0; c < 2 * Dim; ++c) {
        std::swap(a[i][c], a[piv][c]);
      }
      const double d = a[i][i];
      for (int c = 0; c < 2 * Dim; ++c) {
        a[i][c] /= d;
      }
      for (int r = 0; r < Dim; ++r) {
        if (r != i) {
          const double f = a[r][i];
          for (int c = 0; c < 2 * Dim; ++c) {
            a[r][c] -= f * a[i][c];
          }
        }
      }
    }
    // row j of the inverse is grad lambda_{j+1}
    lambda[0] = 1.0;
    std::array<double, Dim> g0{};
    for (int j = 0; j < Dim; ++j) {
      double l = 0, n = 0;
      for (int c = 0; c < Dim; ++c) {
        const double g = a[j][Dim + c];
        l += g * (p[c] - X(e, 0)[c]);
        n += g * g;
        g0[c] -= g;
      }
      lambda[j + 1] = l;
      lambda[0] -= l;
      grad_norm[j + 1] = std::sqrt(n);
    }
    double n0 = 0;
    for (int c = 0; c < Dim; ++c) {
      n0 += g0[c] * g0[c];
    }
    grad_norm[0] = std::sqrt(n0);
  }

  /// The answer over the given elements, in increasing id order
  Answer operator()(const std::array<double, Dim>& p,
                    const std::vector<Omega_h::LO>& candidates) const
  {
    Answer best;
    Omega_h::LO nearest = -1;
    double nearest_dist = -INFINITY;
    bool tied = false;
    for (const auto e : candidates) {
      std::array<double, Dim + 1> lambda, gn;
      Barycentric(e, p, lambda, gn);
      double dist = INFINITY;
      for (int i = 0; i <= Dim; ++i) {
        dist = std::min(dist, lambda[i] / gn[i]);
      }
      if (dist < -tau_) {
        tied = tied || std::fabs(dist - nearest_dist) < 1e-9;
        if (dist > nearest_dist) {
          nearest = e;
          nearest_dist = dist;
        }
        continue;
      }
      const auto a = Classify(e, lambda, gn, false);
      if (best.id < 0 || a.dim < best.dim ||
          (a.dim == best.dim && a.id < best.id)) {
        best = a;
      }
    }
    if (best.id < 0 && nearest >= 0) {
      std::array<double, Dim + 1> lambda, gn;
      Barycentric(nearest, p, lambda, gn);
      best = Classify(nearest, lambda, gn, true);
      best.tied = tied;
    }
    return best;
  }

  std::array<double, Dim> X(Omega_h::LO e, int k) const
  {
    std::array<double, Dim> x;
    for (int d = 0; d < Dim; ++d) {
      x[d] = coords_[elem_verts_[e * (Dim + 1) + k] * Dim + d];
    }
    return x;
  }

  Omega_h::LO Vertex(Omega_h::LO e, int k) const
  {
    return elem_verts_[e * (Dim + 1) + k];
  }

  /// The point with barycentric coordinates `w` in element `e`.
  std::array<double, Dim> Point(Omega_h::LO e,
                                const std::array<double, Dim + 1>& w) const
  {
    std::array<double, Dim> p{};
    for (int k = 0; k <= Dim; ++k) {
      for (int d = 0; d < Dim; ++d) {
        p[d] += w[k] * X(e, k)[d];
      }
    }
    return p;
  }

private:
  Answer Classify(Omega_h::LO e, const std::array<double, Dim + 1>& lambda,
                  const std::array<double, Dim + 1>& gn, bool outside) const
  {
    std::vector<int> v; // local vertices not in a facet band
    int largest = 0;
    for (int i = 0; i <= Dim; ++i) {
      const double eps = tau_ * gn[i];
      const bool on = outside ? lambda[i] <= eps : std::fabs(lambda[i]) <= eps;
      if (!on) {
        v.push_back(i);
      }
      if (lambda[i] > lambda[largest]) {
        largest = i;
      }
    }
    if (v.empty()) {
      v.push_back(largest);
    }
    Answer a;
    a.outside = outside;
    a.dim = static_cast<int>(v.size()) - 1;
    std::vector<Omega_h::LO> ids;
    for (int i : v) {
      ids.push_back(Vertex(e, i));
    }
    // reduced coordinates in the vertex order of the entity
    std::vector<Omega_h::LO> order;
    if (a.dim == Dim) {
      a.id = e;
      for (int k = 0; k <= Dim; ++k) {
        order.push_back(Vertex(e, k));
      }
    } else if (a.dim == 0) {
      a.id = ids[0];
    } else {
      a.id = FindEntity(mesh_, a.dim, ids);
      for (int m = 0; m <= a.dim; ++m) {
        order.push_back(entity_verts_[a.dim][a.id * (a.dim + 1) + m]);
      }
    }
    double sum = 0;
    for (int i : v) {
      sum += std::max(lambda[i], 0.0);
    }
    for (int m = 1; m < static_cast<int>(order.size()); ++m) {
      for (int i : v) {
        if (Vertex(e, i) == order[m]) {
          a.xi[m - 1] = std::max(lambda[i], 0.0) / sum;
        }
      }
    }
    return a;
  }

  Omega_h::Mesh& mesh_;
  double tau_;
  Omega_h::HostRead<Omega_h::Real> coords_;
  std::vector<Omega_h::LO> elem_verts_;
  std::vector<Omega_h::LO> entity_verts_[Dim];
};

/// A box mesh whose interior vertices are moved by up to `shake` times the
/// spacing, so that neighbouring elements differ in size and shape.
template <int Dim>
Omega_h::Mesh BuildShakenBox(Omega_h::Library& lib, int divisions, double shake)
{
  auto mesh =
    Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, Dim == 3 ? 1 : 0,
                       divisions, divisions, Dim == 3 ? divisions : 0, false);
  std::mt19937 gen(7);
  std::uniform_real_distribution<double> uni(-shake / divisions,
                                             shake / divisions);
  const auto x = HostCoords(mesh);
  Omega_h::HostWrite<Omega_h::Real> moved(x.size());
  for (Omega_h::LO v = 0; v < mesh.nverts(); ++v) {
    bool interior = true;
    for (int d = 0; d < Dim; ++d) {
      interior = interior && x[v * Dim + d] > 1e-9 && x[v * Dim + d] < 1 - 1e-9;
    }
    for (int d = 0; d < Dim; ++d) {
      moved[v * Dim + d] = x[v * Dim + d] + (interior ? uni(gen) : 0.0);
    }
  }
  mesh.set_coords(Omega_h::Reals(moved.write()));
  return mesh;
}

template <int Dim, typename Search>
auto SearchAll(const Search& search,
               const std::vector<std::array<double, Dim>>& pts)
{
  Kokkos::View<pcms::Real* [Dim]> points("points", pts.size());
  auto points_h = Kokkos::create_mirror_view(points);
  for (size_t i = 0; i < pts.size(); ++i) {
    for (int d = 0; d < Dim; ++d) {
      points_h(i, d) = pts[i][d];
    }
  }
  Kokkos::deep_copy(points, points_h);
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                             search(points));
}

template <int Dim, typename Search>
void CheckAgainstBandReference(Omega_h::Mesh& mesh, const Search& search,
                               double tau)
{
  const BandReference<Dim> ref(mesh, tau);
  std::vector<Omega_h::LO> all(mesh.nelems());
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    all[e] = e;
  }
  std::mt19937 gen(99);
  std::uniform_real_distribution<double> uni(0.05, 1.0);
  std::uniform_int_distribution<int> pick_elem(0, mesh.nelems() - 1);
  std::uniform_int_distribution<int> pick_sign(0, 1);
  std::vector<std::array<double, Dim>> pts;
  std::vector<std::string> kinds;
  // `bands` facets at a fraction of the tolerance (negative: beyond the
  // facet plane), the others random and positive
  const auto place = [&](int bands, double fraction, const char* kind) {
    for (int n = 0; n < 150; ++n) {
      const auto e = pick_elem(gen);
      std::array<double, Dim + 1> lambda, gn;
      std::array<double, Dim + 1> w;
      std::array<bool, Dim + 1> banded{};
      for (int b = 0; b < bands; ++b) {
        int i;
        do {
          i = std::uniform_int_distribution<int>(0, Dim)(gen);
        } while (banded[i]);
        banded[i] = true;
      }
      std::array<double, Dim> any{};
      for (int d = 0; d < Dim; ++d) {
        any[d] = ref.X(e, 0)[d];
      }
      ref.Barycentric(e, any, lambda, gn);
      double fixed = 0, free_sum = 0;
      for (int i = 0; i <= Dim; ++i) {
        if (banded[i]) {
          const double sign = fraction < 0 || pick_sign(gen) ? 1 : -1;
          w[i] = sign * fraction * tau * gn[i];
          fixed += w[i];
        } else {
          w[i] = uni(gen);
          free_sum += w[i];
        }
      }
      for (int i = 0; i <= Dim; ++i) {
        if (!banded[i]) {
          w[i] *= (1.0 - fixed) / free_sum;
        }
      }
      pts.push_back(ref.Point(e, w));
      kinds.push_back(kind);
    }
  };
  place(0, 0.0, "interior");
  for (int bands = 1; bands <= Dim; ++bands) {
    place(bands, 0.5, "within the band of facets");
  }
  place(1, -1.5, "just outside the closure band");
  std::uniform_real_distribution<double> wide(-0.3, 1.3);
  for (int n = 0; n < 150; ++n) {
    std::array<double, Dim> p;
    for (auto& c : p) {
      c = wide(gen);
    }
    pts.push_back(p);
    kinds.push_back("anywhere");
  }

  const auto r = SearchAll<Dim>(search, pts);
  int dims[Dim + 2] = {};
  int outside = 0;
  for (size_t i = 0; i < pts.size(); ++i) {
    const auto a = ref(pts[i], all);
    CAPTURE(i, kinds[i]);
    if (a.id < 0) {
      REQUIRE(static_cast<int>(r(i).dimensionality) == Dim + 1);
      continue;
    }
    REQUIRE(r(i).IsOutside() == a.outside);
    if (a.outside && kinds[i] == "anywhere") {
      // far outside, the nearest entity is the one among the candidates of
      // the point's grid cell, not among all elements
      continue;
    }
    if (a.tied) {
      // equally near candidates (the facets along a flat boundary) are
      // ordered by rounding
      continue;
    }
    REQUIRE(static_cast<int>(r(i).dimensionality) == a.dim);
    REQUIRE(r(i).EntityId() == a.id);
    for (int j = 0; j < Dim; ++j) {
      REQUIRE(r(i).parametric[j] == Catch::Approx(a.xi[j]).margin(1e-8));
    }
    ++dims[a.dim];
    outside += a.outside;
  }
  // every kind of entity and the outside branch were exercised
  for (int d = 0; d <= Dim; ++d) {
    CAPTURE(d);
    REQUIRE(dims[d] > 10);
  }
  REQUIRE(outside > 10);
}
} // namespace

TEST_CASE("point search matches the facet-band rule on a non-uniform mesh",
          "[point_search][band]")
{
  auto lib = Omega_h::Library{};
  const double tau = 1e-3;
  SECTION("2D")
  {
    auto mesh = BuildShakenBox<2>(lib, 7, 0.3);
    pcms::GridPointSearch2D search{mesh, 8, 8, tau};
    CheckAgainstBandReference<2>(mesh, search, tau);
  }
  SECTION("3D")
  {
    auto mesh = BuildShakenBox<3>(lib, 4, 0.2);
    GridPointSearch3D search{mesh, 6, 6, 6, tau};
    CheckAgainstBandReference<3>(mesh, search, tau);
  }
}

namespace
{
/// Searches `pts` from per-point guessed elements (none for a negative
/// guess), on the host.
template <int Dim, typename Search>
auto SearchGuessed(const Search& search,
                   const std::vector<std::array<double, Dim>>& pts,
                   const std::vector<Omega_h::LO>& guess_values)
{
  const auto n = pts.size();
  Kokkos::View<pcms::Real* [Dim]> points("points", n);
  auto points_h = Kokkos::create_mirror_view(points);
  Kokkos::View<Omega_h::LO*> guesses("guesses", guess_values.size());
  auto guesses_h = Kokkos::create_mirror_view(guesses);
  for (size_t i = 0; i < n; ++i) {
    for (int d = 0; d < Dim; ++d) {
      points_h(i, d) = pts[i][d];
    }
  }
  for (size_t i = 0; i < guess_values.size(); ++i) {
    guesses_h(i) = guess_values[i];
  }
  Kokkos::deep_copy(points, points_h);
  Kokkos::deep_copy(guesses, guesses_h);
  Kokkos::View<typename Search::Result*> results("results", n);
  search(points, guesses, results);
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);
}

/// The sample of a mesh the guessed search is checked at: random points in a
/// box around the mesh (inside and out), and the centroids of all entities
/// (vertices, edge midpoints, face and element centers).
template <int Dim>
std::vector<std::array<double, Dim>> GuessSamplePoints(Omega_h::Mesh& mesh)
{
  std::mt19937 gen(11);
  std::uniform_real_distribution<double> uni(-0.1, 1.1);
  std::vector<std::array<double, Dim>> pts(300);
  for (auto& p : pts) {
    for (int d = 0; d < Dim; ++d) {
      p[d] = uni(gen);
    }
  }
  for (int dim = 0; dim <= Dim; ++dim) {
    const auto c =
      Omega_h::HostRead<Omega_h::Real>(pcms::get_entity_centroids(mesh, dim));
    for (Omega_h::LO e = 0; e < mesh.nents(dim); ++e) {
      std::array<double, Dim> p;
      for (int d = 0; d < Dim; ++d) {
        p[d] = c[e * Dim + d];
      }
      pts.push_back(p);
    }
  }
  return pts;
}

/// Checks that searching from guesses gives the results of the search without
/// guesses (the entity, the outside flag and the reduced coordinates), for
/// the correct elements, for elements next to them, for wrong elements, for
/// no guesses (-1 and an empty view).
template <int Dim, typename Search>
void CheckGuessedSearchMatchesCold(Omega_h::Mesh& mesh, const Search& search)
{
  const auto pts = GuessSamplePoints<Dim>(mesh);
  const auto n = static_cast<Omega_h::LO>(pts.size());
  Kokkos::View<pcms::Real* [Dim]> points("points", n);
  auto points_h = Kokkos::create_mirror_view(points);
  for (Omega_h::LO i = 0; i < n; ++i) {
    for (int d = 0; d < Dim; ++d) {
      points_h(i, d) = pts[i][d];
    }
  }
  Kokkos::deep_copy(points, points_h);
  const auto cold_device = search(points);
  const auto cold =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cold_device);
  // the elements the cold results lift to: -1 where the grid cell has none
  const auto correct = LiftOnHost<Dim>(mesh, cold_device).first;
  const auto dual = mesh.ask_dual();
  const auto offsets = Omega_h::HostRead<Omega_h::LO>(dual.a2ab);
  const auto neighbours = Omega_h::HostRead<Omega_h::LO>(dual.ab2b);
  std::mt19937 gen(5);
  std::uniform_int_distribution<Omega_h::LO> any_element(0, mesh.nelems() - 1);

  std::vector<Omega_h::LO> right(n), next_to(n), wrong(n), none(n, -1);
  for (Omega_h::LO i = 0; i < n; ++i) {
    right[i] = correct(i);
    next_to[i] = right[i];
    if (right[i] >= 0 && offsets[right[i] + 1] > offsets[right[i]]) {
      next_to[i] = neighbours[offsets[right[i]] +
                              i % (offsets[right[i] + 1] - offsets[right[i]])];
    }
    wrong[i] = any_element(gen);
  }
  const auto check = [&](const auto& results) {
    for (Omega_h::LO i = 0; i < n; ++i) {
      CAPTURE(i, pts[i][0], pts[i][1]);
      REQUIRE(results(i).dimensionality == cold(i).dimensionality);
      REQUIRE(results(i).entity == cold(i).entity);
      REQUIRE(results(i).IsOutside() == cold(i).IsOutside());
      for (int j = 0; j < Dim; ++j) {
        REQUIRE(results(i).parametric[j] ==
                Catch::Approx(cold(i).parametric[j]).margin(1e-10));
      }
    }
  };
  check(SearchGuessed<Dim>(search, pts, right));
  check(SearchGuessed<Dim>(search, pts, next_to));
  check(SearchGuessed<Dim>(search, pts, wrong));
  check(SearchGuessed<Dim>(search, pts, none));
  check(SearchGuessed<Dim>(search, pts, {}));
}
} // namespace

TEST_CASE("point search from guessed elements gives the search's results",
          "[point_search][guess]")
{
  auto lib = Omega_h::Library{};
  const double tau = 1e-8;
  SECTION("2D")
  {
    auto mesh = BuildShakenBox<2>(lib, 7, 0.3);
    pcms::GridPointSearch2D search{mesh, 8, 8, tau};
    CheckGuessedSearchMatchesCold<2>(mesh, search);
  }
  SECTION("3D")
  {
    auto mesh = BuildShakenBox<3>(lib, 4, 0.2);
    GridPointSearch3D search{mesh, 6, 6, 6, tau};
    CheckGuessedSearchMatchesCold<3>(mesh, search);
  }
}

TEST_CASE("point search band reaches tau / sin(angle / 2) at an acute corner",
          "[point_search][band]")
{
  auto lib = Omega_h::Library{};
  const double tau = 1e-3;
  auto mesh = BuildShakenBox<2>(lib, 7, 0.3);
  pcms::GridPointSearch2D search{mesh, 8, 8, tau};
  const BandReference<2> ref(mesh, tau);
  // the sharpest corner of the mesh
  double smallest = 10;
  Omega_h::LO elem = -1;
  int corner = 0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    for (int k = 0; k < 3; ++k) {
      const auto v = ref.X(e, k);
      const auto a = ref.X(e, (k + 1) % 3);
      const auto b = ref.X(e, (k + 2) % 3);
      const double ax = a[0] - v[0], ay = a[1] - v[1];
      const double bx = b[0] - v[0], by = b[1] - v[1];
      const double angle = std::acos((ax * bx + ay * by) /
                                     (std::hypot(ax, ay) * std::hypot(bx, by)));
      if (angle < smallest) {
        smallest = angle;
        elem = e;
        corner = k;
      }
    }
  }
  const double reach = tau / std::sin(smallest / 2);
  REQUIRE(reach > 2 * tau);
  const auto v = ref.X(elem, corner);
  const auto a = ref.X(elem, (corner + 1) % 3);
  const auto b = ref.X(elem, (corner + 2) % 3);
  std::array<double, 2> dir{
    (a[0] - v[0]) / std::hypot(a[0] - v[0], a[1] - v[1]) +
      (b[0] - v[0]) / std::hypot(b[0] - v[0], b[1] - v[1]),
    (a[1] - v[1]) / std::hypot(a[0] - v[0], a[1] - v[1]) +
      (b[1] - v[1]) / std::hypot(b[0] - v[0], b[1] - v[1])};
  const double len = std::hypot(dir[0], dir[1]);
  // farther from the vertex than the tolerance, inside the facet bands
  const double r = 0.5 * (tau + reach);
  const std::array<double, 2> p{v[0] + r * dir[0] / len,
                                v[1] + r * dir[1] / len};
  const auto res = SearchAll<2>(search, {p});
  REQUIRE(std::hypot(p[0] - v[0], p[1] - v[1]) > tau);
  REQUIRE(res(0).dimensionality ==
          pcms::GridPointSearch2D::Result::Dimensionality::VERTEX);
  REQUIRE(res(0).entity == ref.Vertex(elem, corner));
  // beyond the reach it is an edge or face
  const double far = 1.5 * reach;
  const auto res_far = SearchAll<2>(
    search, {{v[0] + far * dir[0] / len, v[1] + far * dir[1] / len}});
  REQUIRE(res_far(0).dimensionality !=
          pcms::GridPointSearch2D::Result::Dimensionality::VERTEX);
}

TEST_CASE("point search gives the same entity from either adjacent element",
          "[point_search][band]")
{
  auto lib = Omega_h::Library{};
  const double tau = 1e-3;
  const auto check = [&](auto dim_tag, int entity_dim, auto& mesh,
                         const auto& search) {
    constexpr int Dim = decltype(dim_tag)::value;
    const BandReference<Dim> ref(mesh, tau);
    const auto up = mesh.ask_up(entity_dim, Dim);
    const auto up_a2ab = Omega_h::HostRead<Omega_h::LO>(up.a2ab);
    const auto up_ab2b = Omega_h::HostRead<Omega_h::LO>(up.ab2b);
    const auto ev = Omega_h::HostRead<Omega_h::LO>(
      mesh.ask_down(entity_dim, Omega_h::VERT).ab2b);
    const auto coords = HostCoords(mesh);
    int checked = 0;
    for (Omega_h::LO ent = 0; ent < mesh.nents(entity_dim); ++ent) {
      if (up_a2ab[ent + 1] - up_a2ab[ent] != 2 || checked >= 40) {
        continue;
      }
      // a point on the entity at weights 1/(d+1)+..., placed from each
      // element's coordinates and moved towards that element's centroid by
      // a fraction of the tolerance
      std::array<double, Dim> on{};
      std::vector<double> w{0.5, 0.3, 0.2};
      w.resize(entity_dim + 1);
      double total = 0;
      for (auto x : w) {
        total += x;
      }
      for (int m = 0; m <= entity_dim; ++m) {
        for (int d = 0; d < Dim; ++d) {
          on[d] +=
            w[m] / total * coords[ev[ent * (entity_dim + 1) + m] * Dim + d];
        }
      }
      std::vector<std::array<double, Dim>> pts;
      for (int s = 0; s < 2; ++s) {
        const auto e = up_ab2b[up_a2ab[ent] + s];
        std::array<double, Dim> c{};
        for (int k = 0; k <= Dim; ++k) {
          for (int d = 0; d < Dim; ++d) {
            c[d] += ref.X(e, k)[d] / (Dim + 1);
          }
        }
        double n = 0;
        for (int d = 0; d < Dim; ++d) {
          n += (c[d] - on[d]) * (c[d] - on[d]);
        }
        std::array<double, Dim> p = on;
        for (int d = 0; d < Dim; ++d) {
          p[d] += 0.2 * tau * (c[d] - on[d]) / std::sqrt(n);
        }
        pts.push_back(p);
      }
      const auto r = SearchAll<Dim>(search, pts);
      CAPTURE(ent, entity_dim);
      REQUIRE(static_cast<int>(r(0).dimensionality) == entity_dim);
      REQUIRE(static_cast<int>(r(1).dimensionality) == entity_dim);
      REQUIRE(r(0).entity == ent);
      REQUIRE(r(1).entity == ent);
      for (int j = 0; j < Dim; ++j) {
        REQUIRE(r(0).parametric[j] ==
                Catch::Approx(r(1).parametric[j]).margin(1e-2));
      }
      // the parameters follow the entity's own vertex order
      for (int m = 1; m <= entity_dim; ++m) {
        REQUIRE(r(0).parametric[m - 1] ==
                Catch::Approx(w[m] / total).margin(1e-2));
      }
      ++checked;
    }
    REQUIRE(checked > 0);
  };
  SECTION("2D edges")
  {
    auto mesh = BuildShakenBox<2>(lib, 7, 0.3);
    pcms::GridPointSearch2D search{mesh, 8, 8, tau};
    check(std::integral_constant<int, 2>{}, 1, mesh, search);
  }
  SECTION("3D edges and faces")
  {
    auto mesh = BuildShakenBox<3>(lib, 4, 0.2);
    GridPointSearch3D search{mesh, 6, 6, 6, tau};
    check(std::integral_constant<int, 3>{}, 1, mesh, search);
    check(std::integral_constant<int, 3>{}, 2, mesh, search);
  }
}

TEST_CASE("uniform grid search reports the cell and the position in it",
          "[point_search][uniform_grid]")
{
  using Search = pcms::UniformGridPointSearch<2>;
  // 4 x 2 cells of size 0.5 x 1 over [0, 2] x [0, 2]
  const pcms::UniformGrid<2> grid{
    .edge_length = {2.0, 2.0}, .bot_left = {0.0, 0.0}, .divisions = {4, 2}};
  Kokkos::View<pcms::Real* [2]> points("points", 3);
  auto points_h = Kokkos::create_mirror_view(points);
  // in the cell at (1, 0) in axis order, a quarter and a half along its axes
  points_h(0, 0) = 0.5 + 0.125;
  points_h(0, 1) = 0.5;
  // on the grid's far corner
  points_h(1, 0) = 2.0;
  points_h(1, 1) = 2.0;
  // outside, a cell width beyond the near edge of the first cell
  points_h(2, 0) = -0.5;
  points_h(2, 1) = 0.5;
  Kokkos::deep_copy(points, points_h);

  const auto results = Kokkos::create_mirror_view_and_copy(
    Kokkos::HostSpace{},
    Search(grid)(
      Kokkos::View<const pcms::Real**, pcms::DeviceMemorySpace>(points)));

  const auto inside = results(0);
  REQUIRE_FALSE(inside.IsOutside());
  REQUIRE(inside.cell == grid.ClosestCellID(Omega_h::Vector<2>{0.625, 0.5}));
  REQUIRE(inside.parametric[0] == Catch::Approx(0.25));
  REQUIRE(inside.parametric[1] == Catch::Approx(0.5));

  // the grid is closed: its far corner is in its last cell, at the corner
  const auto corner = results(1);
  REQUIRE_FALSE(corner.IsOutside());
  REQUIRE(corner.cell == 7);
  REQUIRE(corner.parametric[0] == Catch::Approx(1.0));
  REQUIRE(corner.parametric[1] == Catch::Approx(1.0));

  // a point outside gets its nearest cell and its position relative to it
  const auto outside = results(2);
  REQUIRE(outside.IsOutside());
  REQUIRE(outside.CellId() == 0);
  REQUIRE(outside.cell == Search::Result::OutsideId(0));
  REQUIRE(outside.parametric[0] == Catch::Approx(-1.0));
  REQUIRE(outside.parametric[1] == Catch::Approx(0.5));
}
