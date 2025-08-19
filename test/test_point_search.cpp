#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/point_search.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

using pcms::AABBox;
using pcms::barycentric_from_global;
using pcms::Uniform2DGrid;

TEST_CASE("global to local")
{
  Omega_h::Matrix<2, 3> coords{{0.0}, {1, 0}, {0.5, 1}};
  SECTION("check verts")
  {
    for (int i = 0; i < coords.size(); ++i) {
      auto xi = barycentric_from_global(coords[i], coords);
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
    auto xi = barycentric_from_global(point, coords);
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
      auto num_candidates = intersection_map.row_map(i + 1) - intersection_map.row_map(i);
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
//extern Omega_h::Library omega_h_library;

TEST_CASE("construct intersection map")
{
  //auto world = omega_h_library.world();
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  REQUIRE(mesh.dim() == 2);
  SECTION("grid bbox overlap")
  {
    Kokkos::View<Uniform2DGrid[1]> grid_d("uniform grid");
    auto grid_h = Kokkos::create_mirror_view(grid_d);
    grid_h(0) = Uniform2DGrid{.edge_length{1, 1}, .bot_left = {0, 0}, .divisions = {10, 10}};
    Kokkos::deep_copy(grid_d, grid_h);
    auto intersection_map = pcms::detail::construct_intersection_map(mesh, grid_d, grid_h(0).GetNumCells());
    // assert(cudaSuccess == cudaDeviceSynchronize());
    REQUIRE(intersection_map.numRows() == 100);
    REQUIRE(num_candidates_within_range(intersection_map, 2, 16));
  }
  SECTION("fine grid")
  {
    Kokkos::View<Uniform2DGrid[1]> grid_d("uniform grid");
    auto grid_h = Kokkos::create_mirror_view(grid_d);
    grid_h(0) = Uniform2DGrid{.edge_length{1, 1}, .bot_left = {0, 0}, .divisions = {60, 60}};
    Kokkos::deep_copy(grid_d, grid_h);
    // require number of candidates is >=1 and <=6
    auto intersection_map = pcms::detail::construct_intersection_map(mesh, grid_d, grid_h(0).GetNumCells());

    REQUIRE(intersection_map.numRows() == 3600);
    REQUIRE(num_candidates_within_range(intersection_map, 1, 6));
  }
}
TEST_CASE("uniform grid search") {
  using pcms::GridPointSearch;
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  GridPointSearch search{mesh,10,10};
  Kokkos::View<pcms::Real*[2]> points("test_points", 7);
  //Kokkos::View<pcms::Real*[2]> points("test_points", 1);
  auto points_h = Kokkos::create_mirror_view(points);
  points_h(0,0) = 0;
  points_h(0,1) = 0;
  points_h(1,0) = 0.55;
  points_h(1,1) = 0.54;
  points_h(2,0) = 100;
  points_h(2,1) = 100;
  points_h(3,0) = 1;
  points_h(3,1) = 1;
  points_h(4,0) = -1;
  points_h(4,1) = -1;
  points_h(5, 0) = 1.01;
  points_h(5, 1) = 0.95;
  points_h(6, 0) = 0.05;
  points_h(6, 1) = -0.01;
  Kokkos::deep_copy(points, points_h);
  auto results = search(points);
  auto results_h = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(results_h, results);
  SECTION("global coordinate within mesh")
  {
    {
      auto [dim, idx,coords] = results_h(0);
      REQUIRE(dim == GridPointSearch::Result::Dimensionality::FACE);
      REQUIRE(idx == 0);
      REQUIRE(coords[0] == Catch::Approx(1));
      REQUIRE(coords[1] == Catch::Approx(0));
      REQUIRE(coords[2] == Catch::Approx(0));
    }
    {
      auto [dim, idx,coords] = results_h(1);
      REQUIRE(dim == GridPointSearch::Result::Dimensionality::FACE);
      REQUIRE(idx == 91);
      REQUIRE(coords[0] == Catch::Approx(0.5));
      REQUIRE(coords[1] == Catch::Approx(0.1));
      REQUIRE(coords[2] == Catch::Approx(0.4));
    }
  }
  // feature needs to be added
  SECTION("Global coordinate outside mesh", "[!mayfail]") {
    auto out_of_bounds = results_h(2);
    auto top_right = results_h(3);
    REQUIRE(out_of_bounds.dimensionality == GridPointSearch::Result::Dimensionality::VERTEX);
    REQUIRE(-1*out_of_bounds.tri_id == top_right.tri_id);

    out_of_bounds = results_h(4);
    auto bot_left = results_h(0);
    REQUIRE(out_of_bounds.dimensionality == GridPointSearch::Result::Dimensionality::VERTEX);
    REQUIRE(-1*out_of_bounds.tri_id == bot_left.tri_id);

    out_of_bounds = results_h(5);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch::Result::Dimensionality::EDGE);
    REQUIRE(-1 * out_of_bounds.tri_id == top_right.tri_id);

    out_of_bounds = results_h(6);
    REQUIRE(out_of_bounds.dimensionality ==
            GridPointSearch::Result::Dimensionality::EDGE);
    REQUIRE(-1 * out_of_bounds.tri_id == bot_left.tri_id);
  }
}

TEST_CASE("radial_intersects_bbox 2D") {
  using pcms::detail::radial_intersects_bbox;
  using pcms::AABBox;

  AABBox<2> box{.center = {0.0, 0.0}, .half_width = {0.5, 0.5}};
  double cutoff_squared = 0.25;

  SECTION("Point outside cutoff") {
    double pt[2] = {1.5, 0.0};
    REQUIRE_FALSE(radial_intersects_bbox<2>(pt, box, cutoff_squared));
  }

  SECTION("Point on edge of cutoff") {
    double pt[2] = {1.0, 0.0};
    REQUIRE(radial_intersects_bbox<2>(pt, box, 1.0));
  }

  SECTION("Point inside AABB") {
    double pt[2] = {0.0, 0.0};
    REQUIRE(radial_intersects_bbox<2>(pt, box, cutoff_squared));
  }

  SECTION("Point close to corner") {
    double pt[2] = {0.8, 0.8};
    REQUIRE_FALSE(radial_intersects_bbox<2>(pt, box, 0.1));
    REQUIRE(radial_intersects_bbox<2>(pt, box, 0.64));
  }

  SECTION("Radial region partially overlaps box") {
    double pt[2] = {1.0, 0.0};
    REQUIRE(radial_intersects_bbox<2>(pt, box, 0.3)); // cutoff_squared=0.3 (~0.55 radius)
  }

  SECTION("Point on AABB edge") {
    double pt[2] = {0.5, 0.0}; // Exactly on the box edge
    REQUIRE(radial_intersects_bbox<2>(pt, box, 0.25));
  }

  SECTION("Negative coordinate check") {
    double pt[2] = {-1.5, 0.0};
    REQUIRE_FALSE(radial_intersects_bbox<2>(pt, box, 0.25));
  }

  SECTION("Excess in one axis only") {
    double pt[2] = {1.0, 0.0}; // Excess in x-axis but within y-axis
    REQUIRE(radial_intersects_bbox<2>(pt, box, 1.0));
  }
}

TEST_CASE("GridRadialNeighborFunctor 2D") {
  using namespace pcms;
  using detail::GridRadialNeighborFunctor;

  const int dim = 2;
  const int num_sources = 4;
  const int num_targets = 2;
  double cutoff = 1.0;

  // Allocate Views in the default memory space (device if enabled)
  Kokkos::View<double**> sources("sources", num_sources, dim);
  Kokkos::View<double**> targets("targets", num_targets, dim);
  
  // Initialize via host mirrors
  auto hsources = Kokkos::create_mirror_view(sources);
  auto htargets = Kokkos::create_mirror_view(targets);

  // Source points (2x2 grid)
  hsources(0,0) = 0.0; hsources(0,1) = 0.0;
  hsources(1,0) = 1.0; hsources(1,1) = 0.0;
  hsources(2,0) = 0.0; hsources(2,1) = 1.0;
  hsources(3,0) = 1.0; hsources(3,1) = 1.0;

  // Target points
  htargets(0,0) = 0.5; htargets(0,1) = 0.5;  // Center
  htargets(1,0) = 2.0; htargets(1,1) = 2.0;  // Outside

  Kokkos::deep_copy(sources, hsources);
  Kokkos::deep_copy(targets, htargets);

  // Grid setup matching production code
  UniformGrid<2> grid;
  grid.bot_left = {-0.5, -0.5};
  grid.edge_length = {2.0, 2.0};
  grid.divisions = {2, 2};

  Kokkos::View<UniformGrid<2>*> grid_view("grid", 1);
  auto hgrid = Kokkos::create_mirror_view(grid_view);
  hgrid(0) = grid;
  Kokkos::deep_copy(grid_view, hgrid);

  // Compute cell_size and copy to device
  std::array<double, 2> cell_size_host;
  for (int d = 0; d < 2; ++d)
    cell_size_host[d] = grid.edge_length[d] / grid.divisions[d];

  Kokkos::View<double[2]> cell_size_view("cell_size_view");
  auto h_cell_size_view = Kokkos::create_mirror_view(cell_size_view);
  for (int d = 0; d < 2; ++d)
    h_cell_size_view(d) = cell_size_host[d];
  Kokkos::deep_copy(cell_size_view, h_cell_size_view);

  // Cell data initialization (all sources in cell 0)
  Kokkos::View<LO*> cell_ptrs("cell_ptrs", 5);
  Kokkos::View<LO*> cell_indices("cell_indices", 4);
  
  auto hptrs = Kokkos::create_mirror_view(cell_ptrs);
  auto hidx = Kokkos::create_mirror_view(cell_indices);
  hptrs(0) = 0; hptrs(1) = 4; hptrs(2) = 4; hptrs(3) = 4; hptrs(4) = 4;
  for (int i=0; i<4; ++i) hidx(i) = i;
  
  Kokkos::deep_copy(cell_ptrs, hptrs);
  Kokkos::deep_copy(cell_indices, hidx);

  // Updated functor with cell_size_view
  GridRadialNeighborFunctor<2> functor(
    targets, sources, grid_view, cell_ptrs, cell_indices, cutoff, 4, cell_size_view
  );

  SECTION("Target inside should find 4 neighbors") {
    Kokkos::View<LO*> d_neighbors("neighbors", 4);
    Kokkos::View<LO> d_count("count");

    Kokkos::parallel_for("test_inside", 1, KOKKOS_LAMBDA(const int) {
      d_count() = functor(0, d_neighbors.data());
    });

    auto h_count = Kokkos::create_mirror_view(d_count);
    auto h_neighbors = Kokkos::create_mirror_view(d_neighbors);
    Kokkos::deep_copy(h_count, d_count);
    Kokkos::deep_copy(h_neighbors, d_neighbors);

    REQUIRE(h_count() == 4);
  }

  SECTION("Target outside should find no neighbors") {
    Kokkos::View<LO*> d_neighbors("neighbors", 4);
    Kokkos::View<LO> d_count("count");

    Kokkos::parallel_for("test_outside", 1, KOKKOS_LAMBDA(const int) {
      d_count() = functor(1, d_neighbors.data());
    });

    auto h_count = Kokkos::create_mirror_view(d_count);
    Kokkos::deep_copy(h_count, d_count);

    REQUIRE(h_count() == 0);
  }
}