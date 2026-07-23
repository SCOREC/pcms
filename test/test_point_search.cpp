#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/localization/point_search.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

using pcms::AABBox;
using pcms::Uniform2DGrid;

TEST_CASE("global to local")
{
  Omega_h::Matrix<2, 3> coords{{0.0}, {1, 0}, {0.5, 1}};
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
TEST_CASE("uniform grid search")
{
  using pcms::GridPointSearch2D;
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  auto tolerances =
    GridPointSearch2D::PointSearchTolerances("point search 2d tolerances", 2);
  auto tolerances_h = Kokkos::create_mirror_view(tolerances);
  tolerances_h(0) = 0.01;
  tolerances_h(1) = 0.01;
  Kokkos::deep_copy(tolerances, tolerances_h);

  GridPointSearch2D search{mesh, 10, 10, tolerances};

  Kokkos::View<pcms::Real**> points("test_points", 8, 2);
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
  auto results = search.apply(
	pcms::CoordinateView(pcms::CoordinateSystem::Cartesian, 
		pcms::MakeConstRank2View(points)));
  auto result_dims_h = Kokkos::create_mirror_view(results.dimensionalities);
  Kokkos::deep_copy(result_dims_h, results.dimensionalities);
  auto result_ids_h = Kokkos::create_mirror_view(results.element_ids);
  Kokkos::deep_copy(result_ids_h, results.element_ids);
  auto result_coords_h = Kokkos::create_mirror_view(results.parametric_coords);
  Kokkos::deep_copy(result_coords_h, results.parametric_coords);
  SECTION("global coordinate within mesh")
  {
    {
	  auto dim = result_dims_h(0);
	  auto idx = result_ids_h(0);
	  auto coords = Omega_h::Vector<3>{result_coords_h(0,0), result_coords_h(0,1), result_coords_h(0,2)};
      const auto face_idx = pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(result_dims_h(0)), result_ids_h(0));

      CAPTURE(idx);

      REQUIRE(dim == GridPointSearch2D::Results::Dimensionality::VERTEX);
      REQUIRE(idx == 0);
      REQUIRE(face_idx >= 0);
      REQUIRE(coords[0] == Catch::Approx(1));
      REQUIRE(coords[1] == Catch::Approx(0));
      REQUIRE(coords[2] == Catch::Approx(0));
    }
    {
	  auto dim = result_dims_h(1);
	  auto idx = result_ids_h(1);
	  auto coords = Omega_h::Vector<3>{result_coords_h(1,0), result_coords_h(1,1), result_coords_h(1,2)};
      const auto face_idx = pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(result_dims_h(1)), result_ids_h(1));
      REQUIRE(dim == GridPointSearch2D::Results::Dimensionality::EDGE);
      REQUIRE(idx == 156);
      REQUIRE(face_idx >= 0);
      REQUIRE(coords[0] == Catch::Approx(0.5));
      REQUIRE(coords[1] == Catch::Approx(0.1));
      REQUIRE(coords[2] == Catch::Approx(0.4));
    }
    {
	  auto dim = result_dims_h(7);
	  auto idx = result_ids_h(7);
      const auto face_idx = pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(result_dims_h(7)), result_ids_h(7));
      REQUIRE(dim == GridPointSearch2D::Results::Dimensionality::FACE);
      REQUIRE(idx == 0);
      REQUIRE(face_idx >= 0);
    }
  }
  // feature needs to be added
  SECTION("Global coordinate outside mesh", "[!mayfail]")
  {
	auto out_of_bounds_dim = result_dims_h(2);
    auto out_of_bounds_id = result_ids_h(2);
    auto top_right_id = result_ids_h(3);
    REQUIRE(out_of_bounds_dim ==
            GridPointSearch2D::Results::Dimensionality::VERTEX);
    REQUIRE(-1 * out_of_bounds_id == top_right_id);
    REQUIRE(pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(out_of_bounds_dim), -out_of_bounds_id) >= 0);

    out_of_bounds_dim = result_dims_h(4);
    out_of_bounds_id = result_ids_h(4);
    auto bot_left_id = result_ids_h(0);
    REQUIRE(out_of_bounds_dim ==
            GridPointSearch2D::Results::Dimensionality::VERTEX);
    REQUIRE(-1 * out_of_bounds_id == bot_left_id);
    REQUIRE(pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(out_of_bounds_dim), -out_of_bounds_id) >= 0);

	out_of_bounds_dim = result_dims_h(5);
    out_of_bounds_id = result_ids_h(5);
    REQUIRE(out_of_bounds_dim ==
            GridPointSearch2D::Results::Dimensionality::EDGE);
    REQUIRE(out_of_bounds_id == -219);
    REQUIRE(pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(out_of_bounds_dim), -out_of_bounds_id) >= 0);

    out_of_bounds_dim = result_dims_h(6);
    out_of_bounds_id = result_ids_h(6);
    REQUIRE(out_of_bounds_dim ==
            GridPointSearch2D::Results::Dimensionality::EDGE);
    REQUIRE(-1 * out_of_bounds_id == bot_left_id);
    REQUIRE(pcms::GetOwningElementId(mesh, mesh.dim(), static_cast<int>(out_of_bounds_dim), -out_of_bounds_id) >= 0);
  }
  SECTION("point on extension of an edge")
  {
    Kokkos::View<pcms::Real**> ext_points("ext_test_points", 1, 2);
    auto ext_h = Kokkos::create_mirror_view(ext_points);
    ext_h(0, 0) = 1.5;
    ext_h(0, 1) = 0.0;
    Kokkos::deep_copy(ext_points, ext_h);
    auto ext_results = search.apply(pcms::CoordinateView(pcms::CoordinateSystem::Cartesian, 
		pcms::MakeConstRank2View(ext_points)));

	auto result_dims_h = Kokkos::create_mirror_view(ext_results.dimensionalities);
	Kokkos::deep_copy(result_dims_h, ext_results.dimensionalities);
	auto result_ids_h = Kokkos::create_mirror_view(ext_results.element_ids);
	Kokkos::deep_copy(result_ids_h, ext_results.element_ids);
    // Must be out-of-bounds (negative element id)
    REQUIRE(result_ids_h(0) < 0);
    // Dimensionality must be VERTEX — the nearest entity is the
    // rightmost bottom vertex (1.0, 0).
    REQUIRE(result_dims_h(0) ==
            GridPointSearch2D::Results::Dimensionality::VERTEX);
    // Owning element should still resolve to a valid face
    REQUIRE(pcms::GetOwningElementId(mesh, mesh.dim(), 0, -result_ids_h(0)) >= 0);
  }
}
