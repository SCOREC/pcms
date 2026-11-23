#include <catch2/catch_test_macros.hpp>
#include <pcms/localization/candidate_map.hpp>

#include <vector>

using pcms::AABBox;
using pcms::LO;
using pcms::Real;

namespace
{
// The rows of a CRS, on the host.
std::vector<std::vector<LO>> Rows(const pcms::CandidateMap<2>::Bins& crs)
{
  const auto row_map =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, crs.row_map);
  const auto entries =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, crs.entries);
  std::vector<std::vector<LO>> rows(crs.numRows());
  for (size_t r = 0; r < rows.size(); ++r) {
    for (auto k = row_map(r); k < row_map(r + 1); ++k) {
      rows[r].push_back(entries(k));
    }
  }
  return rows;
}

AABBox<2> Box(Real lo_x, Real hi_x, Real lo_y, Real hi_y)
{
  return {.center = {(lo_x + hi_x) / 2, (lo_y + hi_y) / 2},
          .half_width = {(hi_x - lo_x) / 2, (hi_y - lo_y) / 2}};
}

template <typename T>
Kokkos::View<const T*> ToDevice(const std::vector<T>& host)
{
  Kokkos::View<T*> device("device", host.size());
  auto mirror = Kokkos::create_mirror_view(device);
  for (size_t i = 0; i < host.size(); ++i) {
    mirror(i) = host[i];
  }
  Kokkos::deep_copy(device, mirror);
  return device;
}
} // namespace

TEST_CASE("Shape accepts the types with a bounding box")
{
  static_assert(pcms::Shape<pcms::Simplex<2>>);
  static_assert(pcms::Shape<pcms::Simplex<3>>);
  static_assert(pcms::Shape<AABBox<3>>);
  static_assert(!pcms::Shape<pcms::Point<2>>);
  static_assert(!pcms::Shape<pcms::Segment<2>>);
  SUCCEED();
}

TEST_CASE("CandidateMap bins shapes and answers queries")
{
  // 4 x 4 unit cells over [0, 4]^2
  const pcms::UniformGrid2D grid{
    .edge_length = {4, 4}, .bot_left = {0, 0}, .divisions = {4, 4}};
  const auto cell = [&](Real x, Real y) {
    return grid.ClosestCellID(Kokkos::Array<Real, 2>{x, y});
  };
  // shape 0 inside cell (0, 0); shape 1 across cells (1, 1) and (2, 1)
  const auto shapes = ToDevice(
    std::vector{Box(0.25, 0.75, 0.25, 0.75), Box(1.25, 2.75, 1.25, 1.75)});
  const auto map = pcms::CreateCandidateMap(shapes, grid);

  REQUIRE(map.GetGrid().divisions[0] == 4);
  REQUIRE(map.GetGrid().divisions[1] == 4);

  SECTION("bins")
  {
    const auto bins = Rows(map.GetBins());
    REQUIRE(bins.size() == 16);
    REQUIRE(bins[cell(0.5, 0.5)] == std::vector<LO>{0});
    REQUIRE(bins[cell(1.5, 1.5)] == std::vector<LO>{1});
    REQUIRE(bins[cell(2.5, 1.5)] == std::vector<LO>{1});
    REQUIRE(bins[cell(3.5, 3.5)].empty());
  }

  SECTION("points")
  {
    Kokkos::View<Real* [2]> points("points", 4);
    auto points_h = Kokkos::create_mirror_view(points);
    const Real coords[4][2] = {{0.5, 0.5}, {2.6, 1.5}, {3.5, 3.5}, {-1, -1}};
    for (int i = 0; i < 4; ++i) {
      points_h(i, 0) = coords[i][0];
      points_h(i, 1) = coords[i][1];
    }
    Kokkos::deep_copy(points, points_h);
    const auto rows = Rows(map.GetCandidates(points));
    REQUIRE(rows.size() == 4);
    REQUIRE(rows[0] == std::vector<LO>{0});
    REQUIRE(rows[1] == std::vector<LO>{1});
    REQUIRE(rows[2].empty());
    // outside the grid: the nearest cell, (0, 0)
    REQUIRE(rows[3] == std::vector<LO>{0});
  }

  SECTION("boxes")
  {
    const auto boxes = ToDevice(std::vector{
      Box(0, 4, 0, 4), Box(3.2, 3.8, 3.2, 3.8), Box(1.2, 1.8, 1.2, 1.8)});
    const auto rows = Rows(map.GetCandidates(boxes));
    REQUIRE(rows.size() == 3);
    // shape 1 is in two of the cells, and is listed once
    REQUIRE(rows[0] == std::vector<LO>{0, 1});
    REQUIRE(rows[1].empty());
    REQUIRE(rows[2] == std::vector<LO>{1});
  }
}
