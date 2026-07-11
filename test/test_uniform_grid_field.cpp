#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Kokkos_Core.hpp>
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/field/evaluator/uniform_grid.h"
#include "pcms/field/uniform_grid_binary_field.h"
#include "pcms/field/data/simple.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/uniform_grid.h"
#include "Omega_h_library.hpp"
#include "Omega_h_build.hpp"
#include "pcms/transfer/copy.h"
#include "pcms/transfer/interpolator.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/arrays.h"
#include "field_test_utils.h"
#include <cmath>

using pcms::CreateUniformGridFromMesh;

TEST_CASE("UniformGridDiscretization SameEntities: identical grids")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {4, 4};

  pcms::UniformGridFieldLayout<2> layout_a(grid, 1,
                                           pcms::CoordinateSystem::Cartesian);
  pcms::UniformGridFieldLayout<2> layout_b(grid, 1,
                                           pcms::CoordinateSystem::Cartesian);

  auto disc_a = layout_a.GetDiscretization();
  auto disc_b = layout_b.GetDiscretization();

  REQUIRE(disc_a != nullptr);
  REQUIRE(disc_b != nullptr);
  REQUIRE(disc_a->SameEntities(*disc_b));
}

TEST_CASE("UniformGridDiscretization SameEntities: different grids")
{
  pcms::UniformGrid<2> grid_a;
  grid_a.bot_left = {0.0, 0.0};
  grid_a.edge_length = {1.0, 1.0};
  grid_a.divisions = {4, 4};

  pcms::UniformGrid<2> grid_b;
  grid_b.bot_left = {0.0, 0.0};
  grid_b.edge_length = {1.0, 1.0};
  grid_b.divisions = {8, 8};

  pcms::UniformGridFieldLayout<2> layout_a(grid_a, 1,
                                           pcms::CoordinateSystem::Cartesian);
  pcms::UniformGridFieldLayout<2> layout_b(grid_b, 1,
                                           pcms::CoordinateSystem::Cartesian);

  auto disc_a = layout_a.GetDiscretization();
  auto disc_b = layout_b.GetDiscretization();

  REQUIRE_FALSE(disc_a->SameEntities(*disc_b));
}

// Helper to verify ug_field values against f(x,y) = x + 2*y.
void VerifyUniformGridFieldValues(
  const pcms::UniformGrid<2>& grid,
  const pcms::CoordinateView<pcms::HostMemorySpace>& ug_coords,
  const pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>& ug_field_data)
{
  for (int j = 0; j <= grid.divisions[1]; ++j) {
    for (int i = 0; i <= grid.divisions[0]; ++i) {
      int vertex_id = j * (grid.divisions[0] + 1) + i;
      pcms::Real x = ug_coords.GetValues()(vertex_id, 0);
      pcms::Real y = ug_coords.GetValues()(vertex_id, 1);
      pcms::Real expected = x + 2.0 * y;
      pcms::Real actual = ug_field_data[vertex_id];
      REQUIRE(std::abs(expected - actual) <= 1e-10);
    }
  }
}

// Helper to verify binary mask field (all values == 1.0).
void VerifyMaskFieldValues(const pcms::UniformGrid<2>& grid,
                           const pcms::Field<pcms::Real>& mask_field)
{
  auto mask_data = pcms::FlattenToRank1View(mask_field.GetDOFHolderDataHost());
  for (int j = 0; j <= grid.divisions[1]; ++j) {
    for (int i = 0; i <= grid.divisions[0]; ++i) {
      int vertex_id = j * (grid.divisions[0] + 1) + i;
      REQUIRE(mask_data[vertex_id] == 1.0);
    }
  }
}

TEST_CASE("UniformGrid field creation")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {5, 5};

  auto layout = std::make_shared<pcms::UniformGridFieldLayout<2>>(
    grid, 1, pcms::CoordinateSystem::Cartesian);

  REQUIRE(layout->GetNumComponents() == 1);
  REQUIRE(layout->GetNumOwnedDofHolder() == 36); // (5+1)x(5+1) = 36 vertices
  REQUIRE(layout->GetNumGlobalDofHolder() == 36);
  REQUIRE_FALSE(layout->IsDistributed());

  auto field_space = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field = field_space.CreateField<pcms::Real>(pcms::FieldMetadata{});
  REQUIRE(field.GetDOFHolderDataHost().size() ==
          static_cast<size_t>(layout->OwnedSize()));
}

TEST_CASE("UniformGrid order-0 field creation and evaluation")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {2, 2};

  auto layout = std::make_shared<pcms::UniformGridFieldLayout<2>>(
    grid, 1, pcms::CoordinateSystem::Cartesian, 0);
  auto field_space = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian, 0);
  auto field = field_space.CreateField<pcms::Real>(pcms::FieldMetadata{});
  pcms::UniformGridEvaluatorFactory<2> eval_factory(layout);

  REQUIRE(layout->GetOrder() == 0);
  REQUIRE(layout->GetNumOwnedDofHolder() == 4);

  auto coords_device = layout->GetDOFHolderCoordinates().GetValues();
  auto coords = pcms::test::CopyCoordinatesToHost(coords_device, 4, 2);

  REQUIRE(coords(0, 0) == Catch::Approx(2.5));
  REQUIRE(coords(0, 1) == Catch::Approx(2.5));
  REQUIRE(coords(3, 0) == Catch::Approx(7.5));
  REQUIRE(coords(3, 1) == Catch::Approx(7.5));

  std::vector<pcms::Real> data = {1.0, 2.0, 3.0, 4.0};
  field.SetDOFHolderDataHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      data.data(), static_cast<pcms::LO>(data.size()), 1));

  std::vector<pcms::Real> eval_coords = {1.0, 1.0, 9.0, 1.0,
                                         1.0, 9.0, 9.0, 9.0};
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    eval_coords, pcms::CoordinateSystem::Cartesian);
  auto evaluator = eval_factory.CreatePointEvaluator(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<pcms::Real*, pcms::HostMemorySpace> results_host("results_host",
                                                                4);
  Kokkos::View<pcms::Real*, pcms::DeviceMemorySpace> results_device(
    "results_device", 4);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<pcms::DeviceMemorySpace>;
  pcms::Rank2View<pcms::Real, pcms::DeviceMemorySpace, LayoutPolicy> out(
    results_device.data(), 4, 1);
  evaluator->Evaluate(field, out);
  Kokkos::deep_copy(results_host, results_device);

  REQUIRE(results_host(0) == Catch::Approx(1.0));
  REQUIRE(results_host(1) == Catch::Approx(2.0));
  REQUIRE(results_host(2) == Catch::Approx(3.0));
  REQUIRE(results_host(3) == Catch::Approx(4.0));
}

TEST_CASE("UniformGrid field data operations", "[uniform_grid_field]")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {4, 4};

  auto layout = std::make_shared<pcms::UniformGridFieldLayout<2>>(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field_space = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field = field_space.CreateField<pcms::Real>(pcms::FieldMetadata{});

  std::vector<pcms::Real> data(25);
  for (size_t i = 0; i < 25; ++i)
    data[i] = static_cast<pcms::Real>(i);

  field.SetDOFHolderDataHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      data.data(), static_cast<pcms::LO>(data.size()), 1));

  auto retrieved = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  REQUIRE(retrieved.size() == 25);
  for (size_t i = 0; i < 25; ++i)
    REQUIRE(retrieved[i] == static_cast<pcms::Real>(i));
}

TEST_CASE("UniformGrid field evaluation - piecewise constant")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {2, 2};

  auto layout = std::make_shared<pcms::UniformGridFieldLayout<2>>(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field_space = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field = field_space.CreateField<pcms::Real>(pcms::FieldMetadata{});
  pcms::UniformGridEvaluatorFactory<2> eval_factory(layout);

  // Set vertex values for a 2x2 cell grid (3x3 = 9 vertices)
  // Vertex layout:
  //   v6---v7---v8
  //   |  2 |  3 |
  //   v3---v4---v5
  //   |  0 |  1 |
  //   v0---v1---v2
  std::vector<pcms::Real> data = {
    1.0, 1.5, 2.0, // v0, v1, v2 (bottom row, y=0)
    2.0, 2.5, 3.0, // v3, v4, v5 (middle row, y=5)
    3.0, 3.5, 4.0  // v6, v7, v8 (top row, y=10)
  };
  field.SetDOFHolderDataHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      data.data(), static_cast<pcms::LO>(data.size()), 1));

  std::vector<pcms::Real> eval_coords = {
    2.5, 2.5, // Cell 0 center
    7.5, 2.5, // Cell 1 center
    2.5, 7.5, // Cell 2 center
    7.5, 7.5  // Cell 3 center
  };
  auto device_coords = pcms::test::CreateDeviceCoordinateView(
    eval_coords, pcms::CoordinateSystem::Cartesian);
  auto evaluator = eval_factory.CreatePointEvaluator(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<pcms::Real*, pcms::HostMemorySpace> results_host("results_host",
                                                                4);
  Kokkos::View<pcms::Real*, pcms::DeviceMemorySpace> results_device(
    "results_device", 4);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<pcms::DeviceMemorySpace>;
  pcms::Rank2View<pcms::Real, pcms::DeviceMemorySpace, LayoutPolicy> out(
    results_device.data(), 4, 1);
  evaluator->Evaluate(field, out);
  Kokkos::deep_copy(results_host, results_device);

  // Check results - interpolated from vertices
  // Cell 0 center (2.5, 2.5): avg of v0,v1,v3,v4 = (1.0+1.5+2.0+2.5)/4 = 1.75
  // Cell 1 center (7.5, 2.5): avg of v1,v2,v4,v5 = (1.5+2.0+2.5+3.0)/4 = 2.25
  // Cell 2 center (2.5, 7.5): avg of v3,v4,v6,v7 = (2.0+2.5+3.0+3.5)/4 = 2.75
  // Cell 3 center (7.5, 7.5): avg of v4,v5,v7,v8 = (2.5+3.0+3.5+4.0)/4 = 3.25
  REQUIRE(std::abs(results_host(0) - 1.75) < 1e-10);
  REQUIRE(std::abs(results_host(1) - 2.25) < 1e-10);
  REQUIRE(std::abs(results_host(2) - 2.75) < 1e-10);
  REQUIRE(std::abs(results_host(3) - 3.25) < 1e-10);
}

TEST_CASE("UniformGrid field serialization")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {3, 3};

  auto layout = std::make_shared<pcms::UniformGridFieldLayout<2>>(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field_space = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field = field_space.CreateField<pcms::Real>(pcms::FieldMetadata{});

  std::vector<pcms::Real> data(16);
  for (size_t i = 0; i < 16; ++i)
    data[i] = static_cast<pcms::Real>(i * 10);

  field.SetDOFHolderDataHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      data.data(), static_cast<pcms::LO>(data.size()), 1));

  pcms::test::CheckSerializeDeserialize(field);
}

TEST_CASE("UniformGrid field copy")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {2, 2};

  auto layout = std::make_shared<pcms::UniformGridFieldLayout<2>>(
    grid, 1, pcms::CoordinateSystem::Cartesian);

  // Set vertex values with f(x,y) = x + y at 3x3 vertex positions
  // Vertices at: (0,0), (5,0), (10,0), (0,5), (5,5), (10,5), (0,10), (5,10),
  // (10,10)
  std::vector<pcms::Real> data = {
    0.0,  5.0,  10.0, // y=0:  v0(0,0)=0,   v1(5,0)=5,   v2(10,0)=10
    5.0,  10.0, 15.0, // y=5:  v3(0,5)=5,   v4(5,5)=10,  v5(10,5)=15
    10.0, 15.0, 20.0  // y=10: v6(0,10)=10, v7(5,10)=15, v8(10,10)=20
  };

  auto factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<pcms::Real>(pcms::FieldMetadata{});
  field.SetDOFHolderDataHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      data.data(), static_cast<pcms::LO>(data.size()), 1));

  auto field2 = factory.CreateField<pcms::Real>(pcms::FieldMetadata{});
  pcms::Copy<pcms::Real> copy(factory, factory);
  copy.Apply(field, field2);

  auto copied_data = pcms::FlattenToRank1View(field2.GetDOFHolderDataHost());
  REQUIRE(copied_data.size() == data.size());
  for (size_t i = 0; i < data.size(); ++i)
    REQUIRE(copied_data[i] == data[i]);
}

TEST_CASE("Transfer from OmegaH field to UniformGrid field")
{
  Omega_h::Library lib;
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 2,
                                 2, 0, false);

  auto omega_h_factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto omega_h_field =
    omega_h_factory.CreateField<pcms::Real>(pcms::FieldMetadata{});
  pcms::test::SetField(
    omega_h_field,
    OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

  pcms::UniformGrid<2> grid;
  grid.edge_length = {1.0, 1.0};
  grid.bot_left = {0.0, 0.0};
  grid.divisions = {2, 2};
  auto ug_factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto ug_field = ug_factory.CreateField<pcms::Real>(pcms::FieldMetadata{});

  pcms::Interpolator<pcms::Real> interp(omega_h_factory, ug_factory);
  interp.Apply(omega_h_field, ug_field);

  auto transferred_data =
    pcms::FlattenToRank1View(ug_field.GetDOFHolderDataHost());
  auto ug_coords = ug_factory.GetLayout()->GetDOFHolderCoordinates();
  int num_ug_nodes = ug_factory.GetLayout()->GetNumOwnedDofHolder();

  // set up_coords to host
  auto ug_coords_host =
    pcms::test::CopyCoordinatesToHost(ug_coords.GetValues(), num_ug_nodes, 2);

  for (int i = 0; i < num_ug_nodes; ++i) {
    pcms::Real x = ug_coords_host(i, 0);
    pcms::Real y = ug_coords_host(i, 1);
    pcms::Real expected = x + 2.0 * y;
    REQUIRE(std::abs(transferred_data[i] - expected) < 1e-6);
  }
}

TEST_CASE("Create binary field from uniform grid")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  SECTION("Simple 2D box mesh - all vertices inside")
  {
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 10,
                                   10, 0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{5, 5});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 36); // (5+1) * (5+1) = 36 vertices

    pcms::Real sum = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i) {
      REQUIRE((field_data[i] == 0.0 || field_data[i] == 1.0));
      sum += field_data[i];
    }
    REQUIRE(sum == 36.0);
  }

  SECTION("Binary field with custom divisions")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{10, 8});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 99); // (10+1) * (8+1) = 99 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i)
      inside_count += field_data[i];
    REQUIRE(inside_count > 0.0);
    REQUIRE(inside_count <= 99.0);
  }

  SECTION("Verify field values are binary")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{10, 10});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    for (size_t i = 0; i < field_data.size(); ++i) {
      REQUIRE((field_data[i] == 0.0 || field_data[i] == 1.0));
    }
  }

  SECTION("Grid larger than mesh - vertices outside should be marked 0")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 0.5, 0.5, 0.0, 5, 5, 0, false);

    pcms::UniformGrid<2> grid;
    grid.edge_length = {1.0, 1.0};
    grid.bot_left = {0.0, 0.0};
    grid.divisions = {10, 10};

    auto [layout, field] = pcms::CreateUniformGridBinaryField<2>(mesh, grid);
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 121); // (10+1) * (10+1) = 121 vertices

    // Count vertices inside and outside
    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i)
      inside_count += field_data[i];
    pcms::Real outside_count = field_data.size() - inside_count;

    // Should have both inside (1) and outside (0) vertices
    REQUIRE(inside_count > 0.0);
    REQUIRE(outside_count > 0.0);

    int corner_id = 0 * 11 + 0;
    REQUIRE(field_data[corner_id] == 1.0);

    corner_id = 10 * 11 + 10;
    REQUIRE(field_data[corner_id] == 0.0);

    corner_id = 6 * 11 + 6;
    REQUIRE(field_data[corner_id] == 0.0);

    auto center_id = 2 * 11 + 2;
    REQUIRE(field_data[center_id] == 1.0);
  }

  SECTION("Fine grid over coarse mesh")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 2.0, 2.0, 0.0, 4, 4, 0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{20, 20});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 441); // (20+1) * (20+1) = 441 vertices

    // Verify consistency: check some specific vertices
    // Center vertex should be inside (vertex at i=10, j=10)
    int center_idx = 10 * 21 + 10; // 21 vertices per row
    REQUIRE(field_data[center_idx] == 1.0);

    // Corner vertex should be inside
    int corner_idx = 0; // vertex (0 , 0)
    REQUIRE(field_data[corner_idx] == 1.0);
  }

  SECTION("Test with different aspect ratio")
  {
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 3.0, 1.0, 0.0, 12, 4,
                                   0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{30, 10});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 341); // (30+1) * (10+1) = 341 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i)
      inside_count += field_data[i];
    REQUIRE(inside_count > 0.0);

    double inside_percent = 100.0 * inside_count / field_data.size();
    REQUIRE(inside_percent > 50.0);
  }
}

TEST_CASE("Binary field integration with grid methods")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 10, 10, 0, false);

  SECTION("Query field value at specific grid vertex")
  {
    auto grid = CreateUniformGridFromMesh<2>(mesh, std::array{8, 8});
    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{8, 8});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    // Get field value for a specific vertex (middle vertex at i=4, j=4)
    pcms::LO vertex_id = 4 * 9 + 4;
    REQUIRE(field_data[vertex_id] == 1.0);

    pcms::Real dx = grid.edge_length[0] / grid.divisions[0];
    pcms::Real dy = grid.edge_length[1] / grid.divisions[1];
    pcms::Real x = grid.bot_left[0] + 4 * dx;
    pcms::Real y = grid.bot_left[1] + 4 * dy;

    REQUIRE(x == Catch::Approx(0.5).margin(0.1));
    REQUIRE(y == Catch::Approx(0.5).margin(0.1));
  }

  SECTION("Count vertices by region")
  {
    auto grid = CreateUniformGridFromMesh<2>(mesh, std::array{10, 10});
    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{10, 10});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    int q1 = 0, q2 = 0, q3 = 0, q4 = 0;

    pcms::Real dx = grid.edge_length[0] / grid.divisions[0];
    pcms::Real dy = grid.edge_length[1] / grid.divisions[1];

    for (int j = 0; j <= grid.divisions[1]; ++j) {
      for (int i = 0; i <= grid.divisions[0]; ++i) {
        pcms::LO vertex_id = j * (grid.divisions[0] + 1) + i;
        if (field_data[vertex_id] == 1.0) {
          pcms::Real x = grid.bot_left[0] + i * dx;
          pcms::Real y = grid.bot_left[1] + j * dy;

          if (x < 0.5 && y < 0.5)
            q1++;
          else if (x >= 0.5 && y < 0.5)
            q2++;
          else if (x < 0.5 && y >= 0.5)
            q3++;
          else
            q4++;
        }
      }
    }

    REQUIRE(q1 > 0);
    REQUIRE(q2 > 0);
    REQUIRE(q3 > 0);
    REQUIRE(q4 > 0);

    // For 11x11 vertices, boundary at x=0.5, y=0.5 splits asymmetrically:
    // q1 (x<0.5, y<0.5): 5x5=25, q2 (x>=0.5, y<0.5): 6x5=30
    // q3 (x<0.5, y>=0.5): 5x6=30, q4 (x>=0.5, y>=0.5): 6x6=36
    int total = q1 + q2 + q3 + q4;
    REQUIRE(total == 121); // All vertices inside
    REQUIRE(q1 == 25);
    REQUIRE(q2 == 30);
    REQUIRE(q3 == 30);
    REQUIRE(q4 == 36);
  }
}

TEST_CASE("Performance and edge cases")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  SECTION("Very fine grid")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 5, 5, 0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{50, 50});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 2601); // (50+1) * (50+1) = 2601 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i)
      inside_count += field_data[i];
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Coarse grid")
  {
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 10,
                                   10, 0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{2, 2});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 9); // (2+1) * (2+1) = 9 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i)
      inside_count += field_data[i];
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Non-square domain")
  {
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 5.0, 2.0, 0.0, 20, 8,
                                   0, false);

    auto [layout, field] =
      pcms::CreateUniformGridBinaryField<2>(mesh, std::array{25, 10});
    auto field_data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());

    REQUIRE(field_data.size() == 286); // (25+1) * (10+1) = 286 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.size(); ++i)
      inside_count += field_data[i];
    REQUIRE(inside_count > 0.0);
  }
}

TEST_CASE("UniformGrid workflow")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, false);
  auto grid = pcms::CreateUniformGridFromMesh<2>(mesh, {4, 4});

  auto omega_h_factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto omega_h_field =
    omega_h_factory.CreateField<pcms::Real>(pcms::FieldMetadata{});
  pcms::test::SetField(
    omega_h_field,
    OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

  auto ug_factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian);
  auto ug_field = ug_factory.CreateField<pcms::Real>(pcms::FieldMetadata{});

  auto [mask_layout, mask_field] =
    pcms::CreateUniformGridBinaryField<2>(mesh, grid);

  pcms::Interpolator<pcms::Real> interp(omega_h_factory, ug_factory);
  interp.Apply(omega_h_field, ug_field);
  auto ug_coords_device_view =
    ug_factory.GetLayout()->GetDOFHolderCoordinates().GetValues();
  auto ug_coords_host_view =
    pcms::test::CopyCoordinatesToHost(ug_coords_device_view, 25, 2);

  auto ug_field_data_device = ug_field.GetDOFHolderData();
  Kokkos::View<pcms::Real*, pcms::DeviceMemorySpace> ug_field_data_device_view(
    "", 25);
  Kokkos::parallel_for(
    "CopyFieldDataToView", 25, KOKKOS_LAMBDA(int i) {
      ug_field_data_device_view(i) = ug_field_data_device(i, 0);
    });
  auto ug_field_data_host_view =
    Kokkos::View<pcms::Real*, pcms::HostMemorySpace>("", 25);
  Kokkos::deep_copy(ug_field_data_host_view, ug_field_data_device_view);

  pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace> ug_coords(
    ug_coords_host_view.data(), 25, 2);
  pcms::CoordinateView<pcms::HostMemorySpace> ug_coords_view(
    pcms::CoordinateSystem::Cartesian, ug_coords);
  pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace> ug_field_data(
    ug_field_data_host_view.data(), 25);
  VerifyUniformGridFieldValues(grid, ug_coords_view, ug_field_data);

  VerifyMaskFieldValues(grid, mask_field);
}
