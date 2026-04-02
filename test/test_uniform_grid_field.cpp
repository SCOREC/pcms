#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Kokkos_Core.hpp>
#include "pcms/adapter/uniform_grid/uniform_grid_field_layout.h"
#include "pcms/adapter/uniform_grid/uniform_grid_field.h"
#include "pcms/utility/uniform_grid.h"
#include "Omega_h_library.hpp"
#include "Omega_h_build.hpp"
#include "pcms/transfer_field2.h"
#include "pcms/adapter/meshfields/mesh_fields_adapter_layout.h"
#include "pcms/adapter/meshfields/mesh_fields_adapter2.h"
#include "pcms/create_field.h"
#include "pcms/utility/arrays.h"
#include <cmath>

using pcms::CreateUniformGridBinaryField;
using pcms::CreateUniformGridFromMesh;

// Helper function to initialize omega_h field data with f(x,y) = x + 2*y
std::vector<pcms::Real> CreateOmegaHFieldData(
  const pcms::CoordinateView<pcms::HostMemorySpace>& coords, int num_nodes)
{
  auto coords_data = coords.GetCoordinates();
  std::vector<pcms::Real> omega_h_data(num_nodes);
  for (int i = 0; i < num_nodes; ++i) {
    pcms::Real x = coords_data(i, 0);
    pcms::Real y = coords_data(i, 1);
    omega_h_data[i] = x + 2.0 * y; // f(x,y) = x + 2y
  }
  return omega_h_data;
}

// Helper function to verify ug_field values
void VerifyUniformGridFieldValues(
  const pcms::UniformGrid<2>& grid,
  const pcms::CoordinateView<pcms::HostMemorySpace>& ug_coords,
  const pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>& ug_field_data)
{
  fprintf(stderr, "\nVerifying ug_field values:\n");
  for (int j = 0; j <= grid.divisions[1]; ++j) {
    for (int i = 0; i <= grid.divisions[0]; ++i) {
      int vertex_id = j * (grid.divisions[0] + 1) + i;
      pcms::Real x = ug_coords.GetCoordinates()(vertex_id, 0);
      pcms::Real y = ug_coords.GetCoordinates()(vertex_id, 1);
      pcms::Real expected = x + 2.0 * y;
      pcms::Real actual = ug_field_data(vertex_id);
      fprintf(
        stderr,
        "ug_field Vertex (%d, %d) at (%.2f, %.2f): expected %.4f, got %.4f\n",
        i, j, x, y, expected, actual);
      REQUIRE(std::abs(expected - actual) <= 1e-10);
    }
  }
}

// Helper function to verify mask field values
void VerifyMaskFieldValues(const pcms::UniformGrid<2>& grid,
                           const pcms::UniformGridField<2>& mask_field)
{
  fprintf(stderr, "\nVerifying mask_field values:\n");
  auto mask_data = mask_field.GetDOFHolderData();
  for (int j = 0; j <= grid.divisions[1]; ++j) {
    for (int i = 0; i <= grid.divisions[0]; ++i) {
      int vertex_id = j * (grid.divisions[0] + 1) + i;
      pcms::Real mask_value = mask_data(vertex_id);
      fprintf(stderr, "Vertex (%d, %d): mask = %.0f\n", i, j, mask_value);
      REQUIRE(mask_value == 1.0);
    }
  }
}

TEST_CASE("UniformGrid field creation")
{
  // Create a simple 2D uniform grid
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {5, 5};

  // Create field layout with 1 component (scalar field)
  pcms::UniformGridFieldLayout<2> layout(grid, 1,
                                         pcms::CoordinateSystem::Cartesian);

  REQUIRE(layout.GetNumComponents() == 1);
  REQUIRE(layout.GetNumOwnedDofHolder() == 36); // (5+1)x(5+1) = 36 vertices
  REQUIRE(layout.GetNumGlobalDofHolder() == 36);
  REQUIRE_FALSE(layout.IsDistributed());

  // Create field
  auto field = layout.CreateFieldReal();
  REQUIRE(field != nullptr);
}

TEST_CASE("UniformGrid field data operations", "[uniform_grid_field]")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {4, 4}; // 4x4 = 16 cells

  pcms::UniformGridFieldLayout<2> layout(grid, 1,
                                         pcms::CoordinateSystem::Cartesian);
  auto field = layout.CreateFieldReal();

  // Initialize field data with vertex indices (5x5 = 25 vertices for 4x4 cells)
  std::vector<pcms::Real> data(25);
  for (size_t i = 0; i < 25; ++i) {
    data[i] = static_cast<pcms::Real>(i);
  }

  auto data_view = pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(
    data.data(), data.size());
  field->SetDOFHolderData(data_view);

  // Retrieve data
  auto retrieved = field->GetDOFHolderData();
  REQUIRE(retrieved.size() == 25);

  for (size_t i = 0; i < 25; ++i) {
    REQUIRE(retrieved[i] == static_cast<pcms::Real>(i));
  }
}

TEST_CASE("UniformGrid field evaluation - piecewise constant")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {2, 2}; // 2x2 = 4 cells

  pcms::UniformGridFieldLayout<2> layout(grid, 1,
                                         pcms::CoordinateSystem::Cartesian);
  auto field = layout.CreateFieldReal();

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
  auto data_view = pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(
    data.data(), data.size());
  field->SetDOFHolderData(data_view);

  // Evaluate at cell centers
  std::vector<pcms::Real> eval_coords = {
    2.5, 2.5, // Cell 0 center
    7.5, 2.5, // Cell 1 center
    2.5, 7.5, // Cell 2 center
    7.5, 7.5  // Cell 3 center
  };

  auto coords_view = pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
    eval_coords.data(), 4, 2);
  auto coord_view = pcms::CoordinateView<pcms::HostMemorySpace>(
    pcms::CoordinateSystem::Cartesian, coords_view);

  auto hint = field->GetLocalizationHint(coord_view);

  std::vector<pcms::Real> results(4);
  auto results_view =
    pcms::Rank1View<pcms::Real, pcms::HostMemorySpace>(results.data(), 4);
  auto results_field_view =
    pcms::FieldDataView<pcms::Real, pcms::HostMemorySpace>(
      results_view, pcms::CoordinateSystem::Cartesian);

  field->Evaluate(hint, results_field_view);

  // Check results - interpolated from vertices
  // Cell 0 center (2.5, 2.5): avg of v0,v1,v3,v4 = (1.0+1.5+2.0+2.5)/4 = 1.75
  // Cell 1 center (7.5, 2.5): avg of v1,v2,v4,v5 = (1.5+2.0+2.5+3.0)/4 = 2.25
  // Cell 2 center (2.5, 7.5): avg of v3,v4,v6,v7 = (2.0+2.5+3.0+3.5)/4 = 2.75
  // Cell 3 center (7.5, 7.5): avg of v4,v5,v7,v8 = (2.5+3.0+3.5+4.0)/4 = 3.25
  REQUIRE(std::abs(results[0] - 1.75) < 1e-10);
  REQUIRE(std::abs(results[1] - 2.25) < 1e-10);
  REQUIRE(std::abs(results[2] - 2.75) < 1e-10);
  REQUIRE(std::abs(results[3] - 3.25) < 1e-10);
}

TEST_CASE("UniformGrid field serialization")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {3, 3}; // 9 cells

  pcms::UniformGridFieldLayout<2> layout(grid, 1,
                                         pcms::CoordinateSystem::Cartesian);
  auto field = layout.CreateFieldReal();

  // Set field data (4x4 = 16 vertices for 3x3 cells)
  std::vector<pcms::Real> data(16);
  for (size_t i = 0; i < 16; ++i) {
    data[i] = static_cast<pcms::Real>(i * 10);
  }

  auto data_view = pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(
    data.data(), data.size());
  field->SetDOFHolderData(data_view);

  // Create identity permutation
  std::vector<pcms::LO> permutation(16);
  for (size_t i = 0; i < 16; ++i) {
    permutation[i] = i;
  }

  auto perm_view = pcms::Rank1View<const pcms::LO, pcms::HostMemorySpace>(
    permutation.data(), 16);

  // Serialize
  std::vector<pcms::Real> buffer(16);
  auto buffer_view =
    pcms::Rank1View<pcms::Real, pcms::HostMemorySpace>(buffer.data(), 16);
  int size = field->Serialize(buffer_view, perm_view);

  REQUIRE(size == 16);

  // Verify serialized data
  for (size_t i = 0; i < 16; ++i) {
    REQUIRE(buffer[i] == data[i]);
  }

  // Create new field and deserialize
  auto field2 = layout.CreateFieldReal();
  auto buffer_const_view =
    pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(buffer.data(), 16);
  field2->Deserialize(buffer_const_view, perm_view);

  // Verify deserialized data
  auto retrieved = field2->GetDOFHolderData();
  for (size_t i = 0; i < 16; ++i) {
    REQUIRE(retrieved[i] == data[i]);
  }
}

TEST_CASE("UniformGrid field copy")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {2, 2}; // 2x2 grid

  pcms::UniformGridFieldLayout<2> layout(grid, 1,
                                         pcms::CoordinateSystem::Cartesian);
  auto field = layout.CreateFieldReal();

  // Set vertex values with f(x,y) = x + y at 3x3 vertex positions
  // Vertices at: (0,0), (5,0), (10,0), (0,5), (5,5), (10,5), (0,10), (5,10),
  // (10,10)
  std::vector<pcms::Real> data = {
    0.0,  5.0,  10.0, // y=0:  v0(0,0)=0,   v1(5,0)=5,   v2(10,0)=10
    5.0,  10.0, 15.0, // y=5:  v3(0,5)=5,   v4(5,5)=10,  v5(10,5)=15
    10.0, 15.0, 20.0  // y=10: v6(0,10)=10, v7(5,10)=15, v8(10,10)=20
  };
  auto data_view = pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(
    data.data(), data.size());
  field->SetDOFHolderData(data_view);

  // Test 1: Evaluate at cell centers
  std::vector<pcms::Real> eval_coords = {
    2.5, 2.5, // Cell 0 center
    7.5, 2.5, // Cell 1 center
    2.5, 7.5, // Cell 2 center
    7.5, 7.5  // Cell 3 center
  };

  auto coords_view = pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
    eval_coords.data(), 4, 2);
  auto coord_view = pcms::CoordinateView<pcms::HostMemorySpace>(
    pcms::CoordinateSystem::Cartesian, coords_view);

  auto hint = field->GetLocalizationHint(coord_view);

  std::vector<pcms::Real> results(coord_view.GetCoordinates().size() / 2);
  auto results_view =
    pcms::Rank1View<pcms::Real, pcms::HostMemorySpace>(results.data(), 4);
  auto results_field_view =
    pcms::FieldDataView<pcms::Real, pcms::HostMemorySpace>(
      results_view, pcms::CoordinateSystem::Cartesian);

  field->Evaluate(hint, results_field_view);

  REQUIRE(std::abs(results[0] - 5.0) < 1e-10);
  REQUIRE(std::abs(results[1] - 10.0) < 1e-10);
  REQUIRE(std::abs(results[2] - 10.0) < 1e-10);
  REQUIRE(std::abs(results[3] - 15.0) < 1e-10);

  // Test 2: test copy_field2
  auto field2 = layout.CreateFieldReal();
  pcms::copy_field2(*field, *field2);

  auto copied_data = field2->GetDOFHolderData();
  REQUIRE(copied_data.size() == data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    REQUIRE(copied_data[i] == data[i]);
  }
}

TEST_CASE("Transfer from OmegaH field to UniformGrid field")
{
  Omega_h::Library lib;

  // Create a simple omega_h mesh (2x2 box)
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 2,
                                 2, 0, false);

  // Create OmegaH field layout with linear elements
  auto omega_h_layout =
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto omega_h_field = omega_h_layout->CreateFieldReal();

  // Initialize omega_h field with a simple function f(x,y) = x + 2*y
  auto coords = omega_h_layout->GetDOFHolderCoordinates();
  int num_nodes = omega_h_layout->GetNumOwnedDofHolder();
  std::vector<pcms::Real> omega_h_data =
    CreateOmegaHFieldData(coords, num_nodes);

  auto omega_h_data_view =
    pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(
      omega_h_data.data(), omega_h_data.size());
  omega_h_field->SetDOFHolderData(omega_h_data_view);

  // Create a uniform grid field covering the same domain [0,1] x [0,1]
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {2, 2}; // 4x4 grid for finer resolution

  pcms::UniformGridFieldLayout<2> ug_layout(grid, 1,
                                            pcms::CoordinateSystem::Cartesian);
  auto ug_field = ug_layout.CreateFieldReal();

  // Transfer from omega_h field to uniform grid field using interpolation
  auto coords_interpolation = ug_layout.GetDOFHolderCoordinates();
  std::vector<pcms::Real> evaluation(
    coords_interpolation.GetCoordinates().size() / 2);
  auto evaluation_view = pcms::Rank1View<pcms::Real, pcms::HostMemorySpace>(
    evaluation.data(), evaluation.size());
  pcms::FieldDataView<pcms::Real, pcms::HostMemorySpace> data_view{
    evaluation_view, omega_h_field->GetCoordinateSystem()};
  auto locale = omega_h_field->GetLocalizationHint(coords_interpolation);
  omega_h_field->Evaluate(locale, data_view);
  auto evaluation_view_const =
    pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(evaluation.data(),
                                                             evaluation.size());
  ug_field->SetDOFHolderData(evaluation_view_const);

  // Verify the transferred data at uniform grid vertices
  auto transferred_data = ug_field->GetDOFHolderData();
  auto ug_coords = ug_layout.GetDOFHolderCoordinates();
  auto ug_coords_data = ug_coords.GetCoordinates();
  int num_ug_nodes = ug_layout.GetNumOwnedDofHolder(); // 5x5 = 25 vertices

  // Check a few sample points
  for (int i = 0; i < num_ug_nodes; ++i) {
    pcms::Real x = ug_coords_data(i, 0);
    pcms::Real y = ug_coords_data(i, 1);
    pcms::Real expected = x + 2.0 * y;
    pcms::Real actual = transferred_data[i];

    // Allow some tolerance for interpolation
    REQUIRE(std::abs(actual - expected) < 1e-6);
  }
}

TEST_CASE("Create binary field from uniform grid")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  SECTION("Simple 2D box mesh - all vertices inside")
  {
    // Create a mesh that fills the domain [0,1] x [0,1]
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 10,
                                   10, 0, false);

    // Create a 5x5 grid (coarser than mesh)
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {5, 5});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 36); // (5+1) * (5+1) = 36 vertices

    // All grid vertices should be inside the mesh
    pcms::Real sum = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      sum += field_data(i);
    }
    REQUIRE(sum == 36.0);
  }

  SECTION("Binary field with custom divisions")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, false);

    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {10, 8});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 99); // (10+1) * (8+1) = 99 vertices

    // Most vertices should be inside
    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);
    REQUIRE(inside_count <= 99.0);
  }

  SECTION("Binary field with equal divisions convenience function")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 5, 5, 0, false);

    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, 8);
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 81); // (8+1) * (8+1) = 81 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Fine grid over coarse mesh")
  {
    // Create a simple mesh
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 2.0, 2.0, 0.0, 4, 4, 0, false);

    // Create a fine grid (20x20)
    auto grid = CreateUniformGridFromMesh<2>(mesh, {20, 20});
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {20, 20});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 441); // (20+1) * (20+1) = 441 vertices

    // Verify consistency: check some specific vertices
    // Center vertex should be inside (vertex at i=10, j=10)
    int center_idx = 10 * 21 + 10; // 21 vertices per row
    REQUIRE(field_data(center_idx) == 1.0);

    // Corner vertex should be inside
    int corner_idx = 0; // vertex (0, 0)
    REQUIRE(field_data(corner_idx) == 1.0);
  }

  SECTION("Verify field values are binary")
  {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 5, 5, 0, false);

    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, 10);
    auto field_data = field->GetDOFHolderData();

    // All values should be 0 or 1
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      pcms::Real val = field_data(i);
      REQUIRE((val == 0.0 || val == 1.0));
    }
  }

  SECTION("Grid extends beyond mesh - vertices outside should be marked 0")
  {
    // Create a small mesh in the center of a domain
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 0.5, 0.5, 0.0, 5, 5, 0, false);

    // The mesh occupies [0, 0.5] x [0, 0.5]
    // Create a grid that would cover this
    auto grid = CreateUniformGridFromMesh<2>(mesh, {10, 10});
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {10, 10});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 121); // (10+1) * (10+1) = 121 vertices

    // All vertices should be inside since grid is exactly on mesh bbox
    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Test with different aspect ratio")
  {
    // Create a rectangular mesh
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 3.0, 1.0, 0.0, 12, 4,
                                   0, false);

    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {30, 10});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 341); // (30+1) * (10+1) = 341 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);

    // Calculate percentage inside
    double inside_percent = 100.0 * inside_count / field_data.extent(0);
    // Most vertices should be inside
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
    auto grid = CreateUniformGridFromMesh<2>(mesh, {8, 8});
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {8, 8});
    auto field_data = field->GetDOFHolderData();

    // Get field value for a specific vertex (middle vertex at i=4, j=4)
    pcms::LO vertex_id = 4 * 9 + 4;        // 9 = (8+1) vertices per row
    REQUIRE(field_data(vertex_id) == 1.0); // Should be inside

    // Compute vertex position
    pcms::Real dx = grid.edge_length[0] / grid.divisions[0];
    pcms::Real dy = grid.edge_length[1] / grid.divisions[1];
    pcms::Real x = grid.bot_left[0] + 4 * dx;
    pcms::Real y = grid.bot_left[1] + 4 * dy;

    // Verify it's roughly in the middle
    REQUIRE(x == Catch::Approx(0.5).margin(0.1));
    REQUIRE(y == Catch::Approx(0.5).margin(0.1));
  }

  SECTION("Count vertices by region")
  {
    auto grid = CreateUniformGridFromMesh<2>(mesh, 10);
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, 10);
    auto field_data = field->GetDOFHolderData();

    // Count vertices in different quadrants
    int q1 = 0, q2 = 0, q3 = 0, q4 = 0; // quadrants

    pcms::Real dx = grid.edge_length[0] / grid.divisions[0];
    pcms::Real dy = grid.edge_length[1] / grid.divisions[1];

    for (int j = 0; j <= grid.divisions[1]; ++j) {
      for (int i = 0; i <= grid.divisions[0]; ++i) {
        pcms::LO vertex_id = j * (grid.divisions[0] + 1) + i;
        if (field_data(vertex_id) == 1.0) {
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

    // All quadrants should have some inside vertices
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

    // Create a very fine grid
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, 50);
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 2601); // (50+1) * (50+1) = 2601 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Coarse grid")
  {
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 10,
                                   10, 0, false);

    // Very coarse grid
    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {2, 2});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 9); // (2+1) * (2+1) = 9 vertices

    // All vertices should be inside for this configuration
    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Non-square domain")
  {
    auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 5.0, 2.0, 0.0, 20, 8,
                                   0, false);

    auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, {25, 10});
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 286); // (25+1) * (10+1) = 286 vertices

    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    REQUIRE(inside_count > 0.0);
  }

  SECTION("Grid larger than mesh - vertices outside marked as 0")
  {
    // Create a mesh covering [0, 0.5] x [0, 0.5]
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 0.5, 0.5, 0.0, 5, 5, 0, false);

    // Manually create a larger grid covering [0, 1] x [0, 1]
    pcms::UniformGrid<2> grid;
    grid.edge_length = {1.0, 1.0};
    grid.bot_left = {0.0, 0.0};
    grid.divisions = {10, 10};

    // Create binary field on the larger grid
    auto [layout, field] =
      pcms::CreateUniformGridBinaryFieldFromGrid<2>(mesh, grid);
    auto field_data = field->GetDOFHolderData();

    REQUIRE(field_data.extent(0) == 121); // (10+1) * (10+1) = 121 vertices

    // Count vertices inside and outside
    pcms::Real inside_count = 0.0;
    for (size_t i = 0; i < field_data.extent(0); ++i) {
      inside_count += field_data(i);
    }
    pcms::Real outside_count = field_data.extent(0) - inside_count;

    // Should have both inside (1) and outside (0) vertices
    REQUIRE(inside_count > 0.0);
    REQUIRE(outside_count > 0.0);

    // Check specific vertices
    // Bottom-left corner should be inside (mesh covers [0, 0.5])
    int corner_id = 0 * 11 + 0; // vertex (0, 0)
    REQUIRE(field_data(corner_id) == 1.0);

    // Top-right corner should be outside (mesh ends at 0.5)
    corner_id = 10 * 11 + 10; // vertex (10, 10)
    REQUIRE(field_data(corner_id) == 0.0);

    // Vertex at i=6, j=6 (coords ~0.6, ~0.6) should be outside
    corner_id = 6 * 11 + 6;
    REQUIRE(field_data(corner_id) == 0.0);

    // Vertex at i=2, j=2 (coords ~0.2, ~0.2) should be inside
    auto center_id = 2 * 11 + 2;
    REQUIRE(field_data(center_id) == 1.0);
  }
}

TEST_CASE("UniformGrid workflow")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  // Create a simple 2D box mesh: 1.0 x 1.0 domain with 4x4 elements
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, false);
  auto grid = pcms::CreateUniformGridFromMesh<2>(mesh, {4, 4});
  auto [mask_layout, mask_field] =
    pcms::CreateUniformGridBinaryField<2>(mesh, {4, 4});

  // Create OmegaH field layout with linear elements
  auto omega_h_layout =
    pcms::CreateLagrangeLayout(mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto omega_h_field = omega_h_layout->CreateFieldReal();

  // Initialize omega_h field with a simple function f(x,y) = x + 2*y
  auto coords = omega_h_layout->GetDOFHolderCoordinates();
  int num_nodes = omega_h_layout->GetNumOwnedDofHolder();
  std::vector<pcms::Real> omega_h_data =
    CreateOmegaHFieldData(coords, num_nodes);

  auto omega_h_data_view =
    pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace>(
      omega_h_data.data(), omega_h_data.size());
  omega_h_field->SetDOFHolderData(omega_h_data_view);

  // Create uniform grid field layout
  pcms::UniformGridFieldLayout<2> ug_layout(grid, 1,
                                            pcms::CoordinateSystem::Cartesian);
  auto ug_field = ug_layout.CreateFieldReal();

  // Transfer from omega_h field to uniform grid field using interpolation
  pcms::interpolate_field2(*omega_h_field, *ug_field);
  auto ug_coords = ug_layout.GetDOFHolderCoordinates();

  // Verify ug_field values directly from the field object
  auto ug_field_data = ug_field->GetDOFHolderData();
  VerifyUniformGridFieldValues(grid, ug_coords, ug_field_data);

  // Verify mask field values
  VerifyMaskFieldValues(grid, *mask_field);
}
