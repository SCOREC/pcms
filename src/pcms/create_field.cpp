#include "create_field.h"
#include "adapter/meshfields/mesh_fields_adapter2.h"
#include "adapter/meshfields/mesh_fields_adapter_layout.h"
#include "adapter/uniform_grid/uniform_grid_field.h"
#include "adapter/uniform_grid/uniform_grid_field_layout.h"
#include "uniform_grid.h"
#include "point_search.h"

#include <Kokkos_Core.hpp>
#include <utility>

namespace pcms
{

std::unique_ptr<FieldLayout> CreateLagrangeLayout(
  Omega_h::Mesh& mesh, int order, int num_components,
  CoordinateSystem coordinate_system, std::string global_id_name)
{

  std::array<int, 4> nodes_per_dim;

  switch (order) {
    case 1: nodes_per_dim = {1, 0, 0, 0}; break;
    case 2: nodes_per_dim = {1, 1, 0, 0}; break;
    default: throw std::runtime_error("Unimplemented order");
  }

  return std::make_unique<MeshFieldsAdapterLayout>(
    mesh, nodes_per_dim, num_components, coordinate_system, global_id_name);
}

template <>
std::pair<std::unique_ptr<UniformGridFieldLayout<2>>,
          std::unique_ptr<UniformGridField<2>>>
CreateUniformGridBinaryFieldFromGrid<2>(Omega_h::Mesh& mesh,
                                        UniformGrid<2>& grid)
{
  constexpr unsigned dim = 2;

  // Get total number of vertices
  const LO num_vertices = (grid.divisions[0] + 1) * (grid.divisions[1] + 1);

  // Create GridPointSearch for point-in-mesh queries
  GridPointSearch2D point_search(mesh, grid.divisions[0], grid.divisions[1]);

  // Create array of grid vertex points
  Kokkos::View<Real* [dim]> vertices("vertices", num_vertices);
  auto vertices_h = Kokkos::create_mirror_view(vertices);

  // Fill vertex coordinates
  const Real dx = grid.edge_length[0] / grid.divisions[0];
  const Real dy = grid.edge_length[1] / grid.divisions[1];

  for (LO j = 0; j <= grid.divisions[1]; ++j) {
    for (LO i = 0; i <= grid.divisions[0]; ++i) {
      const LO vertex_id = j * (grid.divisions[0] + 1) + i;
      vertices_h(vertex_id, 0) = grid.bot_left[0] + i * dx;
      vertices_h(vertex_id, 1) = grid.bot_left[1] + j * dy;
    }
  }

  // Copy to device
  Kokkos::deep_copy(vertices, vertices_h);

  // Perform point search
  auto results = point_search(vertices);
  fprintf(stderr, "Completed point-in-mesh search for %d vertices.\n",
          num_vertices);

  // Copy results back to host
  auto results_h = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(results_h, results);
  fprintf(stderr, "Copied point search results back to host.\n");

  // Create UniformGridFieldLayout and field
  auto layout = std::make_unique<UniformGridFieldLayout<2>>(
    grid, 1, CoordinateSystem::Cartesian);
  auto field = std::make_unique<UniformGridField<2>>(*layout);

  // Create binary data as Real values (0.0 or 1.0)
  Kokkos::View<Real*, HostMemorySpace> binary_data("binary_data", num_vertices);
  for (LO i = 0; i < num_vertices; ++i) {
    binary_data(i) = (results_h(i).element_id >= 0) ? 1.0 : 0.0;
  }
  fprintf(stderr, "Generated binary inside/outside data for grid vertices.\n");

  // Set the DOF holder data
  Rank1View<const Real, HostMemorySpace> data_view(binary_data.data(),
                                                   binary_data.extent(0));
  field->SetDOFHolderData(data_view);
  fprintf(stderr, "Set binary data on UniformGridField.\n");

  return {std::move(layout), std::move(field)};
}

template <>
std::pair<std::unique_ptr<UniformGridFieldLayout<2>>,
          std::unique_ptr<UniformGridField<2>>>
CreateUniformGridBinaryField<2>(Omega_h::Mesh& mesh,
                                const std::array<LO, 2>& divisions)
{
  constexpr unsigned dim = 2;

  // Create the uniform grid from the mesh
  auto grid = CreateUniformGridFromMesh<dim>(mesh, divisions);

  // Delegate to the grid-based implementation
  return CreateUniformGridBinaryFieldFromGrid<2>(mesh, grid);
}

template <>
std::pair<std::unique_ptr<UniformGridFieldLayout<2>>,
          std::unique_ptr<UniformGridField<2>>>
CreateUniformGridBinaryField<2>(Omega_h::Mesh& mesh, LO cells_per_dim)
{
  std::array<LO, 2> divisions;
  divisions.fill(cells_per_dim);
  return CreateUniformGridBinaryField<2>(mesh, divisions);
}

} // namespace pcms
