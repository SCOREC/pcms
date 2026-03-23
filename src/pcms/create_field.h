#ifndef CREATE_FIELD_H_
#define CREATE_FIELD_H_

#include "adapter/meshfields/mesh_fields_adapter_layout.h"
#include "adapter/uniform_grid/uniform_grid_field.h"
#include "adapter/uniform_grid/uniform_grid_field_layout.h"
#include "field_layout.h"
#include "field.h"
#include "coordinate_system.h"
#include "utility/types.h"
#include "uniform_grid.h"

#include <Omega_h_mesh.hpp>
#include <array>
#include <memory>
#include <vector>

namespace pcms
{
std::unique_ptr<FieldLayout> CreateLagrangeLayout(
  Omega_h::Mesh& mesh, int order, int num_components,
  CoordinateSystem coordinate_system, std::string global_id_name = "global");

/**
 * \brief Create a binary field on a uniform grid indicating inside/outside mesh
 *
 * This function creates a uniform grid from the mesh and generates a binary
 * field where each grid vertex is assigned:
 *   - 1 if the vertex is inside the original mesh
 *   - 0 if the vertex is outside the original mesh
 *
 * Uses GridPointSearch to determine if a point lies within the mesh domain.
 * Currently only supports 2D meshes.
 *
 * \tparam dim Spatial dimension of the mesh (currently only dim=2 is supported)
 * \param mesh The Omega_h mesh to create the grid from
 * \param divisions Array specifying the number of cells in each dimension
 * \return std::pair containing the layout and field (layout must outlive field)
 *
 * \note The returned field has vertex-centered data with values 0.0 or 1.0.
 *       Vertex values are ordered according to the grid's internal indexing.
 *       The layout must be kept alive as long as the field is used.
 */
template <unsigned dim = 2>
std::pair<std::unique_ptr<UniformGridFieldLayout<dim>>,
          std::unique_ptr<UniformGridField<dim>>>
CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                             const std::array<LO, dim>& divisions);

/**
 * \brief Create a binary field with equal divisions in all dimensions
 *
 * Convenience function that creates a uniform grid with the same number of
 * cells in each dimension and generates the binary inside/outside field.
 *
 * \tparam dim Spatial dimension of the mesh (currently only dim=2 is supported)
 * \param mesh The Omega_h mesh to create the grid from
 * \param cells_per_dim Number of cells per dimension (same for all dimensions)
 * \return std::pair containing the layout and field (layout must outlive field)
 */
template <unsigned dim = 2>
std::pair<std::unique_ptr<UniformGridFieldLayout<dim>>,
          std::unique_ptr<UniformGridField<dim>>>
CreateUniformGridBinaryField(Omega_h::Mesh& mesh, LO cells_per_dim);

/**
 * \brief Create a binary field on a given uniform grid
 *
 * This function takes a pre-defined uniform grid and generates a binary field
 * where each grid vertex is assigned:
 *   - 1 if the vertex is inside the mesh
 *   - 0 if the vertex is outside the mesh
 *
 * This allows testing with custom grids that may extend beyond the mesh
 * boundaries.
 *
 * \tparam dim Spatial dimension (currently only dim=2 is supported)
 * \param mesh The Omega_h mesh to test against
 * \param grid The uniform grid to evaluate
 * \return std::pair containing the layout and field (layout must outlive field)
 */
template <unsigned dim = 2>
std::pair<std::unique_ptr<UniformGridFieldLayout<dim>>,
          std::unique_ptr<UniformGridField<dim>>>
CreateUniformGridBinaryFieldFromGrid(Omega_h::Mesh& mesh,
                                     UniformGrid<dim>& grid);

} // namespace pcms

#endif // CREATE_FIELD_H_
