#ifndef PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
#define PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP

// Public entry points of the Omega_h backend: build PCMS objects from an
// Omega_h mesh. This header (and only backend headers) may use Omega_h types.

#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/field/field.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_mesh.hpp>

#include <array>
#include <memory>
#include <optional>
#include <string>

namespace pcms::backend::omega_h
{

/// The geometry of an Omega_h simplex mesh (2D or 3D, classified). Builds the
/// mesh's containing-element search: create one per mesh and share it. The
/// geometry owns the mesh's discretization; spaces are created on it with
/// CreateLagrangeSpace and the other space factories.
/// @param mesh the mesh; must outlive the geometry and everything created
///        from it
/// @param system the mesh's coordinate system, or null for Cartesian
/// @param options search, global id tag and ownership
[[nodiscard]] std::shared_ptr<const OmegaHGeometry> CreateGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options = {});

/// The geometry of an Omega_h mesh that shares ownership of the mesh and its
/// library: they live as long as the geometry and what holds it (the Python
/// bindings).
[[nodiscard]] std::shared_ptr<const OmegaHGeometry> CreateGeometry(
  std::shared_ptr<Omega_h::Mesh> mesh,
  std::shared_ptr<Omega_h::Library> library,
  std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options = {});

/// A point-cloud geometry whose sites are the vertices (`entity_dim` 0) or the
/// element centroids (the mesh dimension) of the mesh of an Omega_h geometry;
/// the support search of its neighbourhoods is the mesh's adjacency search.
/// Spaces made on the same geometry object share its entities, so binding the
/// cloud at their DOF holders starts from those entities. For
/// CreatePolynomialReconstructionSpace.
/// @param geometry an Omega_h geometry
/// @param entity_dim 0 or the mesh dimension
/// @param options the support search of the geometry
/// @throws pcms_error if `entity_dim` is not 0 or the mesh dimension
[[nodiscard]] std::shared_ptr<const point_cloud::PointCloudGeometry>
CreatePointCloudGeometry(std::shared_ptr<const OmegaHGeometry> geometry,
                         int entity_dim, NeighbourhoodOptions options = {});

/// 0/1 field on the vertices of a 2D uniform grid marking the vertices that
/// lie inside mesh. The grid is given, or covers the mesh's bounding box with
/// the given number of cells per axis.
Field<Real> CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                                         const UniformGrid<2>& grid);
Field<Real> CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                                         std::array<LO, 2> divisions);

/// Rank-2 view of an Omega_h array that holds dim values per entity (e.g.
/// mesh.coords()), without copying.
template <typename T>
auto MakeConstRank2View(Omega_h::Read<T> array, int dim)
{
  // Omega_h arrays are layout_right 1D arrays
  return Rank2View<const T, DeviceMemorySpace, Kokkos::layout_right>(
    array.data(), array.size() / dim, dim);
}

/**
 * \brief Create a uniform grid layout from an Omega_h mesh as a bounding box
 *
 * This function computes the bounding box of the given mesh and creates a
 * uniform Cartesian grid that covers the entire mesh domain. The grid cells
 * are uniformly distributed across each dimension.
 *
 * \tparam dim Spatial dimension of the mesh (2 or 3)
 * \param mesh The Omega_h mesh to create the grid from
 * \param divisions Array specifying the number of cells in each dimension
 * \return UniformGrid<dim> A uniform grid covering the mesh bounding box
 *
 */
template <unsigned dim = 2>
UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,
                                           const std::array<LO, dim>& divisions)
{
  auto bbox = Omega_h::get_bounding_box<dim>(&mesh);
  Kokkos::Array<Real, dim> edge_length;
  Kokkos::Array<Real, dim> bot_left;
  Kokkos::Array<LO, dim> divs;
  for (unsigned i = 0; i < dim; ++i) {
    bot_left[i] = bbox.min[i];
    edge_length[i] = bbox.max[i] - bbox.min[i];
    divs[i] = divisions[i];
  }
  return UniformGrid<dim>{
    .edge_length = edge_length, .bot_left = bot_left, .divisions = divs};
}

/**
 * \brief Create a uniform grid with equal divisions in all dimensions
 *
 * This is a convenience function that creates a uniform grid with the same
 * number of cells in each dimension.
 *
 * \tparam dim Spatial dimension of the mesh (2 or 3)
 * \param mesh The Omega_h mesh to create the grid from
 * \param cells_per_dim Number of cells per dimension (same for all dimensions)
 * \return UniformGrid<dim> A uniform grid covering the mesh bounding box
 *
 */
template <unsigned dim = 2>
UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,
                                           LO cells_per_dim)
{
  std::array<LO, dim> divisions;
  for (unsigned i = 0; i < dim; ++i)
    divisions[i] = cells_per_dim;
  return CreateUniformGridFromMesh<dim>(mesh, divisions);
}

/// Host mask for OverlapMask(size, mask): true where mask[i] != 0, for the
/// first `size` entries of mask (all of them if size is not given). Like the
/// OverlapMask constructors it replaces, mask may be longer than size.
inline Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromArray(
  Omega_h::HostRead<Omega_h::I8> mask, std::optional<size_t> size = {})
{
  const size_t n = size.value_or(mask.size());
  PCMS_ALWAYS_ASSERT(n <= static_cast<size_t>(mask.size()));
  Kokkos::View<bool*, HostMemorySpace> out("overlap_info", n);
  for (size_t i = 0; i < n; ++i) {
    out(i) = static_cast<bool>(mask[i]);
  }
  return out;
}

inline Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromArray(
  Omega_h::Read<Omega_h::I8> mask, std::optional<size_t> size = {})
{
  return OverlapMaskFromArray(Omega_h::HostRead<Omega_h::I8>(mask), size);
}

} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
