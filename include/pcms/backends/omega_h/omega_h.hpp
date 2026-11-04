#ifndef PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
#define PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP

// Public entry points of the Omega_h backend: build PCMS objects from an
// Omega_h mesh. This header (and only backend headers) may use Omega_h types.

#include "pcms/field/field.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>

#include <array>
#include <memory>
#include <string>

namespace pcms::backend::omega_h
{

/// Options of an Omega_h geometry.
struct OmegaHGeometryOptions
{
  /// divisions and tolerance of the containing-element search
  SearchOptions search;
  /// the tag holding global entity ids
  std::string global_id_name = "global";
  /// per-vertex ownership override: a vertex the mask is 0 for is not owned,
  /// whatever the mesh says; empty means the mesh's own ownership. Layouts
  /// whose DOF holders are the mesh vertices apply it; layouts of other
  /// entities do not.
  Omega_h::Read<Omega_h::I8> owned_mask;
};

/// The geometry of an Omega_h simplex mesh (2D or 3D, classified). Builds the
/// mesh's containing-element search: create one per mesh and share it. The
/// geometry owns the mesh's discretization; spaces are created on it with
/// CreateLagrangeSpace and the other space factories.
/// @param mesh the mesh; must outlive the geometry and everything created
///        from it
/// @param system the mesh's coordinate system, or null for Cartesian
/// @param options search, global id tag and ownership
[[nodiscard]] std::shared_ptr<const Geometry> CreateGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options = {});

/// The mesh of an Omega_h geometry.
/// @throws pcms_error if `geometry` is not an Omega_h geometry
[[nodiscard]] Omega_h::Mesh& GetMesh(const Geometry& geometry);

/// A point-cloud geometry whose sites are the vertices (`entity_dim` 0) or the
/// element centroids (the mesh dimension) of the mesh of an Omega_h geometry;
/// the support search of its neighbourhoods is the mesh's adjacency search.
/// Spaces made on the same geometry object share its entities, so binding the
/// cloud at their DOF holders starts from those entities. For
/// CreatePolynomialReconstructionSpace.
/// @param geometry an Omega_h geometry
/// @param entity_dim 0 or the mesh dimension
/// @param options the support search of the geometry
/// @throws pcms_error if `geometry` is not an Omega_h geometry or `entity_dim`
///         is not 0 or the mesh dimension
[[nodiscard]] std::shared_ptr<const Geometry> CreatePointCloudGeometry(
  std::shared_ptr<const Geometry> geometry, int entity_dim,
  NeighbourhoodOptions options = {});

/// 0/1 field on the vertices of a 2D uniform grid marking the vertices that
/// lie inside mesh. The grid is given, or covers the mesh's bounding box with
/// the given number of cells per axis.
Field<Real> CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                                         const UniformGrid<2>& grid);
Field<Real> CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                                         std::array<LO, 2> divisions);

} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
