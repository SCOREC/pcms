#ifndef PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
#define PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP

// Public entry points of the Omega_h backend: build PCMS objects from an
// Omega_h mesh. This header (and only backend headers) may use Omega_h types.

#include "pcms/discretization/discretization.hpp"
#include "pcms/field/field.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_mesh.hpp>

#include <array>
#include <memory>
#include <string>

namespace pcms::backend::omega_h
{

/// Discretization on an Omega_h mesh. The mesh is borrowed and must outlive
/// the discretization and everything created from it. global_id_name names
/// the tag holding global entity ids.
std::shared_ptr<const Discretization> CreateDiscretization(
  Omega_h::Mesh& mesh, std::string global_id_name = "global");

/// As above, with per-vertex ownership given by owned_mask (1 = owned)
/// instead of the mesh's partition.
std::shared_ptr<const Discretization> CreateDiscretization(
  Omega_h::Mesh& mesh, Omega_h::Read<Omega_h::I8> owned_mask,
  std::string global_id_name = "global");

/// The Omega_h mesh of a discretization created by CreateDiscretization.
/// Throws pcms_error if the discretization is not Omega_h-based.
Omega_h::Mesh& GetMesh(const Discretization& discretization);

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

/// Uniform grid covering the bounding box of mesh with divisions[i] cells
/// along axis i.
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

/// Uniform grid covering the bounding box of mesh with cells_per_dim cells
/// along every axis.
template <unsigned dim = 2>
UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,
                                           LO cells_per_dim)
{
  std::array<LO, dim> divisions;
  for (unsigned i = 0; i < dim; ++i)
    divisions[i] = cells_per_dim;
  return CreateUniformGridFromMesh<dim>(mesh, divisions);
}

/// Host mask (true where mask[i] != 0) for OverlapMask(size, mask).
inline Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromArray(
  Omega_h::HostRead<Omega_h::I8> mask)
{
  Kokkos::View<bool*, HostMemorySpace> out("overlap_info", mask.size());
  for (Omega_h::LO i = 0; i < mask.size(); ++i) {
    out(i) = static_cast<bool>(mask[i]);
  }
  return out;
}

inline Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromArray(
  Omega_h::Read<Omega_h::I8> mask)
{
  return OverlapMaskFromArray(Omega_h::HostRead<Omega_h::I8>(mask));
}

} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
