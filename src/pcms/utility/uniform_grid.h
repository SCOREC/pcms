#ifndef PCMS_COUPLING_UNIFORM_GRID_H
#define PCMS_COUPLING_UNIFORM_GRID_H
#include "pcms/utility/bounding_box.h"
#include "Omega_h_comm.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_bbox.hpp"
#include "Omega_h_mesh.hpp"
#include <Kokkos_Array.hpp>
#include <algorithm>
#include <array>
#include <cmath>
namespace pcms
{

template <unsigned dim>
struct UniformGrid
{
  // Make private?
  Kokkos::Array<Real, dim> edge_length;
  Kokkos::Array<Real, dim> bot_left;
  Kokkos::Array<LO, dim> divisions;

public:
  [[nodiscard]] LO GetNumCells() const
  {
    LO total = 1;
    for (std::size_t i = 0; i < dim; ++i)
      total *= divisions[i];
    return total;
  }
  /// return the grid cell ID that the input point is inside or closest to if
  /// the point lies outside
  // take the view as a template because it might be a subview type
  // template <typename T>
  //[[nodiscard]] KOKKOS_INLINE_FUNCTION LO ClosestCellID(const T& point) const
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  ClosestCellID(const Omega_h::Vector<dim>& point) const
  {
    Kokkos::Array<LO, dim> indexes;

    for (size_t i = 0; i < dim; ++i) {
      indexes[i] = AxisCellIndex(i, point[i]);
    }
    return CellIndexFromAxisIndices(indexes);
  }

  /// Index of the cell containing `coord` along `axis`, clamped to the grid.
  /// The single definition of the grid's binning rule: keeping the point query
  /// and the cell-range query on this one formula is what guarantees a point
  /// always lands in a cell that its own element was binned into.
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO AxisCellIndex(size_t axis,
                                                        Real coord) const
  {
    const auto index = static_cast<LO>(std::floor(
      (coord - bot_left[axis]) * divisions[axis] / edge_length[axis]));
    return Kokkos::clamp(index, 0, divisions[axis] - 1);
  }

  /// Linear cell ID from per-axis indices given in coordinate order (x,y,z).
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  CellIndexFromAxisIndices(Kokkos::Array<LO, dim> indexes) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(indexes);
    return GetCellIndex(indexes);
  }

  /// Per-axis cell indices in coordinate order (x,y,z) for a linear cell ID.
  ///
  /// The exact inverse of `CellIndexFromAxisIndices`. Provided so callers that
  /// work with linear cell IDs -- a routing table keyed by cell, for example --
  /// do not have to reimplement the grid's internal index ordering, which is
  /// reversed relative to coordinates and easy to get wrong silently.
  [[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<LO, dim>
  AxisIndicesFromCellIndex(LO idx) const
  {
    auto indices = GetDimensionedIndex(idx);
    reverse(indices);
    return indices;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<dim> GetCellBBOX(LO idx) const
  {
    auto index = GetDimensionedIndex(idx);
    reverse(index);

    Kokkos::Array<Real, dim> half_width, center;

    for (size_t i = 0; i < dim; ++i) {
      half_width[i] = edge_length[i] / divisions[i] / 2;
    }

    for (size_t i = 0; i < dim; ++i) {
      center[i] = (2.0 * index[i] + 1.0) * half_width[i] + bot_left[i];
    }

    return {.center = center, .half_width = half_width};
  }

  /// Check if a point is within the grid bounds
  [[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsPointInBounds(
    const Omega_h::Vector<dim>& point) const
  {
    for (size_t i = 0; i < dim; ++i) {
      Real coord = point[i];
      Real grid_min = bot_left[i];
      Real grid_max = bot_left[i] + edge_length[i];
      if (coord < grid_min || coord > grid_max) {
        return false;
      }
    }
    return true;
  }

  /// As `IsPointInBounds`, but a point within `tolerance` of the box on every
  /// side still counts as inside. Used to reject query points that no element
  /// can contain: the tolerance matches the band where the search still accepts
  /// a point just outside a boundary face.
  [[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsPointInBounds(
    const Omega_h::Vector<dim>& point, Real tolerance) const
  {
    for (size_t i = 0; i < dim; ++i) {
      Real coord = point[i];
      Real grid_min = bot_left[i] - tolerance;
      Real grid_max = bot_left[i] + edge_length[i] + tolerance;
      if (coord < grid_min || coord > grid_max) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<LO, dim>
  GetDimensionedIndex(LO idx) const
  {
    LO stride = 1;
    for (std::size_t i = 0; i < dim - 1; ++i) {
      stride *= divisions[i];
    }
    Kokkos::Array<LO, dim> result;

    for (size_t i = 0; i < dim; ++i) {
      result[i] = idx / stride;
      idx -= result[i] * stride;
      stride /= divisions[i];
    }

    return result;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  GetCellIndex(Kokkos::Array<LO, dim> dimensionedIndex) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(dimensionedIndex);

    LO idx = 0;
    LO stride = 1;

    for (size_t i = 0; i < dim; ++i) {
      idx += dimensionedIndex[i] * stride;
      stride *= divisions[i];
    }

    return idx;
  }

private:
  template <typename T, std::size_t N>
  KOKKOS_INLINE_FUNCTION static void reverse(Kokkos::Array<T, N>& arr)
  {
    for (size_t i = 0, j = N - 1; i < j; ++i, --j) {
      auto temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }
};

using Uniform2DGrid = UniformGrid<2>;
using Uniform3DGrid = UniformGrid<3>;

/**
 * \brief Divisions per axis giving roughly one element per grid cell.
 *
 */
template <unsigned dim>
LO DivisionsPerAxisForMesh(LO nelems, LO min_divisions = 10)
{
  const auto per_axis = static_cast<LO>(std::ceil(
    std::pow(static_cast<double>(nelems), 1.0 / static_cast<double>(dim))));
  return std::max(per_axis, min_divisions);
}

/**
 * \brief Number of elements this rank owns, counted on the device.
 *
 * Used instead of `Omega_h::Mesh::nents_owned`, which is unusable in a
 * Kokkos+CUDA build without OpenMP: its reduction functor reads a device view
 * from a host loop (`Omega_h_reduce.hpp`, the non-OpenMP `transform_reduce`),
 * which trips Kokkos' memory-space access check.
 *
 * \param mesh A mesh that carries ownership information, i.e. one that has been
 * partitioned. Callers must establish that before calling, typically from
 * `nghost_layers() > 0`.
 * \return The number of elements this rank owns
 */
inline LO OwnedElementCount(Omega_h::Mesh& mesh)
{
  const auto ent_dim = mesh.dim();
  auto owned = mesh.owned(ent_dim);
  LO count = 0;
  Kokkos::parallel_reduce(
    "pcms::OwnedElementCount", owned.size(),
    KOKKOS_LAMBDA(LO i, LO & update) { update += owned[i] ? 1 : 0; }, count);
  return count;
}

/**
 * \brief Number of elements in the whole decomposition.
 *
 * This is the quantity a grid resolution must be derived from so that every
 * rank constructs the *same* grid. `mesh.nelems()` cannot be used for that
 * purpose on its own: it counts the ghost copies a rank holds as well as its
 * owned elements, so it varies with the partition. Summing it does not vary,
 * because every rank reduces the same total, which is the property that
 * matters.
 *
 * The count is exact. On a mesh with no ghost layers each element lives on
 * exactly one rank, so the reduced local count is already the global count. On
 * a ghosted mesh the ghost copies make that local count an overestimate, so the
 * owned elements are counted instead -- ghosting implies the mesh was
 * partitioned and therefore carries the ownership information that needs.
 *
 * \param mesh The mesh to count elements for
 * \return The global number of elements, identical on every rank
 */
inline LO GlobalElementCount(Omega_h::Mesh& mesh)
{
  const auto local_count =
    mesh.nghost_layers() > 0 ? OwnedElementCount(mesh) : mesh.nelems();
  return mesh.comm()->allreduce(local_count, Omega_h_Op::OMEGA_H_SUM);
}

/**
 * \brief Divisions per axis giving roughly one element per grid cell, derived
 * from the *global* element count.
 *
 * The distributed counterpart of `DivisionsPerAxisForMesh(mesh.nelems())`,
 * which silently produces a different grid on every rank because `nelems()`
 * includes ghosts. Reach for this function whenever the resulting grid has to
 * describe the whole decomposition rather than one partition.
 *
 * \tparam dim Spatial dimension of the grid
 * \param mesh The mesh the grid is built over
 * \param min_divisions Lower bound on the number of cells per axis
 * \return The same value on every rank
 */
template <unsigned dim>
LO DivisionsPerAxisForGlobalMesh(Omega_h::Mesh& mesh, LO min_divisions = 10)
{
  return DivisionsPerAxisForMesh<dim>(GlobalElementCount(mesh), min_divisions);
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
  // Get the bounding box of the mesh
  auto bbox = Omega_h::get_bounding_box<dim>(&mesh);

  // Calculate edge lengths and bottom-left corner
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

} // namespace pcms

#endif // PCMS_COUPLING_UNIFORM_GRID_H
