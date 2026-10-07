#ifndef PCMS_COUPLING_POINT_SEARCH_H
#define PCMS_COUPLING_POINT_SEARCH_H
#include <cassert>

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>

#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"
#include "pcms/utility/bounding_box.h"
#include "pcms/field/coordinate_system.h"

namespace pcms
{

/**
 * Get the owning element (face for 2D, region for 3D) for a given entity.
 *
 * Given an entity of dimension `entity_dim` with ID `element_id`, returns the
 * smallest ID of the owning element of dimension `mesh_dim` (the mesh's
 * top-level dimension).
 */
LO GetOwningElementId(Omega_h::Mesh& mesh, int mesh_dim, int entity_dim,
                      LO element_id);

template <typename DimEnum>
KOKKOS_INLINE_FUNCTION LO
GetOwningElementIdFromAdj(const Omega_h::Adj& upward_adj, DimEnum entity_dim,
                          DimEnum target_dim, LO element_id)
{
  if (element_id < 0)
    return -1;

  // If entity is already at target dimension, return it directly
  if (static_cast<int>(entity_dim) == static_cast<int>(target_dim))
    return element_id;

  const auto begin = upward_adj.a2ab[element_id];
  const auto end = upward_adj.a2ab[element_id + 1];
  if (begin >= end)
    return -1;

  // Find the smallest owning element ID
  LO owner = upward_adj.ab2b[begin];
  for (auto i = begin + 1; i < end; ++i) {
    const LO candidate = upward_adj.ab2b[i];
    if (candidate < owner)
      owner = candidate;
  }
  return owner;
}
/**
 * Point search base class
 */
class PointSearch
{
public:
  using ExecSpace = Omega_h::ExecSpace;
  using MemorySpace = ExecSpace::memory_space;
  /**
   * Results type gives dimensionalities of point intersection, the intersected
   * element IDs, and the barycentric coordinate mappings of each point as
   * follows: the `i`th input point has intersection dimensionality of
   * `dimensionalities(i)`, intersected element ID of `element_ids(i)`, and
   * likewise parametric coordinate mapping of `parametric_coords(i,0..d)`
   * where `d` is the dimension of the space
   */
  struct Results
  {
    enum class Dimensionality
    {
      VERTEX = 0,
      EDGE = 1,
      FACE = 2,
      REGION = 3,
      NO_INTERSECT = 4
    };

    Kokkos::View<Dimensionality*, MemorySpace> dimensionalities;
    Kokkos::View<LO*, MemorySpace> element_ids;
    Kokkos::View<Real**, MemorySpace> parametric_coords;
  };

  using PointSearchTolerances = Kokkos::View<Real*>;

  explicit PointSearch(const PointSearchTolerances& tolerances)
    : tolerances_(tolerances)
  {
    assert(tolerances_.is_allocated());
  }

  ~PointSearch() = default;

  virtual Results Apply(const CoordinateView<MemorySpace>& coords) const = 0;
  /**
   * This function provides a temporary solution to retrieve the original
   * behavior of the point search, which previously returned the face ID
   * regardless of which entity the search result belonged to. Many parts of
   * the codebase still rely on this legacy behavior. After updating the point
   * search to return the exact entity, this function can be called to retrieve
   * the old behavior and ensure correctness. Long term, the implementation
   * should be updated to properly handle the new result.
   */
  [[nodiscard]] virtual LO GetOwningElementId(const Results& results,
                                              int i) = 0;
  [[nodiscard]] virtual Kokkos::View<LO*> GetOwningElementIds(
    const Results& results) = 0;

protected:
  PointSearchTolerances tolerances_;
};

} // namespace pcms
#endif // PCMS_COUPLING_POINT_SEARCH_H
