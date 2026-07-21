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

//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
namespace detail
{
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform2DGrid[1]> grid,
                              int num_grid_cells);
}

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

template <int dim>
class PointLocalizationSearch
{
public:
  struct Result
  {
    enum class Dimensionality
    {
      VERTEX = 0,
      EDGE = 1,
      FACE = 2,
      REGION = 3
    };

    Dimensionality dimensionality;
    LO element_id; // ID of the actual entity (vertex/edge/face/region)
    Omega_h::Vector<dim + 1> parametric_coords;
  };

  static constexpr auto DIM = dim;

  using PointSearchTolerances = Kokkos::View<Real[DIM]>;

  explicit PointLocalizationSearch(const PointSearchTolerances& tolerances)
    : tolerances_(tolerances)
  {
    assert(tolerances_.is_allocated());
  }

  virtual Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [dim]> point) const = 0;

  /**
   * This function provides a temporary solution to retrieve the original
   * behavior of the point search, which previously returned the face ID
   * regardless of which entity the search result belonged to. Many parts of
   * the codebase still rely on this legacy behavior. After updating the point
   * search to return the exact entity, this function can be called to retrieve
   * the old behavior and ensure correctness. Long term, the implementation
   * should be updated to properly handle the new result.
   */
  [[nodiscard]] virtual LO GetOwningElementId(const Result& result) = 0;
  [[nodiscard]] virtual Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) = 0;
  virtual ~PointLocalizationSearch() = default;

protected:
  PointSearchTolerances tolerances_;
};

using PointLocalizationSearch2D = PointLocalizationSearch<2>;
using PointLocalizationSearch3D = PointLocalizationSearch<3>;

class GridPointSearch2D : public PointLocalizationSearch2D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch2D::Result;

  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny);
  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                    const PointSearchTolerances& tolerances);

  /**
   *  given a point in global coordinates give the id of the triangle that the
   * point lies within and the parametric coordinate of the point within the
   * triangle. If the point does not lie within any triangle element. Then the
   * id will be a negative number and (TODO) will return a negative id of the
   * closest element
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;
  [[nodiscard]] LO GetOwningElementId(const Result& result) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) override;

private:
  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj edges2faces_up_;
  Omega_h::Adj verts2faces_up_;
  Kokkos::View<Uniform2DGrid[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

class GridPointSearch3D : public PointLocalizationSearch3D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch3D::Result;

  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz);
  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                    const PointSearchTolerances& tolerances);

  /**
   *  Given a point in global coordinates, returns the id of the tetrahedron (3D
   * element) that the point lies within and the parametric coordinate of the
   * point within the tetrahedron. If the point does not lie within any
   * tetrahedron element, then the id will be a negative number and (TODO) will
   * return a negative id of the closest element.
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;
  [[nodiscard]] LO GetOwningElementId(const Result& result) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) override;

private:
  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj verts2regions_up_;
  Omega_h::Adj edges2regions_up_;
  Omega_h::Adj faces2regions_up_;
  Kokkos::View<UniformGrid<DIM>[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
  Real fuzz_;
};

} // namespace pcms
#endif // PCMS_COUPLING_POINT_SEARCH_H
