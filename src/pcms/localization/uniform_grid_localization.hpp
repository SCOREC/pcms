#ifndef PCMS_COUPLING_UNIFORM_GRID_LOCALIZATION_HPP
#define PCMS_COUPLING_UNIFORM_GRID_LOCALIZATION_HPP

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>

#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"
#include "pcms/utility/bounding_box.h"
#include "pcms/field/coordinate_system.h"
#include "pcms/localization/localization.h"

namespace pcms
{
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

class GridPointSearch2D : public PointSearch
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  static constexpr int DIM = 2;

  using Results = PointSearch::Results;
  using Dimensionality = Results::Dimensionality;

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
  Results Apply(const CoordinateView<MemorySpace>& coords) const override;
  [[nodiscard]] LO GetOwningElementId(const Results& result, int i) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    const Results& results) override;

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

class GridPointSearch3D : public PointSearch
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  static constexpr int DIM = 3;

  using Results = PointSearch::Results;
  using Dimensionality = Results::Dimensionality;

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
  Results Apply(const CoordinateView<MemorySpace>& coords) const override;
  [[nodiscard]] LO GetOwningElementId(const Results& result, int i) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    const Results& results) override;

private:
  Omega_h::Mesh mesh_;
  Omega_h::Adj tets2faces_adj_;
  Omega_h::Adj tets2edges_adj_;
  Omega_h::Adj tets2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj verts2regions_up_;
  Omega_h::Adj edges2regions_up_;
  Omega_h::Adj faces2regions_up_;
  Kokkos::View<UniformGrid<DIM>[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tets2verts_;
  Omega_h::Reals coords_;
  Real fuzz_;
};
} // namespace pcms
#endif
