#ifndef PCMS_COUPLING_POINT_SEARCH_H
#define PCMS_COUPLING_POINT_SEARCH_H
#include <cassert>

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>

#include "types.h"
#include "pcms/uniform_grid.h"
#include "pcms/bounding_box.h"

namespace pcms
{
//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
namespace detail {
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh, Kokkos::View<Uniform2DGrid[1]> grid, int num_grid_cells);
}
KOKKOS_FUNCTION

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

template <int dim>
class PointLocalizationSearch
{
public:
  struct Result {
    enum class Dimensionality
    {
      VERTEX = 0,
      EDGE = 1,
      FACE = 2,
      REGION = 3
    };

    Dimensionality dimensionality;
    LO element_id;
    Omega_h::Vector<dim + 1> parametric_coords;
  };

  static constexpr auto DIM = dim;

  using PointSearchTolerances = Kokkos::View<Real[DIM]>;

  explicit PointLocalizationSearch(const PointSearchTolerances& tolerances) : tolerances_(tolerances)
  {
    assert(tolerances_.is_allocated());
  }

  virtual Kokkos::View<Result*> operator()(Kokkos::View<Real*[dim] > point) const = 0;
  virtual ~PointLocalizationSearch() = default;

protected:
  PointSearchTolerances tolerances_;
};

using PointLocalizationSearch2D = PointLocalizationSearch<2>;
using PointLocalizationSearch3D = PointLocalizationSearch<3>;

class GridPointSearch2D : public PointLocalizationSearch2D
{
  using CandidateMapT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch2D::Result;

  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny);
  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny, const PointSearchTolerances& tolerances);

  /**
   *  given a point in global coordinates give the id of the triangle that the
   * point lies within and the parametric coordinate of the point within the
   * triangle. If the point does not lie within any triangle element. Then the
   * id will be a negative number and (TODO) will return a negative id of the
   * closest element
   */
  Kokkos::View<Result*> operator()(Kokkos::View<Real* [DIM]> point) const override;

private:
  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Kokkos::View<Uniform2DGrid[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

class GridPointSearch3D : public PointLocalizationSearch3D
{
  using CandidateMapT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch3D::Result;

  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz);
  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz, const PointSearchTolerances& tolerances);

  /**
   *  given a point in global coordinates give the id of the triangle that the
   * point lies within and the parametric coordinate of the point within the
   * triangle. If the point does not lie within any triangle element. Then the
   * id will be a negative number and (TODO) will return a negative id of the
   * closest element
   */
  Kokkos::View<Result*> operator()(Kokkos::View<Real* [DIM]> point) const override;

private:
  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Kokkos::View<UniformGrid<DIM>[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

} // namespace detail
#endif // PCMS_COUPLING_POINT_SEARCH_H
