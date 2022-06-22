#ifndef WDM_COUPLING_POINT_SEARCH_H
#define WDM_COUPLING_POINT_SEARCH_H
#include <unordered_map>
#include <Kokkos_Crs.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include "types.h"
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>
#include "wdmcpl/uniform_grid.h"

namespace wdmcpl
{

namespace detail
{
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh, const UniformGrid& grid);
} // namespace detail
Omega_h::Vector<3> barycentric_from_global(
  const Omega_h::Vector<2>& point, const Omega_h::Matrix<2, 3>& vertex_coords);

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

class GridPointSearch
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  static constexpr auto dim = 2;

public:
  GridPointSearch(Omega_h::Mesh& mesh, LO Nx, LO Ny);
  /**
   *  given a point in global coordinates give the id of the triangle that the
   * point lies within and the parametric coordinate of the point within the
   * triangle. If the point does not lie within any triangle element. Then the
   * id will be a negative number and (TODO) will return a negative id of the
   * closest element
   */
  std::pair<LO, Omega_h::Vector<dim + 1>> operator()(
    Omega_h::Vector<dim> point);

private:
  Omega_h::Mesh mesh_;
  UniformGrid grid_{};
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

} // namespace wdmcpl

#endif // WDM_COUPLING_POINT_SEARCH_H
