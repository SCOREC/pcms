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
//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh, const UniformGrid& grid);

Omega_h::Vector<3> barycentric_from_global(
  const Omega_h::Vector<2>& point, const Omega_h::Matrix<2, 3>& vertex_coords);


[[nodiscard]]
KOKKOS_FUNCTION
bool triangle_intersects_bbox(const Omega_h::Matrix<2, 3>& coords,
                              const AABBox<2>& bbox);

// PUBLIC INTERFACE BELOW...
// functors needed
// get parametric inversion. Given set of candidate elements do parametric
// inversion
//

// given array of coordinate data, find the associated element that that point
// lies within (or closest element) and the parametric coordinate
class GridPointSearch
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  static constexpr auto dim = 2;
  // LO closest_cell_id_from_point(const BBOX&, const n_divisions&, const
public:
  GridPointSearch(Omega_h::Mesh& mesh, LO Nx, LO Ny)

  {
    auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
    // get mesh bounding box
    grid_ = {.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                             mesh_bbox.max[1] - mesh_bbox.min[1]},
             .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
             .divisions = {Nx, Ny}};
    candidate_map_ = construct_intersection_map(mesh, grid_);
    coords_ = mesh.coords();
    tris2verts_ = mesh.ask_elem_verts();
    // TODO intiaialize coords and tri2vert and mesh!
    // TODO pass coords& tri2vert to construct map by references!
  }
  // TODO use a functor for parallel for over a list of points
  // TODO return result struct rather than pair
  std::pair<LO, Omega_h::Vector<dim + 1>> operator()(Omega_h::Vector<dim> point)
  {
    auto cell_id = grid_.ClosestCellID(point);
    assert(cell_id < candidate_map_.numRows() && cell_id >= 0);
    auto candidates_begin = candidate_map_.row_map(cell_id);
    auto candidates_end = candidate_map_.row_map(cell_id + 1);
    // create array that's size of number of candidates x num coords to store
    // parametric inversion
    for (auto i = candidates_begin; i < candidates_end; ++i) {
      const auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts_, candidate_map_.entries(i));
      // 2d mesh with 2d coords, but 3 triangles
      const auto vertex_coords = Omega_h::gather_vectors<3, 2>(coords_, elem_tri2verts);
      auto parametric_coords = barycentric_from_global(point, vertex_coords);
      if (Omega_h::is_barycentric_inside(parametric_coords)) {
        return {candidate_map_.entries(i), parametric_coords};
      }
    }
    return {-1, {}};
  }

private:
  Omega_h::Mesh mesh_;
  UniformGrid grid_;
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

} // namespace wdmcpl

#endif // WDM_COUPLING_POINT_SEARCH_H
