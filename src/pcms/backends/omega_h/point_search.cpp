#include "pcms/backends/omega_h/point_search.hpp"
#include "pcms/localization/overlap.hpp"
#include <Omega_h_mesh.hpp>
#include <Omega_h_simplex.hpp>

namespace pcms::omega_h
{

/**
 * Check if a triangle element represented by 3 coordinates in two dimensions
 * intersects with a bounding box.
 */
[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox)
{
  return overlaps(detail::ToSimplex<2>(coords), bbox);
}

namespace detail
{

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_2d(Omega_h::Mesh& mesh,
                                        Kokkos::View<UniformGrid2D[1]> grid,
                                        int num_grid_cells)
{
  return localization::detail::construct_intersection_map_reference_2d(
    mesh.coords().view(), mesh.ask_elem_verts().view(), grid, num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_3d(Omega_h::Mesh& mesh,
                                        Kokkos::View<UniformGrid3D[1]> grid,
                                        int num_grid_cells)
{
  return localization::detail::construct_intersection_map_reference_3d(
    mesh.coords().view(), mesh.ask_elem_verts().view(), grid, num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<UniformGrid2D[1]> grid,
                              int num_grid_cells, Real inflation)
{
  return localization::detail::construct_intersection_map_2d(
    mesh.coords().view(), mesh.ask_elem_verts().view(), grid, num_grid_cells,
    inflation);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Omega_h::Mesh& mesh,
                              Kokkos::View<UniformGrid3D[1]> grid,
                              int num_grid_cells, Real inflation)
{
  return localization::detail::construct_intersection_map_3d(
    mesh.coords().view(), mesh.ask_elem_verts().view(), grid, num_grid_cells,
    inflation);
}

} // namespace detail

namespace
{
// Per entity dimension 1..Dim-1 and local entity of an element, its local
// vertices as a bitmask, in Omega_h's order of the element's entities.
template <int Dim>
Kokkos::Array<Kokkos::Array<unsigned, 6>, Dim> EntityVertexMasks()
{
  Kokkos::Array<Kokkos::Array<unsigned, 6>, Dim> masks{};
  for (int dim = 1; dim < Dim; ++dim) {
    for (int j = 0; j < Omega_h::simplex_degree(Dim, dim); ++j) {
      for (int v = 0; v <= dim; ++v) {
        masks[dim][j] |= 1u << Omega_h::simplex_down_template(Dim, dim, j, v);
      }
    }
  }
  return masks;
}
} // namespace

GridPointSearch2D MakeGridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                        Real tolerance)
{
  auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
  const UniformGrid2D grid{.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                                           mesh_bbox.max[1] - mesh_bbox.min[1]},
                           .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
                           .divisions = {Nx, Ny}};
  Kokkos::Array<Kokkos::View<const LO*>, 2> down{};
  Kokkos::Array<Kokkos::View<const LO*>, 2> down_verts{};
  down[1] = mesh.ask_down(Omega_h::FACE, Omega_h::EDGE).ab2b.view();
  down_verts[1] = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT).ab2b.view();
  const auto dual = mesh.ask_dual();
  return GridPointSearch2D(mesh.coords().view(), mesh.ask_elem_verts().view(),
                           down, down_verts, EntityVertexMasks<2>(),
                           dual.a2ab.view(), dual.ab2b.view(), grid, tolerance);
}

GridPointSearch3D MakeGridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                        LO Nz, Real tolerance)
{
  auto mesh_bbox = Omega_h::get_bounding_box<3>(&mesh);
  Kokkos::Array<Real, 3> edge_lengths{};
  Kokkos::Array<Real, 3> bot_left{};
  for (int i = 0; i < 3; ++i) {
    edge_lengths[i] = mesh_bbox.max[i] - mesh_bbox.min[i];
    bot_left[i] = mesh_bbox.min[i];
  }
  const UniformGrid3D grid{.edge_length = edge_lengths,
                           .bot_left = bot_left,
                           .divisions = {Nx, Ny, Nz}};
  Kokkos::Array<Kokkos::View<const LO*>, 3> down{};
  Kokkos::Array<Kokkos::View<const LO*>, 3> down_verts{};
  down[1] = mesh.ask_down(Omega_h::REGION, Omega_h::EDGE).ab2b.view();
  down[2] = mesh.ask_down(Omega_h::REGION, Omega_h::FACE).ab2b.view();
  down_verts[1] = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT).ab2b.view();
  down_verts[2] = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b.view();
  const auto dual = mesh.ask_dual();
  return GridPointSearch3D(mesh.coords().view(), mesh.ask_elem_verts().view(),
                           down, down_verts, EntityVertexMasks<3>(),
                           dual.a2ab.view(), dual.ab2b.view(), grid, tolerance);
}

} // namespace pcms::omega_h
