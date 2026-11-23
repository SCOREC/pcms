#include "pcms/backends/omega_h/point_search.hpp"
#include "pcms/localization/overlap.hpp"
#include <Omega_h_mesh.hpp>
#include <Omega_h_simplex.hpp>
#include <stdexcept>

namespace pcms::omega_h
{

/**
 * Check if a triangle element represented by 3 coordinates in two dimensions
 * intersects with a bounding box.
 */
[[nodiscard]] KOKKOS_FUNCTION bool TriangleIntersectsBBox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox)
{
  return Overlaps(detail::ToSimplex<2>(coords), bbox);
}

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

template <int Dim>
GridPointSearch<Dim> CreateGridPointSearchOfDim(
  Omega_h::Mesh& mesh, const std::vector<LO>& divisions, Real tolerance)
{
  auto mesh_bbox = Omega_h::get_bounding_box<Dim>(&mesh);
  UniformGrid<Dim> grid{};
  for (int i = 0; i < Dim; ++i) {
    grid.edge_length[i] = mesh_bbox.max[i] - mesh_bbox.min[i];
    grid.bot_left[i] = mesh_bbox.min[i];
    grid.divisions[i] = divisions[i];
  }
  Kokkos::Array<Kokkos::View<const LO*>, Dim> down{};
  Kokkos::Array<Kokkos::View<const LO*>, Dim> down_verts{};
  for (int d = 1; d < Dim; ++d) {
    down[d] = mesh.ask_down(Dim, d).ab2b.view();
    down_verts[d] = mesh.ask_down(d, Omega_h::VERT).ab2b.view();
  }
  const auto dual = mesh.ask_dual();
  return GridPointSearch<Dim>(mesh.coords().view(),
                              mesh.ask_elem_verts().view(), down, down_verts,
                              EntityVertexMasks<Dim>(), dual.a2ab.view(),
                              dual.ab2b.view(), grid, tolerance);
}
} // namespace

GridPointSearchVariant CreateGridPointSearch(Omega_h::Mesh& mesh,
                                             const std::vector<LO>& divisions,
                                             Real tolerance)
{
  if (mesh.dim() != 2 && mesh.dim() != 3) {
    throw std::invalid_argument(
      "CreateGridPointSearch: only 2D and 3D meshes are supported");
  }
  if (divisions.size() != static_cast<std::size_t>(mesh.dim())) {
    throw std::invalid_argument(
      "CreateGridPointSearch: divisions must have one entry per axis");
  }
  if (mesh.dim() == 2) {
    return CreateGridPointSearchOfDim<2>(mesh, divisions, tolerance);
  }
  return CreateGridPointSearchOfDim<3>(mesh, divisions, tolerance);
}

} // namespace pcms::omega_h
