#include "point_search.h"
#include <Omega_h_mesh.hpp>
#include <bitset>

namespace wdmcpl
{
KOKKOS_INLINE_FUNCTION
AABBox<2> triangle_bbox(const Omega_h::Matrix<2, 3>& coords)
{
  std::array<Real, 2> max{coords(0, 0), coords(1, 0)};
  std::array<Real, 2> min{coords(0, 0), coords(1, 0)};
  for (int i = 1; i < 3; ++i) {
    max[0] = Kokkos::fmax(max[0], coords(0, i));
    max[1] = Kokkos::fmax(max[1], coords(1, i));
    min[0] = Kokkos::fmin(min[0], coords(0, i));
    min[1] = Kokkos::fmin(min[1], coords(1, i));
  }
  return {.center = {(max[0] + min[0]) / 2.0, (max[1] + min[1]) / 2.0},
          .half_width = {(max[0] - min[0]) / 2.0, (max[1] - min[1]) / 2.0}};
}
// Liang, You-Dong, and B. A. Barsky. “A New Concept and Method for Line
// Clipping.” ACM Transactions on Graphics 3, no. 1 (January 1984): 1–22.
// https://doi.org/10.1145/357332.357333.
KOKKOS_INLINE_FUNCTION
bool clipt(Real p, Real q, Real& t0, Real& t1)
{
  if (p < 0) {
    auto r = q / p;
    if (r > t1) {
      return false;
    }
    if (r > t0) {
      t0 = r;
    }
  } else if (p > 0) {
    auto r = q / p;
    if (r < t0) {
      return false;
    }
    if (r < t1) {
      t1 = r;
    }
  }
  // p==0
  else {
    if (q < 0) {
      return false;
    }
  }
  return true;
}
// Liang, You-Dong, and B. A. Barsky. “A New Concept and Method for Line
// Clipping.” ACM Transactions on Graphics 3, no. 1 (January 1984): 1–22.
// https://doi.org/10.1145/357332.357333.
KOKKOS_INLINE_FUNCTION
bool line_intersects_bbox(const Omega_h::Vector<2>& p0,
                          const Omega_h::Vector<2>& p1, const AABBox<2>& bbox)
{
  Real t0 = 0;
  Real t1 = 1;
  auto xleft = bbox.center[0] - bbox.half_width[0];
  auto deltax = p1[0] - p0[0];
  if (clipt(-deltax, p0[0] - xleft, t0, t1)) {
    auto xright = bbox.center[0] + bbox.half_width[0];
    if (clipt(deltax, xright - p0[0], t0, t1)) {
      auto deltay = p1[1] - p0[1];
      auto ybottom = bbox.center[1] - bbox.half_width[1];
      if (clipt(-deltay, p0[1] - ybottom, t0, t1)) {
        auto ytop = bbox.center[1] + bbox.half_width[1];
        if (clipt(deltay, ytop - p0[1], t0, t1)) {
          // full liang-barksy algorithm computes new x-y coordinates here, but
          // we just need to check intersections
          return true;
        }
      }
    }
  }
  return false;
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION
bool within_bbox(const Omega_h::Vector<2> coord, const AABBox<2> & bbox) noexcept {
  auto left = bbox.center[0] - bbox.half_width[0];
  auto right = bbox.center[0] + bbox.half_width[0];
  auto bot = bbox.center[1] - bbox.half_width[1];
  auto top = bbox.center[1] + bbox.half_width[1];
  return (coord[0]>=left) && (coord[0]<=right) && (coord[1] >= bot) && (coord[1]<=top);
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION
bool bbox_verts_within_triangle(const AABBox<2>& bbox, const Omega_h::Matrix<2,3>& coords) {
  auto left = bbox.center[0] - bbox.half_width[0];
  auto right = bbox.center[0] + bbox.half_width[0];
  auto bot = bbox.center[1] - bbox.half_width[1];
  auto top = bbox.center[1] + bbox.half_width[1];
  auto xi = barycentric_from_global({left,bot}, coords);
  if(Omega_h::is_barycentric_inside(xi)){ return true;}
  xi = barycentric_from_global({left,top}, coords);
  if(Omega_h::is_barycentric_inside(xi)){ return true;}
  xi = barycentric_from_global({right,top}, coords);
  if(Omega_h::is_barycentric_inside(xi)){ return true;}
  xi = barycentric_from_global({right,bot}, coords);
  if(Omega_h::is_barycentric_inside(xi)){ return true;}
  return false;
}

/**
 * Check if a triangle element represented by 3 coordinates in two dimensions
 * intersects with a bounding box
 */
[[nodiscard]] bool triangle_intersects_bbox(const Omega_h::Matrix<2, 3>& coords,
                                            const AABBox<2>& bbox)
{
  // triangle and grid cell bounding box intersect
  if (intersects(triangle_bbox(coords), bbox)) {
    // if any of the triangle verts inside of bbox
    if (within_bbox(coords[0], bbox) ||
        within_bbox(coords[1], bbox) ||
        within_bbox(coords[2], bbox)) {
      return true;
    }
    // if any of the bbox verts are within the triangle
    if (bbox_verts_within_triangle(bbox, coords)) {
      return true;
    }
    // if any of the triangle's edges intersect with the bounding box
    if (line_intersects_bbox(coords[0], coords[1], bbox) ||
        line_intersects_bbox(coords[1], coords[2], bbox) ||
        line_intersects_bbox(coords[2], coords[0], bbox)) {
      return true;
    }
  }
  return false;
}

namespace detail
{
/**
 *  Functor for constructing the mapping from grid cells to intersecting
 * triangles. Each row of the resulting CSR structure represents a grid cell and
 * each row entry corresponds to an ID of an element that intersects that grid
 * cell.
 * \Warning since this works on Omega_h meshes, we currently assume each element
 * is a 2D simplex (triangle)
 * \Warning since this uses Omega_h data which is only available in the
 * "Default" Execution space, the should not be used in an alternative EXE space
 */
struct GridTriIntersectionFunctor
{
  GridTriIntersectionFunctor(Omega_h::Mesh& mesh, UniformGrid grid)
    : mesh_(mesh),
      tris2verts_(mesh_.ask_elem_verts()),
      coords_(mesh_.coords()),
      grid_(grid),
      nelems_(mesh_.nelems())
  {
    if (mesh_.dim() != 2) {
      std::cerr << "GridTriIntersection currently only developed for 2D "
                   "triangular meshes\n";
      std::terminate();
    }
  }
  /// Two-pass functor. On the first pass we set the number of grid/triangle
  /// intersections. On the second pass we fill the CSR array with the indexes
  /// to the triangles intersecting with the current grid cell (row)
  KOKKOS_INLINE_FUNCTION
  LO operator()(LO row, LO* fill) const
  {
    const auto grid_cell_bbox = grid_.GetCellBBOX(row);

    LO num_intersections = 0;

    // hierarchical parallel may make be very beneficial here...
    for (LO elem_idx = 0; elem_idx < nelems_; ++elem_idx) {
      const auto elem_tri2verts = Omega_h::gather_verts<3>(tris2verts_, elem_idx);
      // 2d mesh with 2d coords, but 3 triangles
      const auto vertex_coords = Omega_h::gather_vectors<3, 2>(coords_, elem_tri2verts);
      if (triangle_intersects_bbox(vertex_coords, grid_cell_bbox)) {
        if (fill) {
          fill[num_intersections] = elem_idx;
        }
        ++num_intersections;
      }
    }
    return num_intersections;
  }

private:
  Omega_h::Mesh& mesh_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
  UniformGrid grid_;
public:
  LO nelems_;
};

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh, const UniformGrid& grid)
{
  Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO> intersection_map{};
  auto f = detail::GridTriIntersectionFunctor{mesh, grid};
  Kokkos::count_and_fill_crs(intersection_map, grid.GetNumCells(), f);
  return intersection_map;
}
} // namespace detail

Omega_h::Vector<3> barycentric_from_global(
  const Omega_h::Vector<2>& point, const Omega_h::Matrix<2, 3>& vertex_coords)
{
  const auto inverse_basis =
    Omega_h::invert(Omega_h::simplex_basis<2, 2>(vertex_coords));
  auto xi = inverse_basis * (point - vertex_coords[0]);
  // note omega_h form_barycentric is currently broken.
  // see https://github.com/sandialabs/omega_h/issues/389
  return {1 - xi[0] - xi[1], xi[0], xi[1]};
}

std::pair<LO, Omega_h::Vector<GridPointSearch::dim + 1>>
GridPointSearch::operator()(Omega_h::Vector<dim> point) const
{
  // TODO use a functor for parallel for over a list of points
  // TODO return result struct rather than pair
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
    const auto vertex_coords =
      Omega_h::gather_vectors<3, 2>(coords_, elem_tri2verts);
    auto parametric_coords = barycentric_from_global(point, vertex_coords);
    constexpr auto fuzz = 1E-8;
    if (Omega_h::is_barycentric_inside(parametric_coords, fuzz)) {
      return {candidate_map_.entries(i), parametric_coords};
    }
  }
  return {-1, {}};
}
GridPointSearch::GridPointSearch(Omega_h::Mesh& mesh, LO Nx, LO Ny)
{
  auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
  // get mesh bounding box
  grid_ = {.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                           mesh_bbox.max[1] - mesh_bbox.min[1]},
    .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
    .divisions = {Nx, Ny}};
  candidate_map_ = detail::construct_intersection_map(mesh, grid_);
  coords_ = mesh.coords();
  tris2verts_ = mesh.ask_elem_verts();
  // TODO intiaialize coords and tri2vert and mesh!
  // TODO pass coords& tri2vert to construct map by references!
}
} // namespace wdmcpl
