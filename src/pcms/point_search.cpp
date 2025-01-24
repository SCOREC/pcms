#include "point_search.h"
#include <Omega_h_mesh.hpp>
#include <bitset>

// From https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
KOKKOS_INLINE_FUNCTION
double distance_from_line(const double x0, const double y0, const double x1, const double y1, const double x2, const double y2)
{
  const Omega_h::Vector<2> p1 = { x1, y1 };
  const Omega_h::Vector<2> p2 = { x2, y2 };
  auto disp = p2 - p1;

  return std::abs(disp[1]*x0 - disp[0]*y0 + x2*y1 - y2*x1) / Omega_h::norm(disp);
}

// Law of Cosines, where a, b, c and gamma are defined here: https://en.wikipedia.org/wiki/Law_of_cosines#Use_in_solving_triangles
KOKKOS_INLINE_FUNCTION
double angle_from_side_lengths(const double a, const double b, const double c)
{
  return std::acos((a*a + b*b - c*c) / 2*a*b);
}

KOKKOS_INLINE_FUNCTION
bool normal_intersects_segment(const Omega_h::Few<double, 2> a, const Omega_h::Few<double, 2> b, const Omega_h::Few<double, 2> c)
{
  const auto ab_len = Omega_h::norm(a - b);
  const auto bc_len = Omega_h::norm(b -c);
  const auto ac_len = Omega_h::norm(a - c);

  const double angle1 = angle_from_side_lengths(bc_len, ac_len, ab_len);
  const double angle2 = angle_from_side_lengths(bc_len, ab_len, ac_len);

  return angle1 <= (Omega_h::PI / 2) && angle2 <= (Omega_h::PI / 2);
}

namespace pcms
{
constexpr Real fuzz = 1E-6;

KOKKOS_INLINE_FUNCTION
AABBox<2> triangle_bbox(const Omega_h::Matrix<2, 3>& coords)
{
  std::array<Real, 2> max{coords(0, 0), coords(1, 0)};
  std::array<Real, 2> min{coords(0, 0), coords(1, 0)};
  for (int i = 1; i < 3; ++i) {
    max[0] = std::fmax(max[0], coords(0, i));
    max[1] = std::fmax(max[1], coords(1, i));
    min[0] = std::fmin(min[0], coords(0, i));
    min[1] = std::fmin(min[1], coords(1, i));
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
  if(Omega_h::is_barycentric_inside(xi, fuzz)){ return true;}
  xi = barycentric_from_global({left,top}, coords);
  if(Omega_h::is_barycentric_inside(xi, fuzz)){ return true;}
  xi = barycentric_from_global({right,top}, coords);
  if(Omega_h::is_barycentric_inside(xi, fuzz)){ return true;}
  xi = barycentric_from_global({right,bot}, coords);
  if(Omega_h::is_barycentric_inside(xi, fuzz)){ return true;}
  return false;
}

/**
 * Check if a triangle element represented by 3 coordinates in two dimensions
 * intersects with a bounding box
 */
[[nodiscard]]
KOKKOS_FUNCTION
bool triangle_intersects_bbox(const Omega_h::Matrix<2, 3>& coords,
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
  GridTriIntersectionFunctor(Omega_h::Mesh& mesh, Kokkos::View<UniformGrid[1]> grid)
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
    const auto grid_cell_bbox = grid_(0).GetCellBBOX(row);
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
  Kokkos::View<UniformGrid[1]> grid_;
public:
  LO nelems_;
};

// num_grid_cells should be result of grid.GetNumCells(), take as argument to avoid extra copy
// of grid from gpu to cpu
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh, Kokkos::View<UniformGrid[1]> grid, int num_grid_cells)
{
  Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO> intersection_map{};
  auto f = detail::GridTriIntersectionFunctor{mesh, grid};
  Kokkos::count_and_fill_crs(intersection_map, num_grid_cells, f);
  return intersection_map;
}
} // namespace detail

KOKKOS_FUNCTION
Omega_h::Vector<3> barycentric_from_global(
  const Omega_h::Vector<2>& point, const Omega_h::Matrix<2, 3>& vertex_coords)
{
  const auto inverse_basis =
    Omega_h::pseudo_invert(Omega_h::simplex_basis<2, 2>(vertex_coords));
  auto xi = inverse_basis * (point - vertex_coords[0]);
  // note omega_h form_barycentric is currently broken.
  // see https://github.com/sandialabs/omega_h/issues/389
  return {1 - xi[0] - xi[1], xi[0], xi[1]};
}

template <int n,  typename Op>
OMEGA_H_INLINE double myreduce(const Omega_h::Vector<n> & x, Op op) OMEGA_H_NOEXCEPT {
  auto out = x[0];
  for (int i = 1; i < n; ++i) out = op(out, x[i]);
  return out;
}         

Kokkos::View<GridPointSearch::Result*> GridPointSearch::operator()(Kokkos::View<Real*[dim] > points) const
{
  static_assert(dim == 2, "point search assumes dim==2");
  Kokkos::View<GridPointSearch::Result*> results("point search result", points.extent(0));
  auto num_rows = candidate_map_.numRows();
  // needed so that we don't capture this ptr which will be memory error on cuda
  auto grid = grid_; 
  auto candidate_map = candidate_map_;
  auto tris2verts = tris2verts_;
  auto tris2verts_adj = tris2verts_adj_;
  auto tris2edges_adj = tris2edges_adj_;
  auto edges2verts_adj = edges2verts_adj_;
  auto coords = coords_;
  Kokkos::parallel_for(points.extent(0), KOKKOS_LAMBDA(int p) {
    Omega_h::Vector<2> point(std::initializer_list<double>{points(p,0), points(p,1)});
    auto cell_id = grid(0).ClosestCellID(point);
    assert(cell_id < num_rows && cell_id >= 0);
    auto candidates_begin = candidate_map.row_map(cell_id);
    auto candidates_end = candidate_map.row_map(cell_id + 1);
    bool found = false;

    auto nearest_triangle = candidates_begin;
    auto dimensionality = GridPointSearch::Result::Dimensionality::EDGE;
    Omega_h::Real distance_to_nearest { INFINITY };
    Omega_h::Vector<3> parametric_coords_to_nearest;
    // create array that's size of number of candidates x num coords to store
    // parametric inversion
    for (auto i = candidates_begin; i < candidates_end; ++i) {
      int triangleID = candidate_map.entries(i);
      auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts, triangleID);
      // 2d mesh with 2d coords, but 3 triangles
     auto vertex_coords = Omega_h::gather_vectors<3, 2>(coords, elem_tri2verts);
     auto parametric_coords = barycentric_from_global(point, vertex_coords);

     for (int j = 0; j < 3; ++j) {
       // Every triangle (face) is connected to 3 edges
       const int edgeID = tris2edges_adj.ab2b[triangleID * 3 + j];

       auto vertex_a_id = edges2verts_adj.ab2b[edgeID * 2];
       auto vertex_b_id = edges2verts_adj.ab2b[edgeID * 2 + 1];

       auto vertex_a = Omega_h::get_vector<2>(coords, vertex_a_id);
       auto vertex_b = Omega_h::get_vector<2>(coords, vertex_b_id);

       if (!normal_intersects_segment(point, vertex_a, vertex_b)) continue;

       const auto xa = vertex_a[0];
       const auto ya = vertex_a[1];
       const auto xb = vertex_b[0];
       const auto yb = vertex_b[1];

       const auto xp = point[0];
       const auto yp = point[1];

       const auto distance_to_ab = distance_from_line(xp, yp, xa, ya, xb, yb);

       if (distance_to_ab >= distance_to_nearest) { continue; }

       dimensionality = GridPointSearch::Result::Dimensionality::EDGE;
       nearest_triangle = i;
       distance_to_nearest = distance_to_ab;
       parametric_coords_to_nearest = parametric_coords;
      }

      // Every triangle (face) is connected to 3 vertices
      for (int j = 0; j < 3; ++j) {
        // Get the vertex ID from the connectivity array
        const int vertexID = tris2verts_adj_.ab2b[triangleID * 3 + j];
        // Get the vertex coordinates from the mesh using vertexID
        const Omega_h::Few<double, 2> vertex =
          Omega_h::get_vector<2>(coords, vertexID);

        if (const auto distance = Omega_h::norm(point - vertex);distance < distance_to_nearest) {
         dimensionality = GridPointSearch::Result::Dimensionality::VERTEX;
         nearest_triangle = i;
         distance_to_nearest = distance;
         parametric_coords_to_nearest = parametric_coords;
       }
      }

     if (Omega_h::is_barycentric_inside(parametric_coords, fuzz)) {
       results(p) = GridPointSearch::Result{GridPointSearch::Result::Dimensionality::FACE, triangleID, parametric_coords};
       found = true;
       break;
     }
    }
    if(!found)
    {
      results(p) = GridPointSearch::Result{dimensionality, -1 * candidate_map.entries(nearest_triangle), parametric_coords_to_nearest};
    }
  });

  return results;
}

GridPointSearch::GridPointSearch(Omega_h::Mesh& mesh, LO Nx, LO Ny)
{
  auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
  auto grid_h = Kokkos::create_mirror_view(grid_);
  grid_h(0) = UniformGrid{.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                           mesh_bbox.max[1] - mesh_bbox.min[1]},
    .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
    .divisions = {Nx, Ny}};
  Kokkos::deep_copy(grid_, grid_h);
  candidate_map_ = detail::construct_intersection_map(mesh, grid_, grid_h(0).GetNumCells());
  coords_ = mesh.coords();
  tris2verts_ = mesh.ask_elem_verts();
  tris2edges_adj_ = mesh.ask_down(Omega_h::FACE, Omega_h::EDGE);
  tris2verts_adj_ = mesh.ask_down(Omega_h::FACE, Omega_h::VERT);
  edges2verts_adj_ = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT);
}
} // namespace pcms
