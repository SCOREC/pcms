#ifndef PCMS_UTILITY_OVERLAP_HPP
#define PCMS_UTILITY_OVERLAP_HPP
#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstddef>

#include "pcms/utility/shapes.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

namespace pcms
{

/// Inclusive per-axis range of grid cells overlapped by a bounding box.
template <int Dim>
struct CellRange
{
  Kokkos::Array<LO, Dim> lo;
  Kokkos::Array<LO, Dim> hi;

  /// The range no cell belongs to, as published by a rank with no geometry.
  [[nodiscard]] static CellRange Empty() noexcept
  {
    CellRange range{};
    for (int d = 0; d < Dim; ++d) {
      range.lo[d] = 0;
      range.hi[d] = -1;
    }
    return range;
  }

  /// True when no cell is in the range, so no query can match it.
  [[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsEmpty() const noexcept
  {
    for (int d = 0; d < Dim; ++d) {
      if (hi[d] < lo[d]) {
        return true;
      }
    }
    return false;
  }

  /// True when the cell at per-axis indices `axis_indices` falls in the range.
  [[nodiscard]] KOKKOS_INLINE_FUNCTION bool Contains(
    const Kokkos::Array<LO, Dim>& axis_indices) const noexcept
  {
    // An empty range has hi < lo on every axis, so this is false for it too.
    for (int d = 0; d < Dim; ++d) {
      if (axis_indices[d] < lo[d] || axis_indices[d] > hi[d]) {
        return false;
      }
    }
    return true;
  }
};

namespace detail
{

// The vertices of a simplex, vertices[k][d] the coordinate d of vertex k.
template <int Dim>
using SimplexVertices = Kokkos::Array<Kokkos::Array<Real, Dim>, Dim + 1>;

// A Dim x Dim matrix stored as its columns, m[j][i] the entry in row i of
// column j, as Omega_h::Matrix.
template <int Dim>
using SquareMatrix = Kokkos::Array<Kokkos::Array<Real, Dim>, Dim>;

[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, 3> cross(
  const Kokkos::Array<Real, 3>& u, const Kokkos::Array<Real, 3>& v)
{
  return {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
          u[0] * v[1] - u[1] * v[0]};
}

// The edge vectors p[k + 1] - p[0] of a simplex, as the columns of a matrix.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION SquareMatrix<Dim> simplex_basis(
  const SimplexVertices<Dim>& p)
{
  SquareMatrix<Dim> b;
  for (int k = 0; k < Dim; ++k) {
    for (int d = 0; d < Dim; ++d) {
      b[k][d] = p[k + 1][d] - p[0][d];
    }
  }
  return b;
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION Real determinant(const SquareMatrix<2>& m)
{
  const Real a = m[0][0];
  const Real b = m[1][0];
  const Real c = m[0][1];
  const Real d = m[1][1];
  return a * d - b * c;
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION Real determinant(const SquareMatrix<3>& m)
{
  const Real a = m[0][0];
  const Real b = m[1][0];
  const Real c = m[2][0];
  const Real d = m[0][1];
  const Real e = m[1][1];
  const Real f = m[2][1];
  const Real g = m[0][2];
  const Real h = m[1][2];
  const Real i = m[2][2];
  return (a * e * i) + (b * f * g) + (c * d * h) - (c * e * g) - (b * d * i) -
         (a * f * h);
}

// The inverses of 1x1 to 3x3 matrices, as Omega_h::invert.
[[nodiscard]] KOKKOS_INLINE_FUNCTION SquareMatrix<1> invert(
  const SquareMatrix<1>& m)
{
  return {{{1.0 / m[0][0]}}};
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION SquareMatrix<2> invert(
  const SquareMatrix<2>& m)
{
  const Real a = m[0][0];
  const Real b = m[1][0];
  const Real c = m[0][1];
  const Real d = m[1][1];
  const Real det = determinant(m);
  // columns of [[d, -b], [-c, a]] / det
  return {{{d / det, -c / det}, {-b / det, a / det}}};
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION SquareMatrix<3> invert(
  const SquareMatrix<3>& m)
{
  // row k of the inverse is the cross product of the other two columns, over
  // the determinant
  const SquareMatrix<3> rows{cross(m[1], m[2]), cross(m[2], m[0]),
                             cross(m[0], m[1])};
  const Real det = determinant(m);
  SquareMatrix<3> inverse;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      inverse[j][i] = rows[i][j] / det;
    }
  }
  return inverse;
}

// m * v, accumulated column by column.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim> multiply(
  const SquareMatrix<Dim>& m, const Kokkos::Array<Real, Dim>& v)
{
  Kokkos::Array<Real, Dim> out;
  for (int i = 0; i < Dim; ++i) {
    out[i] = m[0][i] * v[0];
  }
  for (int j = 1; j < Dim; ++j) {
    for (int i = 0; i < Dim; ++i) {
      out[i] = out[i] + m[j][i] * v[j];
    }
  }
  return out;
}

// The barycentric coordinates from the last Dim of them: the first is 1 minus
// their sum.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim + 1>
form_barycentric(const Kokkos::Array<Real, Dim>& c)
{
  Kokkos::Array<Real, Dim + 1> xi;
  xi[0] = 1.0;
  for (int k = 1; k < Dim + 1; ++k) {
    xi[k] = c[k - 1];
    xi[0] -= c[k - 1];
  }
  return xi;
}

// The barycentric coordinates of `x` in the simplex `p`.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim + 1>
barycentric_from_global(const Kokkos::Array<Real, Dim>& x,
                        const SimplexVertices<Dim>& p)
{
  static_assert(Dim >= 1 && Dim <= 3,
                "barycentric_from_global: only 1 to 3 dimensions");
  Kokkos::Array<Real, Dim> r;
  for (int d = 0; d < Dim; ++d) {
    r[d] = x[d] - p[0][d];
  }
  const auto inverse_basis = invert(simplex_basis<Dim>(p));
  return form_barycentric<Dim>(multiply<Dim>(inverse_basis, r));
}

// Whether barycentric coordinates are inside the simplex, within `fuzz`, as
// Omega_h::is_barycentric_inside.
template <std::size_t N>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool is_barycentric_inside(
  const Kokkos::Array<Real, N>& xi, Real fuzz = 0)
{
  Real lo = xi[0];
  Real hi = xi[0];
  for (std::size_t k = 1; k < N; ++k) {
    lo = Kokkos::min(lo, xi[k]);
    hi = Kokkos::max(hi, xi[k]);
  }
  return (0.0 - fuzz) <= lo && hi <= (1.0 + fuzz);
}

KOKKOS_INLINE_FUNCTION
AABBox<2> triangle_bbox(const SimplexVertices<2>& coords)
{
  Kokkos::Array<Real, 2> max{coords[0][0], coords[0][1]};
  Kokkos::Array<Real, 2> min{coords[0][0], coords[0][1]};
  for (int i = 1; i < 3; ++i) {
    max[0] = std::fmax(max[0], coords[i][0]);
    max[1] = std::fmax(max[1], coords[i][1]);
    min[0] = std::fmin(min[0], coords[i][0]);
    min[1] = std::fmin(min[1], coords[i][1]);
  }
  return {.center = {(max[0] + min[0]) / 2.0, (max[1] + min[1]) / 2.0},
          .half_width = {(max[0] - min[0]) / 2.0, (max[1] - min[1]) / 2.0}};
}

template <unsigned dim>
KOKKOS_INLINE_FUNCTION AABBox<dim> simplex_bbox(
  const SimplexVertices<dim>& coords)
{
  Kokkos::Array<Real, dim> max;
  Kokkos::Array<Real, dim> min;
  for (int j = 0; j < dim; ++j) {
    max[j] = coords[0][j];
    min[j] = coords[0][j];
  }
  for (int i = 1; i < dim + 1; ++i) {
    for (int j = 0; j < dim; ++j) {
      max[j] = std::fmax(max[j], coords[i][j]);
      min[j] = std::fmin(min[j], coords[i][j]);
    }
  }

  Kokkos::Array<Real, dim> center;
  Kokkos::Array<Real, dim> half_width;

  for (int j = 0; j < dim; ++j) {
    center[j] = (max[j] + min[j]) / 2.0;
    half_width[j] = (max[j] - min[j]) / 2.0;
  }

  return {.center = center, .half_width = half_width};
}

// Liang, You-Dong, and B. A. Barsky. "A New Concept and Method for Line
// Clipping." ACM Transactions on Graphics 3, no. 1 (January 1984): 1-22.
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

template <int dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool bbox_verts_within_simplex(
  const AABBox<dim>& bbox, const SimplexVertices<dim>& coords)
{
  // each dimension has a pair of opposing "walls"
  // 2D: { [left, right], [top, bottom] } -> { left, right, top, bottom }
  // 3D: { [left, right], [top, bottom], [front, back] } -> { left, ..., back }
  Kokkos::Array<Real, dim * 2ul> bbox_walls{};
  for (int i = 0; i < dim; i++) {
    bbox_walls[i * 2] = bbox.center[i] - bbox.half_width[i];
    bbox_walls[i * 2 + 1] = bbox.center[i] + bbox.half_width[i];
  }

  // 1 << dim == 2 ** dim == num vertices in ndim bounding box / hypercube
  constexpr unsigned num_verts = 1 << dim;

  // Each vertex is a just a unique combination of walls
  // eg [left, bottom] (2D) or [right, top, front] (3)
  for (unsigned i = 0; i < num_verts; ++i) {
    // conveniently, i acts a bit field representing the current combination
    Kokkos::Array<Real, dim> vert;
    for (unsigned j = 0; j < dim; ++j) {
      // eg 110 = 6 = [left, top, back]
      // eg 01 = 1 = [right, top]
      vert[j] = (i >> j) & 1 ? bbox_walls[j * 2] : bbox_walls[j * 2 + 1];
    }
    auto xi = barycentric_from_global<dim>(vert, coords);
    if (is_barycentric_inside(xi)) {
      return true;
    }
  }
  return false;
}

} // namespace detail

/**
 * The smallest axis-aligned box containing a shape.
 */
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<Dim> bounding_box(
  const Simplex<Dim>& simplex)
{
  return detail::simplex_bbox<Dim>(simplex.vertices);
}

template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<Dim> bounding_box(
  const AABBox<Dim>& box)
{
  return box;
}

template <typename TagA, typename TagB, typename A, typename B>
struct Overlap;

/**
 * True when `a` and `b` overlaps.
 */
template <typename A, typename B>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool overlaps(const A& a, const B& b)
{
  static_assert(dimension_v<A> == dimension_v<B>,
                "overlaps: operands have different dimensions");
  using Forward = Overlap<tag_t<A>, tag_t<B>, A, B>;
  using Reverse = Overlap<tag_t<B>, tag_t<A>, B, A>;
  if constexpr (requires { Forward::apply(a, b); }) {
    return Forward::apply(a, b);
  } else {
    static_assert(
      requires { Reverse::apply(b, a); },
      "overlaps: no overlap test for this pair of shapes");
    return Reverse::apply(b, a);
  }
}

template <typename A, typename B>
struct Overlap<BoxTag, BoxTag, A, B>
{
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool apply(const A& a, const B& b)
  {
    for (int i = 0; i < dimension_v<A>; ++i) {
      if (std::abs(a.center[i] - b.center[i]) >
          (a.half_width[i] + b.half_width[i]))
        return false;
    }
    return true;
  }
};

template <typename A, typename B>
struct Overlap<PointTag, BoxTag, A, B>
{
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool apply(const A& point,
                                                         const B& bbox)
  {
    for (int i = 0; i < dimension_v<A>; ++i) {
      if (point.coords[i] < bbox.center[i] - bbox.half_width[i])
        return false;
      if (point.coords[i] > bbox.center[i] + bbox.half_width[i])
        return false;
    }
    return true;
  }
};

template <typename A, typename B>
struct Overlap<PointTag, SimplexTag, A, B>
{
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool apply(const A& point,
                                                         const B& simplex)
  {
    constexpr int dim = dimension_v<A>;
    const auto xi =
      detail::barycentric_from_global<dim>(point.coords, simplex.vertices);
    return detail::is_barycentric_inside(xi);
  }
};

// Liang, You-Dong, and B. A. Barsky. "A New Concept and Method for Line
// Clipping." ACM Transactions on Graphics 3, no. 1 (January 1984): 1-22.
// https://doi.org/10.1145/357332.357333.
template <typename A, typename B>
struct Overlap<SegmentTag, BoxTag, A, B>
{
  static_assert(dimension_v<A> == 2,
                "segment-box overlap is implemented for 2D only");

  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool apply(const A& seg,
                                                         const B& bbox)
  {
    Real t0 = 0;
    Real t1 = 1;
    auto xleft = bbox.center[0] - bbox.half_width[0];
    auto deltax = seg.end[0] - seg.start[0];
    if (detail::clipt(-deltax, seg.start[0] - xleft, t0, t1)) {
      auto xright = bbox.center[0] + bbox.half_width[0];
      if (detail::clipt(deltax, xright - seg.start[0], t0, t1)) {
        auto deltay = seg.end[1] - seg.start[1];
        auto ybottom = bbox.center[1] - bbox.half_width[1];
        if (detail::clipt(-deltay, seg.start[1] - ybottom, t0, t1)) {
          auto ytop = bbox.center[1] + bbox.half_width[1];
          if (detail::clipt(deltay, ytop - seg.start[1], t0, t1)) {
            // full liang-barksy algorithm computes new x-y coordinates here,
            // but we just need to check intersections
            return true;
          }
        }
      }
    }
    return false;
  }
};

template <typename A, typename B>
struct Overlap<SimplexTag, BoxTag, A, B>
{
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool apply(const A& simplex,
                                                         const B& bbox)
  {
    if constexpr (dimension_v<A> == 2) {
      const auto& coords = simplex.vertices;
      // triangle and grid cell bounding box intersect
      if (overlaps(detail::triangle_bbox(coords), bbox)) {
        // if any of the triangle verts inside of bbox
        if (overlaps(simplex.Vertex(0), bbox) ||
            overlaps(simplex.Vertex(1), bbox) ||
            overlaps(simplex.Vertex(2), bbox)) {
          return true;
        }
        // if any of the bbox verts are within the triangle
        if (detail::bbox_verts_within_simplex<2>(bbox, coords)) {
          return true;
        }
        // if any of the triangle's edges intersect with the bounding box
        if (overlaps(simplex.Edge(0, 1), bbox) ||
            overlaps(simplex.Edge(1, 2), bbox) ||
            overlaps(simplex.Edge(2, 0), bbox)) {
          return true;
        }
      }
      return false;
    } else {
      // TODO: Add refined cases from the 2D test
      return overlaps(detail::simplex_bbox<dimension_v<A>>(simplex.vertices),
                      bbox);
    }
  }
};

/// Low-level cell range: cells overlapped by `box` on `grid`, widened by `pad`
/// cells per side and clamped to the grid.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION CellRange<Dim> overlapping_cells(
  const AABBox<Dim>& box, const UniformGrid<Dim>& grid, LO pad = 1)
{
  CellRange<Dim> range;
  for (int d = 0; d < Dim; ++d) {
    range.lo[d] = Kokkos::max(
      grid.AxisCellIndex(d, box.center[d] - box.half_width[d]) - pad, 0);
    range.hi[d] = Kokkos::min(
      grid.AxisCellIndex(d, box.center[d] + box.half_width[d]) + pad,
      grid.divisions[d] - 1);
  }
  return range;
}

// Cells on `grid` that a shape inside `bound` could overlap.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION CellRange<Dim> candidate_cells(
  const AABBox<Dim>& bound, const UniformGrid<Dim>& grid)
{
  return overlapping_cells<Dim>(bound, grid, 1);
}

} // namespace pcms

#endif // PCMS_UTILITY_OVERLAP_HPP
