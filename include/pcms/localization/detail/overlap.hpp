#ifndef PCMS_LOCALIZATION_DETAIL_OVERLAP_HPP
#define PCMS_LOCALIZATION_DETAIL_OVERLAP_HPP
#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstddef>

#include "pcms/localization/shapes.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/detail/small_matrix.hpp"
#include "pcms/utility/types.hpp"

namespace pcms::detail
{

// The vertices of a simplex, vertices[k][d] the coordinate d of vertex k.
template <int Dim>
using SimplexVertices = Kokkos::Array<Kokkos::Array<Real, Dim>, Dim + 1>;

// The edge vectors p[k + 1] - p[0] of a simplex, as the columns of a matrix.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION SquareMatrix<Dim> SimplexBasis(
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

// The barycentric coordinates from the last Dim of them: the first is 1 minus
// their sum.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim + 1>
FormBarycentric(const Kokkos::Array<Real, Dim>& c)
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
BarycentricFromGlobal(const Kokkos::Array<Real, Dim>& x,
                      const SimplexVertices<Dim>& p)
{
  static_assert(Dim >= 1 && Dim <= 3,
                "BarycentricFromGlobal: only 1 to 3 dimensions");
  Kokkos::Array<Real, Dim> r;
  for (int d = 0; d < Dim; ++d) {
    r[d] = x[d] - p[0][d];
  }
  const auto inverse_basis = Invert<Dim>(SimplexBasis<Dim>(p));
  return FormBarycentric<Dim>(Multiply<Dim>(inverse_basis, r));
}

// Whether barycentric coordinates are inside the simplex, within `fuzz`, as
// Omega_h::is_barycentric_inside.
template <std::size_t N>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsBarycentricInside(
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
AABBox<2> TriangleBBox(const SimplexVertices<2>& coords)
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
KOKKOS_INLINE_FUNCTION AABBox<dim> SimplexBBox(
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
bool ClipT(Real p, Real q, Real& t0, Real& t1)
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
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool BBoxVertsWithinSimplex(
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
    auto xi = BarycentricFromGlobal<dim>(vert, coords);
    if (IsBarycentricInside(xi)) {
      return true;
    }
  }
  return false;
}

// The overlap test of a pair of shapes, by their tags.
template <typename TagA, typename TagB, typename A, typename B>
struct Overlap;

template <typename A, typename B>
struct Overlap<BoxTag, BoxTag, A, B>
{
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool Apply(const A& a, const B& b)
  {
    for (int i = 0; i < ShapeDimension<A>; ++i) {
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
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool Apply(const A& point,
                                                         const B& bbox)
  {
    for (int i = 0; i < ShapeDimension<A>; ++i) {
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
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool Apply(const A& point,
                                                         const B& simplex)
  {
    constexpr int dim = ShapeDimension<A>;
    const auto xi = BarycentricFromGlobal<dim>(point.coords, simplex.vertices);
    return IsBarycentricInside(xi);
  }
};

// Liang, You-Dong, and B. A. Barsky. "A New Concept and Method for Line
// Clipping." ACM Transactions on Graphics 3, no. 1 (January 1984): 1-22.
// https://doi.org/10.1145/357332.357333.
template <typename A, typename B>
struct Overlap<SegmentTag, BoxTag, A, B>
{
  static_assert(ShapeDimension<A> == 2,
                "segment-box overlap is implemented for 2D only");

  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool Apply(const A& seg,
                                                         const B& bbox)
  {
    Real t0 = 0;
    Real t1 = 1;
    auto xleft = bbox.center[0] - bbox.half_width[0];
    auto deltax = seg.end[0] - seg.start[0];
    if (ClipT(-deltax, seg.start[0] - xleft, t0, t1)) {
      auto xright = bbox.center[0] + bbox.half_width[0];
      if (ClipT(deltax, xright - seg.start[0], t0, t1)) {
        auto deltay = seg.end[1] - seg.start[1];
        auto ybottom = bbox.center[1] - bbox.half_width[1];
        if (ClipT(-deltay, seg.start[1] - ybottom, t0, t1)) {
          auto ytop = bbox.center[1] + bbox.half_width[1];
          if (ClipT(deltay, ytop - seg.start[1], t0, t1)) {
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
  [[nodiscard]] KOKKOS_INLINE_FUNCTION static bool Apply(const A& simplex,
                                                         const B& bbox)
  {
    if constexpr (ShapeDimension<A> == 2) {
      const auto& coords = simplex.vertices;
      // triangle and grid cell bounding box intersect
      if (Overlaps(TriangleBBox(coords), bbox)) {
        // if any of the triangle verts inside of bbox
        if (Overlaps(simplex.Vertex(0), bbox) ||
            Overlaps(simplex.Vertex(1), bbox) ||
            Overlaps(simplex.Vertex(2), bbox)) {
          return true;
        }
        // if any of the bbox verts are within the triangle
        if (BBoxVertsWithinSimplex<2>(bbox, coords)) {
          return true;
        }
        // if any of the triangle's edges intersect with the bounding box
        if (Overlaps(simplex.Edge(0, 1), bbox) ||
            Overlaps(simplex.Edge(1, 2), bbox) ||
            Overlaps(simplex.Edge(2, 0), bbox)) {
          return true;
        }
      }
      return false;
    } else {
      // TODO: Add refined cases from the 2D test
      return Overlaps(SimplexBBox<ShapeDimension<A>>(simplex.vertices), bbox);
    }
  }
};

} // namespace pcms::detail

#endif // PCMS_LOCALIZATION_DETAIL_OVERLAP_HPP
