#ifndef PCMS_LOCALIZATION_OVERLAP_H
#define PCMS_LOCALIZATION_OVERLAP_H
#include <Kokkos_Core.hpp>
#include <Omega_h_fail.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_shape.hpp>
#include <Omega_h_vector.hpp>
#include <r3d.hpp>

#include <cmath>

#include "pcms/utility/bounding_box.h"
#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"

namespace pcms
{

/// Predicate tags. Each tag names the operand pair with the exact type words
/// from the original function name (Triangle/Bbox, Simplex/Bbox,
/// Simplex/Simplex). A `_Method` suffix is added only when two tags would
/// otherwise collide.
///
/// Adding a new predicate is a three-step, local change: define a tag here,
/// declare its padding if it bins by bounding box (see below), and add an
/// `Overlap<Dim, NewTag>` specialization.
struct TriangleVsBbox
{
  /// Slack in grid cells that `candidate_cells` must add around the box.
  ///
  /// The predicate accepts a shared boundary plane as an intersection, so a
  /// cell that merely *touches* the box is a candidate. One cell is enough: a
  /// cell two or more away cannot touch the box at all.
  static constexpr LO boundary_pad = 1;
};
struct SimplexVsBbox
{
  /// \see TriangleVsBbox::boundary_pad
  static constexpr LO boundary_pad = 1;
};
struct SimplexVsSimplex
{
  // No `boundary_pad`
};

/// Inclusive per-axis range of grid cells overlapped by a bounding box.
///
/// The rules that define what a range *means* live here with the type rather
/// than in a consumer. Note that `CellRange{}` is not the empty range: its
/// zero initialization denotes the single cell at axis index 0, so use
/// `Empty()` when no cell should belong to it.
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

KOKKOS_INLINE_FUNCTION
AABBox<2> triangle_bbox(const Omega_h::Matrix<2, 3>& coords)
{
  Kokkos::Array<Real, 2> max{coords(0, 0), coords(1, 0)};
  Kokkos::Array<Real, 2> min{coords(0, 0), coords(1, 0)};
  for (int i = 1; i < 3; ++i) {
    max[0] = std::fmax(max[0], coords(0, i));
    max[1] = std::fmax(max[1], coords(1, i));
    min[0] = std::fmin(min[0], coords(0, i));
    min[1] = std::fmin(min[1], coords(1, i));
  }
  return {.center = {(max[0] + min[0]) / 2.0, (max[1] + min[1]) / 2.0},
          .half_width = {(max[0] - min[0]) / 2.0, (max[1] - min[1]) / 2.0}};
}

template <unsigned dim>
KOKKOS_INLINE_FUNCTION AABBox<dim> simplex_bbox(
  const Omega_h::Matrix<dim, dim + 1>& coords)
{
  Kokkos::Array<Real, dim> max;
  Kokkos::Array<Real, dim> min;
  for (int j = 0; j < dim; ++j) {
    max[j] = coords(j, 0);
    min[j] = coords(j, 0);
  }
  for (int i = 1; i < dim + 1; ++i) {
    for (int j = 0; j < dim; ++j) {
      max[j] = std::fmax(max[j], coords(j, i));
      min[j] = std::fmin(min[j], coords(j, i));
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

// Liang, You-Dong, and B. A. Barsky. "A New Concept and Method for Line
// Clipping." ACM Transactions on Graphics 3, no. 1 (January 1984): 1-22.
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
          // full liang-barksy algorithm computes new x-y coordinates here,
          // but we just need to check intersections
          return true;
        }
      }
    }
  }
  return false;
}

template <unsigned dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool within_bbox(
  const Omega_h::Vector<dim> coord, const AABBox<dim>& bbox) noexcept
{
  for (int i = 0; i < dim; ++i) {
    if (coord[i] < bbox.center[i] - bbox.half_width[i])
      return false;
    if (coord[i] > bbox.center[i] + bbox.half_width[i])
      return false;
  }
  return true;
}

template <int dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool bbox_verts_within_simplex(
  const AABBox<dim>& bbox, const Omega_h::Matrix<dim, dim + 1>& coords)
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
    Omega_h::Vector<dim> vert;
    for (unsigned j = 0; j < dim; ++j) {
      // eg 110 = 6 = [left, top, back]
      // eg 01 = 1 = [right, top]
      vert[j] = (i >> j) & 1 ? bbox_walls[j * 2] : bbox_walls[j * 2 + 1];
    }
    auto xi = Omega_h::barycentric_from_global<dim, dim>(vert, coords);
    if (Omega_h::is_barycentric_inside(xi)) {
      return true;
    }
  }
  return false;
}

} // namespace detail

/// An overlap predicate, selected as `Overlap<Dim, Tag>` and invoked through
/// its static `apply()`.
///
/// This is the single entry point for overlap tests: the primary template is
/// intentionally undefined, so each `Dim x Tag` combination that exists has an
/// explicit specialization whose `apply()` holds the original function body
/// verbatim, and asking for a combination that does not exist is a compile
/// error rather than a silent wrong-dimension call.
///
/// `Tag` names the operand pair (see the tags above); `Dim` fixes the
/// dimension. Because the tag is a template argument rather than a value, a tag
/// that carries no data is never passed as a dummy parameter.
template <int Dim, typename Tag>
struct Overlap;

template <>
struct Overlap<2, TriangleVsBbox>
{
  /**
   * Check if a triangle element represented by 3 coordinates in two dimensions
   * intersects with a bounding box
   */
  [[nodiscard]] KOKKOS_FUNCTION static bool apply(
    const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox)
  {
    // triangle and grid cell bounding box intersect
    if (intersects(detail::triangle_bbox(coords), bbox)) {
      // if any of the triangle verts inside of bbox
      if (detail::within_bbox<2>(coords[0], bbox) ||
          detail::within_bbox<2>(coords[1], bbox) ||
          detail::within_bbox<2>(coords[2], bbox)) {
        return true;
      }
      // if any of the bbox verts are within the triangle
      if (detail::bbox_verts_within_simplex(bbox, coords)) {
        return true;
      }
      // if any of the triangle's edges intersect with the bounding box
      if (detail::line_intersects_bbox(coords[0], coords[1], bbox) ||
          detail::line_intersects_bbox(coords[1], coords[2], bbox) ||
          detail::line_intersects_bbox(coords[2], coords[0], bbox)) {
        return true;
      }
    }
    return false;
  }
};

template <int Dim>
struct Overlap<Dim, SimplexVsBbox>
{
  [[nodiscard]] KOKKOS_FUNCTION static bool apply(
    const Omega_h::Matrix<Dim, Dim + 1>& coords, const AABBox<Dim>& bbox)
  {
    return intersects(detail::simplex_bbox<Dim>(coords), bbox);
    // TODO: Add refined cases from triangle_intersects_bbox
  }
};

template <int Dim>
struct Overlap<Dim, SimplexVsSimplex>
{
  /**
   * @brief Conservative test for a zero-volume overlap between two simplices.
   *
   * If verticies of one simplex lies in the closed outer half-space of some
   * face plane of the other, their intersection is contained in that plane and
   * so has no volume.
   *
   * @return true when the overlap is provably degenerate and the clip can be
   * skipped
   */
  [[nodiscard]] OMEGA_H_INLINE static bool apply(
    const r3d::Few<r3d::Vector<Dim>, Dim + 1>& a,
    const r3d::Few<r3d::Vector<Dim>, Dim + 1>& b)
  {
    // Each iteration treats one simplex as the half-space owner. Both are
    // needed: a face of `a` may separate the pair while no face of `b` does.
    for (int side = 0; side < 2; ++side) {
      const auto& owner = (side == 0) ? a : b;
      const auto& other = (side == 0) ? b : a;

      // Face `f` is the face opposite vertex `f`.
      for (int f = 0; f < Dim + 1; ++f) {
        // Vertices spanning the face, in index order with `f` skipped.
        int face_verts[Dim];
        int n_face_verts = 0;
        for (int v = 0; v < Dim + 1; ++v) {
          if (v != f) {
            face_verts[n_face_verts++] = v;
          }
        }
        const auto& origin = owner[face_verts[0]];

        // Face normal, oriented to point away from the opposite vertex so the
        // owner's interior lies on the negative side.
        r3d::Vector<Dim> normal;
        if constexpr (Dim == 2) {
          const auto edge = owner[face_verts[1]] - origin;
          normal[0] = edge[1];
          normal[1] = -edge[0];
        } else {
          const auto e1 = owner[face_verts[1]] - origin;
          const auto e2 = owner[face_verts[2]] - origin;
          normal[0] = e1[1] * e2[2] - e1[2] * e2[1];
          normal[1] = e1[2] * e2[0] - e1[0] * e2[2];
          normal[2] = e1[0] * e2[1] - e1[1] * e2[0];
        }

        Omega_h::Real opposite_side = 0.0;
        Omega_h::Real normal_sq = 0.0;
        for (int d = 0; d < Dim; ++d) {
          opposite_side += normal[d] * (owner[f][d] - origin[d]);
          normal_sq += normal[d] * normal[d];
        }
        // A degenerate face gives no usable axis.
        if (normal_sq == 0.0) {
          continue;
        }
        const Omega_h::Real orientation = (opposite_side > 0.0) ? -1.0 : 1.0;

        bool all_outside = true;
        for (int v = 0; v < Dim + 1 && all_outside; ++v) {
          Omega_h::Real signed_distance = 0.0;
          for (int d = 0; d < Dim; ++d) {
            signed_distance += normal[d] * (other[v][d] - origin[d]);
          }
          all_outside = (orientation * signed_distance) >= 0.0;
        }
        if (all_outside) {
          return true;
        }
      }
    }
    return false;
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

/// Cells on `grid` that could satisfy a `Pred` test against `simplex`: the
/// cells its bounding box spans, widened by that predicate's own `boundary_pad`
/// so that a cell which merely touches the box is included.
template <int Dim, typename Pred>
[[nodiscard]] KOKKOS_INLINE_FUNCTION CellRange<Dim> candidate_cells(
  const Omega_h::Matrix<Dim, Dim + 1>& simplex, const UniformGrid<Dim>& grid)
{
  return overlapping_cells<Dim>(detail::simplex_bbox<Dim>(simplex), grid,
                                Pred::boundary_pad);
}

/// As above, for a bare `box`: the range a closed-boundary predicate `Pred`
/// could accept within, for callers such as routing a query point to the ranks
/// whose published ranges it could fall in.
template <int Dim, typename Pred>
[[nodiscard]] KOKKOS_INLINE_FUNCTION CellRange<Dim> candidate_cells(
  const AABBox<Dim>& box, const UniformGrid<Dim>& grid)
{
  return overlapping_cells<Dim>(box, grid, Pred::boundary_pad);
}

} // namespace pcms

#endif // PCMS_LOCALIZATION_OVERLAP_H
