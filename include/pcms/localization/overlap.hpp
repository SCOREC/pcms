#ifndef PCMS_LOCALIZATION_OVERLAP_HPP
#define PCMS_LOCALIZATION_OVERLAP_HPP
#include <Kokkos_Core.hpp>

#include "pcms/localization/detail/overlap.hpp"
#include "pcms/localization/shapes.hpp"
#include "pcms/utility/bounding_box.hpp"
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

/**
 * The smallest axis-aligned box containing a shape.
 */
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<Dim> BoundingBox(
  const Simplex<Dim>& simplex)
{
  return detail::SimplexBBox<Dim>(simplex.vertices);
}

template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<Dim> BoundingBox(
  const AABBox<Dim>& box)
{
  return box;
}

/**
 * True when `a` and `b` overlap.
 */
template <typename A, typename B>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool Overlaps(const A& a, const B& b)
{
  static_assert(ShapeDimension<A> == ShapeDimension<B>,
                "Overlaps: operands have different dimensions");
  using TagA = typename detail::ShapeTraits<A>::tag;
  using TagB = typename detail::ShapeTraits<B>::tag;
  using Forward = detail::Overlap<TagA, TagB, A, B>;
  using Reverse = detail::Overlap<TagB, TagA, B, A>;
  if constexpr (requires { Forward::Apply(a, b); }) {
    return Forward::Apply(a, b);
  } else {
    static_assert(
      requires { Reverse::Apply(b, a); },
      "Overlaps: no overlap test for this pair of shapes");
    return Reverse::Apply(b, a);
  }
}

/**
 * The cells of `grid` overlapped by `box`, widened by `pad` cells per side
 * and clamped to the grid. The default pad of one cell covers a shape inside
 * `box` that touches a cell boundary.
 */
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION CellRange<Dim> OverlappingCells(
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

} // namespace pcms

#endif // PCMS_LOCALIZATION_OVERLAP_HPP
