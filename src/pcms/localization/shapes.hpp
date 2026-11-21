#ifndef PCMS_LOCALIZATION_SHAPES_HPP
#define PCMS_LOCALIZATION_SHAPES_HPP

#include <Kokkos_Core.hpp>

#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/types.hpp"

namespace pcms
{

// Shape tags: what kind of shape a type is, for choosing an overlap test.
struct PointTag
{};
struct BoxTag
{};
struct SimplexTag
{};
struct SegmentTag
{};

/** A point in `Dim`-dimensional space. */
template <int Dim>
struct Point
{
  Kokkos::Array<Real, Dim> coords;
};

/** A line segment in `Dim`-dimensional space, from `start` to `end`. */
template <int Dim>
struct Segment
{
  Kokkos::Array<Real, Dim> start;
  Kokkos::Array<Real, Dim> end;
};

/**
 * A simplex in `Dim`-dimensional space: a triangle in 2D, a tetrahedron in 3D.
 */
template <int Dim>
struct Simplex
{
  // vertices[k][d]: coordinate d of vertex k
  Kokkos::Array<Kokkos::Array<Real, Dim>, Dim + 1> vertices;

  /** Vertex `k`, for `k` in `[0, Dim]`. */
  [[nodiscard]] KOKKOS_INLINE_FUNCTION Point<Dim> Vertex(int k) const
  {
    return {vertices[k]};
  }

  /** The edge from vertex `i` to vertex `j`. */
  [[nodiscard]] KOKKOS_INLINE_FUNCTION Segment<Dim> Edge(int i, int j) const
  {
    return {vertices[i], vertices[j]};
  }
};

/**
 * The tag and dimension of a shape type.
 */
template <typename Shape>
struct GeometryTraits;

template <int Dim>
struct GeometryTraits<Point<Dim>>
{
  using tag = PointTag;
  static constexpr int dimension = Dim;
};

template <int Dim>
struct GeometryTraits<AABBox<Dim>>
{
  using tag = BoxTag;
  static constexpr int dimension = Dim;
};

template <int Dim>
struct GeometryTraits<Segment<Dim>>
{
  using tag = SegmentTag;
  static constexpr int dimension = Dim;
};

template <int Dim>
struct GeometryTraits<Simplex<Dim>>
{
  using tag = SimplexTag;
  static constexpr int dimension = Dim;
};

template <typename Shape>
using tag_t = typename GeometryTraits<Shape>::tag;

template <typename Shape>
inline constexpr int dimension_v = GeometryTraits<Shape>::dimension;

} // namespace pcms

#endif // PCMS_LOCALIZATION_SHAPES_HPP
