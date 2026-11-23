#ifndef PCMS_LOCALIZATION_DETAIL_SHAPES_HPP
#define PCMS_LOCALIZATION_DETAIL_SHAPES_HPP

#include <Kokkos_Core.hpp>

#include "pcms/localization/shapes.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/types.hpp"

namespace pcms::detail
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

// The tag of a shape type.
template <typename Shape>
struct ShapeTraits;

template <int Dim>
struct ShapeTraits<Point<Dim>>
{
  using tag = PointTag;
};

template <int Dim>
struct ShapeTraits<AABBox<Dim>>
{
  using tag = BoxTag;
};

template <int Dim>
struct ShapeTraits<Segment<Dim>>
{
  using tag = SegmentTag;
};

template <int Dim>
struct ShapeTraits<Simplex<Dim>>
{
  using tag = SimplexTag;
};

} // namespace pcms::detail

#endif // PCMS_LOCALIZATION_DETAIL_SHAPES_HPP
