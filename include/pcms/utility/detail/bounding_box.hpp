#ifndef PCMS_UTILITY_DETAIL_BOUNDING_BOX_HPP
#define PCMS_UTILITY_DETAIL_BOUNDING_BOX_HPP

#include "pcms/utility/bounding_box.hpp"

namespace pcms::detail
{

template <int dim>
KOKKOS_INLINE_FUNCTION bool intersects(const AABBox<dim>& a,
                                       const AABBox<dim>& b)
{
  for (int i = 0; i < dim; ++i) {
    if (std::abs(a.center[i] - b.center[i]) >
        (a.half_width[i] + b.half_width[i]))
      return false;
  }
  return true;
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_BOUNDING_BOX_HPP
