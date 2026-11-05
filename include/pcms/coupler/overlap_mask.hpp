#ifndef PCMS_COUPLER_OVERLAP_MASK_HPP
#define PCMS_COUPLER_OVERLAP_MASK_HPP

#include "pcms/field/field_layout.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <functional>

namespace pcms
{

struct OverlapMask
{
  Kokkos::View<bool*, HostMemorySpace> is_overlap_;
  std::function<int8_t(int, int)> in_overlap_func_;

  // Default constructor: all DOF holders are in overlap
  OverlapMask(size_t size) : is_overlap_("overlap_info", size)
  {
    Kokkos::deep_copy(is_overlap_, true);
  }

  // Construct from a function that determines overlap from a DOF holder's
  // classification (dimension, id). Holders without classification are passed
  // as (-1, -1).
  explicit OverlapMask(std::function<int8_t(int, int)> in_overlap)
    : in_overlap_func_(std::move(in_overlap))
  {
  }

  // Construct from a precomputed per-DOF-holder host mask (e.g. an MFEM
  // subdomain selected by element attribute). Indexed by local DOF holder.
  OverlapMask(size_t size, Kokkos::View<bool*, HostMemorySpace> is_overlap_host)
    : is_overlap_("overlap_info", size)
  {
    PCMS_ALWAYS_ASSERT(is_overlap_host.extent(0) == size);
    Kokkos::deep_copy(is_overlap_, is_overlap_host);
  }

  // Get the mask, evaluating the function against the layout if one was given
  Kokkos::View<const bool*, HostMemorySpace> GetMask(
    const FieldLayout& layout) const
  {
    if (!in_overlap_func_) {
      return is_overlap_;
    }
    auto class_dims = layout.GetDOFHolderClassificationDimensionsHost();
    auto class_ids = layout.GetDOFHolderClassificationIdsHost();
    Kokkos::View<bool*, HostMemorySpace> is_overlap("overlap_info",
                                                    class_dims.size());
    for (size_t i = 0; i < is_overlap.extent(0); ++i) {
      is_overlap[i] =
        static_cast<bool>(in_overlap_func_(class_dims[i], class_ids[i]));
    }
    return is_overlap;
  }
};

} // namespace pcms

#endif // PCMS_COUPLER_OVERLAP_MASK_HPP
