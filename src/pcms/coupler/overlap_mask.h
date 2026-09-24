#ifndef PCMS_OVERLAP_MASK_H
#define PCMS_OVERLAP_MASK_H

#include "pcms/field/field_layout.h"
#include "pcms/utility/arrays.h"
#include <Omega_h_array.hpp>
#include <functional>

namespace pcms
{

struct OverlapMask
{
  Kokkos::View<bool*, HostMemorySpace> is_overlap_;
  std::function<int8_t(int, int)> in_overlap_func_;

  // Default constructor: all DOFs are in overlap
  OverlapMask(size_t size) : is_overlap_("overlap_info", size)
  {
    Kokkos::deep_copy(is_overlap_, true);
  }

  // Construct from a function that determines overlap based on classification
  OverlapMask(size_t size, std::function<int8_t(int, int)> in_overlap)
    : is_overlap_("overlap_info", size), in_overlap_func_(std::move(in_overlap))
  {
    Kokkos::deep_copy(is_overlap_, true);
  }

  // Construct from Omega_h host array
  OverlapMask(size_t size, Omega_h::HostRead<Omega_h::I8> is_overlap_host)
    : is_overlap_("overlap_info", size)
  {
    for (size_t i = 0; i < size; ++i) {
      is_overlap_[i] = static_cast<bool>(is_overlap_host[i]);
    }
  }

  // Construct from Omega_h device array
  OverlapMask(size_t size, Omega_h::Read<Omega_h::I8> is_overlap_device)
    : is_overlap_("overlap_info", size)
  {
    auto is_overlap_host = Omega_h::HostRead(is_overlap_device);
    for (size_t i = 0; i < size; ++i) {
      is_overlap_[i] = static_cast<bool>(is_overlap_host[i]);
    }
  }

  // Get the mask, evaluating the function if needed. The mask is owned-indexed
  // (one entry per owned DOF holder), matching the owned-indexed exchange plan.
  Rank1View<const bool, HostMemorySpace> GetMask(
    const FieldLayout& layout) const
  {
    if (in_overlap_func_) {
      auto class_dims = layout.GetDOFHolderClassificationDimensionsHost();
      auto class_ids = layout.GetDOFHolderClassificationIdsHost();
      auto owned_to_local = layout.GetOwnedToLocalHost();
      const size_t n = is_overlap_.extent(0);
      if (owned_to_local.size() == 0) {
        // Non-distributed layout: owned == local.
        for (size_t i = 0; i < n; ++i) {
          is_overlap_[i] =
            static_cast<bool>(in_overlap_func_(class_dims[i], class_ids[i]));
        }
      } else {
        for (size_t o = 0; o < n; ++o) {
          const LO local = owned_to_local(o);
          is_overlap_[o] = static_cast<bool>(
            in_overlap_func_(class_dims[local], class_ids[local]));
        }
      }
    }
    return make_const_array_view(is_overlap_);
  }
};

} // namespace pcms

#endif // PCMS_OVERLAP_MASK_H
