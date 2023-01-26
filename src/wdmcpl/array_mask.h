#ifndef WDM_COUPLING_ARRAY_MASK_H
#define WDM_COUPLING_ARRAY_MASK_H
#include "wdmcpl/arrays.h"
#include <Kokkos_Core.hpp>
namespace wdmcpl
{

namespace detail
{

template <typename MemorySpace>
struct ComputeMaskAV
{
  explicit ComputeMaskAV(ScalarArrayView<LO, MemorySpace> index_mask,
                         ScalarArrayView<const int8_t, MemorySpace> mask)
    : mask_(mask), index_mask_(index_mask)
  {
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(LO i, LO& update, bool final) const noexcept
  {
    update += (mask_(i) > 0);
    if (final) {
      index_mask_(i) = update;
    }
  }
  ScalarArrayView<LO, MemorySpace> index_mask_;
  ScalarArrayView<const int8_t, MemorySpace> mask_;
};
template <typename MemorySpace>
struct ScaleAV
{
  explicit ScaleAV(ScalarArrayView<LO, MemorySpace> arr,
                   ScalarArrayView<const int8_t, MemorySpace> s)
    : arr_(arr), s_(s)
  {
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const noexcept { arr_[i] *= s_[i]; }

  ScalarArrayView<LO, MemorySpace> arr_;
  ScalarArrayView<const int8_t, MemorySpace> s_;
};

} // namespace detail

// TODO replace mask/ filter_array with ArrayMask in Omega_h_field
template <typename MemorySpace>
class ArrayMask
{
public:
  using execution_space = typename MemorySpace::execution_space;
  ArrayMask() = default;
  // takes a mask where each entry is 1 for including the entry and 0 for
  // excluding the entry
  explicit ArrayMask(ScalarArrayView<const int8_t, MemorySpace> mask)
    : num_active_entries_(0)
  {
    // we use a parallel scan to construct the mask mapping so that filtering
    // can happen in parallel. This method gives us the index to fill into the
    // filtered array
    Kokkos::View<LO*, MemorySpace> index_mask("mask", mask.size());
    auto policy = Kokkos::RangePolicy<execution_space>(0, mask.size());
    auto index_mask_view = make_array_view(index_mask);
    Kokkos::parallel_scan(policy,detail::ComputeMaskAV<MemorySpace>{index_mask_view, mask},
                          num_active_entries_);
    //// set index mask to 0 anywhere that the original mask is 0
    Kokkos::parallel_for(policy, detail::ScaleAV<MemorySpace>{index_mask_view, mask});
    // does a shallow copy
    mask_ = index_mask;
  }
  template <typename T>
  auto Apply(ScalarArrayView<const T, MemorySpace> data,
             ScalarArrayView<T, MemorySpace> filtered_data,
             ScalarArrayView<const wdmcpl::LO, MemorySpace> permutation = {})
    const -> void
  {
    // it doesn't make sense to call this function when the mask is empty!
    WDMCPL_ALWAYS_ASSERT(!empty());
    WDMCPL_ALWAYS_ASSERT(data.size() == mask_.size());
    WDMCPL_ALWAYS_ASSERT(filtered_data.size() ==
                         static_cast<size_t>(num_active_entries_));

    auto policy = Kokkos::RangePolicy<execution_space>(0, mask_.size());
    // make local copy of the mask_ view to avoid problem with passing in "this"
    // ptr to the KOKKOS_LAMBDA
    auto mask = mask_;
    Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(LO i) {
        if (mask[i]) {
          const auto idx = mask[i] - 1;
          if (permutation.empty()) {
            filtered_data[idx] = data[i];
          } else {
            filtered_data[permutation[idx]] = data[i];
          }
        }
      });
  }
  template <typename T>
  [[nodiscard]] auto Apply(const ScalarArrayView<T, MemorySpace> data) const
    -> Kokkos::View<T, MemorySpace>
  {
    // it doesn't make sense to call this function when the mask is empty!
    WDMCPL_ALWAYS_ASSERT(!empty());
    WDMCPL_ALWAYS_ASSERT(data.size() == mask_.size());
    Kokkos::View<T, MemorySpace> filtered_data("filtered data",
                                               num_active_entries_);
    Apply(data, make_array_view(filtered_data));
    return filtered_data;
  }

  /// This function takes an array that has been filtered (size of
  /// num_active_entries) and sets the appropriate entries in output_array based
  /// on the active entries in the filter
  template <typename T>
  auto ToFullArray(
    ScalarArrayView<const T, MemorySpace> filtered_data,
    ScalarArrayView<T, MemorySpace> output_array,
    ScalarArrayView<const wdmcpl::LO, MemorySpace> permutation = {}) const
    -> void
  {
    if (empty()) {
      if (filtered_data.data_handle() != output_array.data_handle()) {
        WDMCPL_ALWAYS_ASSERT(output_array.size() == filtered_data.size());
        Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, filtered_data.size()),
          KOKKOS_LAMBDA(LO i) { output_array(i) = filtered_data(i); });
      }
    } else {

      WDMCPL_ALWAYS_ASSERT((LO)output_array.size() == mask_.size());
      WDMCPL_ALWAYS_ASSERT((LO)filtered_data.size() == num_active_entries_);
      REDEV_ALWAYS_ASSERT(filtered_data.size() == permutation.size() ||
                          permutation.empty());
      auto mask = mask_;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<execution_space>(0, mask_.size()),
        KOKKOS_LAMBDA(LO i) {
          if (mask[i]) {
            const auto idx = mask[i] - 1;
            output_array[i] = (!permutation.empty())
                                ? filtered_data[permutation[idx]]
                                : filtered_data[idx];
          }
        });
    }
  }
  [[nodiscard]] bool empty() const noexcept { return num_active_entries_ == 0; }
  // mask is true if the mask exists
  explicit operator bool() const noexcept { return !empty(); }

  [[nodiscard]] LO Size() const noexcept { return num_active_entries_; }
  // returns a view where each entry that's greater than 0 is active in the
  // filtered array and the value of that entry -1 is the local order of the
  // data
  [[nodiscard]] auto GetMap() const { return make_const_array_view(mask_); }

private:
  Kokkos::View<LO*, MemorySpace> mask_;
  LO num_active_entries_;
};
} // namespace wdmcpl

#endif // WDM_COUPLING_ARRAY_MASK_H
