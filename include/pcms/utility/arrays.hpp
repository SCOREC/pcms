#ifndef PCMS_UTILITY_ARRAYS_HPP
#define PCMS_UTILITY_ARRAYS_HPP
#include "mdspan/mdspan.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/array_traits.hpp"

namespace pcms
{

template <int Rank, typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>,
          typename IndexType = LO>
using View =
  Kokkos::mdspan<ElementType, Kokkos::dextents<IndexType, Rank>, LayoutPolicy,
                 detail::memory_space_accessor<
                   std::remove_reference_t<ElementType>, MemorySpace>>;

template <typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
using Rank1View = View<1, ElementType, MemorySpace, LayoutPolicy>;

template <typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
using Rank2View = View<2, ElementType, MemorySpace, LayoutPolicy>;

template <typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
using Rank3View = View<3, ElementType, MemorySpace, LayoutPolicy>;

template <typename MemorySpace>
using GlobalIDView = View<1, const GO, MemorySpace, Kokkos::layout_right, GO>;

template <typename T, typename... Properties>
auto MakeConstRank2View(const Kokkos::View<T**, Properties...>& view)
{
  using KokkosView = Kokkos::View<T**, Properties...>;
  using MemorySpace = typename KokkosView::memory_space;
  using LayoutPolicy = detail::layout_selector_t<KokkosView>;
  return Rank2View<const T, MemorySpace, LayoutPolicy>(
    view.data(), view.extent(0), view.extent(1));
}

// utility function to deep copy between layout incompatible views
// layout imcompatible view can't be deep_copied directly between different
// memory spaces, the workaround is provided from
// https://kokkos.org/kokkos-core-wiki/API/core/view/deep_copy.html#how-to-get-layout-incompatible-views-copied
template <typename DestView, typename SrcView>
auto DeepCopyMismatchLayouts(DestView& dest, const SrcView& src)
{
  static_assert(Kokkos::is_view<DestView>::value &&
                  Kokkos::is_view<SrcView>::value,
                "Both arguments must be Kokkos::View types");

  using DestMemSpace = typename DestView::memory_space;
  using SrcMemSpace = typename SrcView::memory_space;

  if constexpr (std::is_same_v<DestMemSpace, HostMemorySpace> &&
                !std::is_same_v<SrcMemSpace, HostMemorySpace>) {
    // Device to Host: mirror the device view to host
    auto src_tmp = Kokkos::create_mirror_view_and_copy(HostMemorySpace(), src);
    Kokkos::deep_copy(dest, src_tmp);
  } else if constexpr (!std::is_same_v<DestMemSpace, HostMemorySpace> &&
                       std::is_same_v<SrcMemSpace, HostMemorySpace>) {
    // Host to Device: create mirror on device from host view
    auto dest_tmp = Kokkos::create_mirror_view(dest);
    Kokkos::deep_copy(dest_tmp, src);
    Kokkos::deep_copy(dest, dest_tmp);
  } else {
    // Same memory space
    Kokkos::deep_copy(dest, src);
  }
}
/// Views over the storage of a Kokkos View with the layout of the Kokkos View:
/// MakeRank1View for a rank-1 Kokkos View, MakeRank2View for a rank-2 one
/// (const and non-const overloads), like MakeConstRank2View.
template <typename T, typename... Properties>
auto MakeRank1View(Kokkos::View<T*, Properties...>& view)
{
  using KokkosView = Kokkos::View<T*, Properties...>;
  using MemorySpace = typename KokkosView::memory_space;
  using LayoutPolicy = detail::layout_selector_t<KokkosView>;
  return Rank1View<T, MemorySpace, LayoutPolicy>(view.data(), view.extent(0));
}

template <typename T, typename... Properties>
auto MakeRank1View(const Kokkos::View<T*, Properties...>& view)
{
  using KokkosView = Kokkos::View<T*, Properties...>;
  using MemorySpace = typename KokkosView::memory_space;
  using LayoutPolicy = detail::layout_selector_t<KokkosView>;
  return Rank1View<const T, MemorySpace, LayoutPolicy>(view.data(),
                                                       view.extent(0));
}

template <typename T, typename... Properties>
auto MakeRank2View(Kokkos::View<T**, Properties...>& view)
{
  using KokkosView = Kokkos::View<T**, Properties...>;
  using MemorySpace = typename KokkosView::memory_space;
  using LayoutPolicy = detail::layout_selector_t<KokkosView>;
  return Rank2View<T, MemorySpace, LayoutPolicy>(view.data(), view.extent(0),
                                                 view.extent(1));
}

template <typename T, typename... Properties>
auto MakeRank2View(const Kokkos::View<T**, Properties...>& view)
{
  using KokkosView = Kokkos::View<T**, Properties...>;
  using MemorySpace = typename KokkosView::memory_space;
  using LayoutPolicy = detail::layout_selector_t<KokkosView>;
  return Rank2View<const T, MemorySpace, LayoutPolicy>(
    view.data(), view.extent(0), view.extent(1));
}

/// Copies a Rank2View into a Kokkos View of the same extents element by
/// element, whatever the layouts are, on the device of `dest`.
template <typename DestView, typename SrcView>
void ConvertMismatchLayoutView2D(DestView& dest, const SrcView& src)
{
  Kokkos::parallel_for(
    "ConvertMismatchLayoutView2D",
    Kokkos::RangePolicy<typename DestView::execution_space>(0, src.extent(0)),
    KOKKOS_LAMBDA(int i) {
      for (int j = 0; j < src.extent(1); ++j) {
        dest(i, j) = src(i, j);
      }
    });
}

} // namespace pcms

#include "pcms/utility/detail/arrays.hpp"

#endif // PCMS_UTILITY_ARRAYS_HPP
