#ifndef PCMS_COUPLING_ARRAYS_H
#define PCMS_COUPLING_ARRAYS_H
#include "mdspan/mdspan.hpp"
#include "pcms/utility/types.h"
#include "pcms/utility/memory_spaces.h"

namespace pcms
{

namespace detail
{

template <typename ElementType, typename MemorySpace>
struct memory_space_accessor : public Kokkos::default_accessor<ElementType>
{
  using memory_space = MemorySpace;
};

} // namespace detail

template <typename ContainerType, typename ElementType, typename Extents,
          typename LayoutPolicy, typename AccessorPolicy>
auto make_mdspan(const ContainerType& /* unused */)
  -> Kokkos::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>
{
  static_assert(detail::dependent_always_false<ContainerType>::type,
                "creating mdspan is not implemented for type");
}

// TODO make_mdspan

template <int Rank, typename ElementType, typename MemorySpace,
          typename IndexType = LO>
using View = Kokkos::mdspan<
  ElementType, Kokkos::dextents<IndexType, Rank>, Kokkos::layout_right,
  detail::memory_space_accessor<std::remove_reference_t<ElementType>,
                                MemorySpace>>;

template <typename ElementType, typename MemorySpace>
using Rank1View = View<1, ElementType, MemorySpace>;

template <typename ElementType, typename MemorySpace>
using Rank2View = View<2, ElementType, MemorySpace>;

template <typename ElementType, typename MemorySpace>
using Rank3View = View<3, ElementType, MemorySpace>;

template <typename ElementType, typename MemorySpace>
using Rank4View = View<4, ElementType, MemorySpace>;

template <typename MemorySpace>
using GlobalIDView = View<1, const GO, MemorySpace, GO>;

namespace detail
{
template <typename T, typename = std::void_t<>>
struct HasValueType : std::false_type
{};
template <typename T>
struct HasValueType<T, std::void_t<typename T::value_type>> : std::true_type
{};

template <typename T, typename = std::void_t<>>
struct memory_space_selector
{
  using type = HostMemorySpace;
};
template <typename T>
struct memory_space_selector<T, std::void_t<typename T::memory_space>>
{
  using type = typename T::memory_space;
};

template <typename T>
using memory_space_selector_t = typename memory_space_selector<T>::type;

template <typename T, bool = false>
struct arr_trait;
template <typename T>
struct arr_trait<T, true>
{
  using type = typename T::value_type;
};

template <typename T, size_t N>
struct arr_trait<T[N]>
{
  using type = T;
};

template <typename T, size_t N>
struct arr_trait<std::array<T, N>>
{
  using type = T;
};

template <typename T>
using element_type_t = typename arr_trait<T, HasValueType<T>::value>::type;

} // namespace detail
  // default implementation of make_array_view
template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_array_view(const T& array)
  -> Rank1View<const ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return Rank1View<const ElementType, MemorySpace>{data(array), size(array)};
}
template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_array_view(T& array) -> Rank1View<ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return Rank1View<ElementType, MemorySpace>{data(array), size(array)};
}

template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_const_array_view(T& array)
  -> Rank1View<const ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return Rank1View<const ElementType, MemorySpace>{data(array), size(array)};
}

// utility function to deep copy between layout incompatible views
// layout imcompatible view can't be deep_copied directly between different
// memory spaces, the workaround is provided from
// https://kokkos.org/kokkos-core-wiki/API/core/view/deep_copy.html#how-to-get-layout-incompatible-views-copied
template <typename DestView, typename SrcView>
auto deep_copy_mismatch_layouts(DestView& dest, const SrcView& src)
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

// utility function to fill a view with sequentially increasing values
template <typename T>
void iota_view(Kokkos::View<T*> view, T start = 0)
{
  Kokkos::parallel_for(
    "iota_view", view.extent(0), KOKKOS_LAMBDA(LO i) { view[i] = start + i; });
}
} // namespace pcms
#endif // PCMS_COUPLING_ARRAYS_H
