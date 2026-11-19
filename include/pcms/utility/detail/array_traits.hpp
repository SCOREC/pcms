#ifndef PCMS_UTILITY_DETAIL_ARRAY_TRAITS_HPP
#define PCMS_UTILITY_DETAIL_ARRAY_TRAITS_HPP

#include <Kokkos_Core.hpp>
#include "mdspan/mdspan.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <array>
#include <type_traits>

namespace pcms::detail
{

template <typename ElementType, typename MemorySpace>
struct memory_space_accessor : public Kokkos::default_accessor<ElementType>
{
  using memory_space = MemorySpace;

  // Default constructor
  memory_space_accessor() = default;

  // Converting constructor to allow adding const qualification
  template <typename OtherElementType,
            typename = std::enable_if_t<
              std::is_convertible_v<OtherElementType (*)[], ElementType (*)[]>>>
  constexpr memory_space_accessor(
    const memory_space_accessor<OtherElementType, MemorySpace>&) noexcept
  {
  }
};
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

template <typename T, size_t N>
struct arr_trait<Kokkos::Array<T, N>>
{
  using type = T;
};

template <typename T>
using element_type_t = typename arr_trait<T, HasValueType<T>::value>::type;

// Layout selector: detects layout from source type
template <typename T, typename = std::void_t<>>
struct layout_selector
{
  using type = Kokkos::layout_right; // Default for host arrays and std::array
};

// Specialization for types with array_layout (Kokkos Views)
template <typename T>
struct layout_selector<T, std::void_t<typename T::array_layout>>
{
  using type = std::conditional_t<
    std::is_same_v<typename T::array_layout, Kokkos::LayoutLeft>,
    Kokkos::layout_left, Kokkos::layout_right>;
};

template <typename T>
using layout_selector_t = typename layout_selector<T>::type;

// Default layout selector for memory spaces
// GPU memory spaces typically use LayoutLeft for better coalescing
template <typename MemorySpace, typename = std::void_t<>>
struct default_layout_for_memory_space
{
  using type = Kokkos::layout_right; // Default for host
};

#ifdef PCMS_HAS_DISTINCT_DEVICE_MEMORY_SPACE
// Specialization for device memory space - use LayoutLeft for GPU
template <>
struct default_layout_for_memory_space<DeviceMemorySpace>
{
  using type = Kokkos::layout_left;
};
#endif

template <typename MemorySpace>
using default_layout_for_memory_space_t =
  typename default_layout_for_memory_space<MemorySpace>::type;

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_ARRAY_TRAITS_HPP
