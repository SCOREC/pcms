#ifndef PCMS_COUPLING_ARRAYS_H
#define PCMS_COUPLING_ARRAYS_H
#include "mdspan/mdspan.hpp"
#include "pcms/types.h"
#include "pcms/memory_spaces.h"

namespace pcms
{

namespace detail
{

template <typename ElementType, typename MemorySpace>
struct memory_space_accessor
  : public Kokkos::default_accessor<ElementType>
{
  using memory_space = MemorySpace;
};

} // namespace detail

template <typename ContainerType, typename ElementType, typename Extents,
          typename LayoutPolicy, typename AccessorPolicy>
auto make_mdspan(const ContainerType& /* unused */)
  -> Kokkos::mdspan<ElementType, Extents, LayoutPolicy,
                               AccessorPolicy>
{
  static_assert(detail::dependent_always_false<ContainerType>::type,
                "creating mdspan is not implemented for type");
}

// TODO make_mdspan

template <typename ElementType, typename MemorySpace>
using ScalarArrayView = Kokkos::mdspan<
  ElementType, Kokkos::dextents<LO, 1>,
  Kokkos::layout_right,
  detail::memory_space_accessor<std::remove_reference_t<ElementType>,
                                MemorySpace>>;

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
  -> ScalarArrayView<const ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return ScalarArrayView<const ElementType, MemorySpace>{data(array),
                                                         size(array)};
}
template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_array_view(T& array) -> ScalarArrayView<ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return ScalarArrayView<ElementType, MemorySpace>{data(array), size(array)};
}

template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
  typename ElementType = detail::element_type_t<T>>
auto make_const_array_view(T& array) -> ScalarArrayView<const ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return ScalarArrayView<const ElementType, MemorySpace>{data(array), size(array)};
}
} // namespace pcms
#endif // PCMS_COUPLING_ARRAYS_H
