#ifndef WDM_COUPLING_ARRAYS_H
#define WDM_COUPLING_ARRAYS_H
#include "wdmcpl/external/mdspan.hpp"
#include "wdmcpl/types.h"
#include "wdmcpl/coordinate.h"
#include "wdmcpl/external/span.h"
#include "wdmcpl/memory_spaces.h"

namespace wdmcpl
{

namespace detail
{

template <typename ElementType, typename MemorySpace>
struct memory_space_accessor
  : public std::experimental::default_accessor<ElementType>
{
  using memory_space = MemorySpace;
};

template <typename ElementType, typename CoordinateSystem, typename MemorySpace>
struct coordinate_accessor
  : public memory_space_accessor<ElementType, MemorySpace>
{
  using coordinate_system = CoordinateSystem;
  using memory_space = MemorySpace;
};

} // namespace detail

template <typename ArrayType>
nonstd::span<typename ArrayType::value_type> GetSpan(ArrayType /* unused */)
{
  static_assert(detail::dependent_always_false<ArrayType>::type,
                "creating span is not implemented for type");
}

template <typename ContainerType, typename ElementType, typename Extents,
          typename LayoutPolicy, typename AccessorPolicy>
auto make_mdspan(const ContainerType& /* unused */)
  -> std::experimental::mdspan<ElementType, Extents, LayoutPolicy,
                               AccessorPolicy>
{
  static_assert(detail::dependent_always_false<ContainerType>::type,
                "creating mdspan is not implemented for type");
}

template <typename CoordinateSystemT, typename ArrayType, typename DataType,
          typename ExeSpace>
class CoordinateArray
{
public:
  CoordinateArray(ArrayType array) : data_(std::move(array)) {}
  using coordinate_system = CoordinateSystemT;
  using execution_space = ExeSpace;
  using ArrayType::const_pointer;
  using ArrayType::const_reference;
  using ArrayType::pointer;
  using ArrayType::reference;
  using ArrayType::value_type;

  // TODO convert this to mdspan
  auto GetSpan() const
  {
    // return make_mdspan
    return wdmcpl::GetSpan<ArrayType, DataType>(data_);
  }

private:
  ArrayType data_;
};

// Array wrapper
// TODO convert instances of this to mdarray
template <typename Container, typename ExeSpace>
class ScalarArray
{
public:
  ScalarArray(Container array) : data_(std::move(array)) {}
  using execution_space = ExeSpace;
  using value_type = typename Container::value_type;
  using pointer = typename Container::pointer;
  using reference = typename Container::reference;
  using const_pointer = typename Container::const_pointer;
  using const_reference = typename Container::const_reference;
  // nonstd::span<DataType> GetSpan() { return
  // wdmcpl::GetSpan<Container,DataType>(data_); }
  nonstd::span<const value_type> GetSpan() const
  {
    return wdmcpl::GetSpan<Container>(data_);
  }

private:
  Container data_;
};

template <typename ElementType, typename CoordinateSystem,
          typename LayoutPolicy = std::experimental::layout_right,
          typename Container =
            std::vector<CoordinateElement<CoordinateSystem, ElementType>>,
          size_t N = 3>

using CoordinateMDArray = std::experimental::mdarray<
  CoordinateElement<CoordinateSystem, ElementType>,
  std::experimental::extents<int, std::experimental::dynamic_extent, N>,
  LayoutPolicy, Container>;

template <typename ElementType,
          typename LayoutPolicy = std::experimental::layout_right,
          typename Container = std::vector<ElementType>>
using MDArray = std::experimental::mdarray<
  ElementType,
  std::experimental::extents<int, std::experimental::dynamic_extent, 1>,
  LayoutPolicy, Container>;

template <typename ElementType, typename CoordinateSystem, typename MemorySpace,
          size_t N = 1>
using CoordinateArrayView = std::experimental::mdspan<
  CoordinateElement<ElementType, CoordinateSystem>,
  std::experimental::extents<LO, std::experimental::dynamic_extent, N>,
  std::experimental::layout_right,
  detail::memory_space_accessor<ElementType, MemorySpace>>;

// TODO make_mdspan

template <typename ElementType, typename MemorySpace>
using ScalarArrayView = std::experimental::mdspan<
  ElementType, std::experimental::dextents<LO, 1>,
  std::experimental::layout_right,
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
} // namespace wdmcpl
#endif // WDM_COUPLING_ARRAYS_H
