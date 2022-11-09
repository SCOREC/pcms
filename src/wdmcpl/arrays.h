#ifndef WDM_COUPLING_ARRAYS_H
#define WDM_COUPLING_ARRAYS_H
#include "wdmcpl/external/mdspan.hpp"
#include "wdmcpl/types.h"
#include "wdmcpl/coordinate.h"
#include "wdmcpl/external/span.h"

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
  detail::memory_space_accessor<ElementType, MemorySpace>>;

// default implementation of make_array_view
 template <typename T, typename MemorySpace= typename T::memory_space,
 typename ElementType = typename T::value_type> ScalarArrayView<ElementType,
 MemorySpace> make_array_view(const T& array) {
   using std::data;
   using std::size;
   return {data(array), size(array)};
 }

} // namespace wdmcpl

#endif // WDM_COUPLING_ARRAYS_H
