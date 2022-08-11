#ifndef WDM_COUPLING_ARRAYS_H
#define WDM_COUPLING_ARRAYS_H
#include "wdmcpl/external/mdspan.hpp"
#include "wdmcpl/types.h"

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


template <typename ArrayType, typename DataType>
nonstd::span<DataType> GetSpan(ArrayType array);

template <typename CoordinateSystemT, typename ArrayType, typename DataType,
  typename ExeSpace>
class CoordinateArray
{
public:
  CoordinateArray(ArrayType array) : data_(std::move(array)) {}
  using CoordinateSystem = CoordinateSystemT;
  using ExecutionSpace = ExeSpace;
  nonstd::span<const DataType> GetSpan() const
  {
    return wdmcpl::GetSpan<ArrayType, DataType>(data_);
  }

private:
  ArrayType data_;
};

template <typename Container, typename DataType, typename ExeSpace>
class ScalarArray
{
public:
  ScalarArray(Container array) : data_(std::move(array)) {}
  using ExecutionSpace = ExeSpace;
  // nonstd::span<DataType> GetSpan() { return
  // wdmcpl::GetSpan<Container,DataType>(data_); }
  nonstd::span<const DataType> GetSpan() const
  {
    return wdmcpl::GetSpan<Container, DataType>(data_);
  }

private:
  Container data_;
};



template <typename ElementType, typename CoordinateSystem, typename MemorySpace>
using CoordinateArrayView = std::experimental::mdspan<
  ElementType, std::experimental::dextents<LO, 1>,
  std::experimental::layout_right,
  detail::coordinate_accessor<ElementType, CoordinateSystem, MemorySpace>>;

template <typename ElementType, typename MemorySpace>
using ScalarArrayView = std::experimental::mdspan<
  ElementType, std::experimental::dextents<LO, 1>,
  std::experimental::layout_right,
  detail::memory_space_accessor<ElementType, MemorySpace>>;


} // namespace wdmcpl

#endif // WDM_COUPLING_ARRAYS_H
