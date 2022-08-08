#ifndef WDM_COUPLING_FIELD_H
#define WDM_COUPLING_FIELD_H
#include "wdmcpl/external/span.h"
#include "wdmcpl/types.h"
namespace wdmcpl
{

namespace detail
{
template <typename Field>
constexpr void check_field()
{
  static_assert(
    std::is_same_v<typename Field::CoordinateSystem,
                   typename Field::CoordinateArrayType::CoordinateSystem>);
}
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

template <typename ArrayType, typename DataType, typename ExeSpace>
class ScalarArray
{
public:
  ScalarArray(ArrayType array) : data_(std::move(array)) {}
  using ExecutionSpace = ExeSpace;
  // nonstd::span<DataType> GetSpan() { return
  // wdmcpl::GetSpan<ArrayType,DataType>(data_); }
  nonstd::span<const DataType> GetSpan() const
  {
    return wdmcpl::GetSpan<ArrayType, DataType>(data_);
  }

private:
  ArrayType data_;
};

// should we use kokkos views and template on execution space?
template <typename NodeHandle, typename Coordinate>
struct NodalCoordinates
{
  nonstd::span<NodeHandle> node_hanldes;
  // Strong type for coordinates should have span of double/arithmatic type
  // internally. This allows for more efficient conversion/copying
  nonstd::span<Coordinate> node_coordinates;
};

template <typename Field>
ScalarArray<typename Field::ArrayType, typename Field::DataType,
            typename Field::ExecutionSpace>
evaluate(const CoordinateArray<
         typename Field::CoordinateSystem, typename Field::ArrayType,
         typename Field::DataType, typename Field::ExecutionSpace>&)
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of evaluate");
};
template <typename Field>
void set(Field& field,
         const ScalarArray<typename Field::ArrayType, typename Field::DataType,
                           typename Field::ExecutionSpace>&)
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of set");
}
// TODO
template <typename Field>
typename Field::CoordinateArrayType get_nodal_coordinates(const Field&)
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of get_nodal_coordinates");
}

template <typename Field>
typename Field::NodeHandleArray get_node_handles(const Field&)
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of get_node_handles");
}
} // namespace wdmcpl

#endif // WDM_COUPLING_FIELD_H
