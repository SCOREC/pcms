#ifndef WDM_COUPLING_FIELD_H
#define WDM_COUPLING_FIELD_H
#include "wdmcpl/external/span.h"
#include "wdmcpl/types.h"
#include "wdmcpl/arrays.h"
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
         ScalarArrayView<const typename Field::DataType,
                           typename Field::ExecutionSpace>)
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
