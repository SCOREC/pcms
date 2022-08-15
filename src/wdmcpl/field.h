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

/*
// should we use kokkos views and template on execution space?
template <typename NodeHandle, typename Coordinate>
struct NodalCoordinates
{
  nonstd::span<NodeHandle> node_hanldes;
  // Strong type for coordinates should have span of double/arithmatic type
  // internally. This allows for more efficient conversion/copying
  nonstd::span<Coordinate> node_coordinates;
};


// note value_type should/could have a coordinate system associated with it
template <typename Field>
auto evaluate(ScalarArrayView<const typename Field::value_type,
                              typename Field::ExecutionSpace>)
  -> ScalarArray<typename Field::ArrayType, typename Field::ExecutionSpace>
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of evaluate");
}
template <typename Field>
auto set(Field& field, ScalarArrayView<const typename Field::value_type,
                                       typename Field::execution_space>) -> void
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of set");
}
// TODO
template <typename Field>
auto get_nodal_coordinates(const Field&) -> typename Field::CoordinateArrayType
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of get_nodal_coordinates");
}

template <typename Field>
auto get_node_handles(const Field&) -> typename Field::NodeHandleArray
{
  static_assert(detail::dependent_always_false<Field>::value,
                "No Matching implementation of get_node_handles");
}
 */


} // namespace wdmcpl

#endif // WDM_COUPLING_FIELD_H
