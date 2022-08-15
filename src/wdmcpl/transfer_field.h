#ifndef WDM_COUPLING_TRANSFER_FIELD_H
#define WDM_COUPLING_TRANSFER_FIELD_H
#include <utility>
#include "wdmcpl/arrays.h"

namespace wdmcpl
{

/**
 * Pointwise copy of data in source_field to target field. Source and target
 * fields must have the same size and implicit field iteration order
 * @tparam Field
 * @param source_field
 * @param target_field
 */
template <typename Field>
void copy_field(const Field& source_field, Field& target_field)
{
  const auto source_data = get_nodal_data(source_field);
  set(target_field, make_array_view(source_data));
}

// Field evaluation methods
template <int o>
struct Lagrange
{
  static constexpr int order{o};
};
// variable order lagrange interpolation
template <>
struct Lagrange<0>
{
  explicit Lagrange(int order) : order{order} {};
  int order;
};

struct NearestNeighbor
{};

// end field evaluation methods

namespace detail
{
template <typename...>
using VoidT = void;

template <typename, typename = void>
struct HasCoordinateSystem : public std::false_type
{};

template <typename T>
struct HasCoordinateSystem<T, VoidT<typename T::coordinate_system>>
  : public std::true_type
{};

} // namespace detail

template <typename SourceField, typename TargetField,
          typename EvaluationMethod = Lagrange<1>>
void interpolate_field(const SourceField& source_field,
                       TargetField& target_field, EvaluationMethod method = {})
{
  // TODO: same_topology
  // if (same_topology(source_field, target_field)) {
  //  copy_field(source_field,target_field);
  //}

  // FIXME rather than getting coordinate system directly from Field get
  // coordinate system from result type of get_nodal_coordinates
  // this alleviates the need to have coordinate_system typedef in field class
  // proxy
  // One downside of this option is that it requires the source has to implement
  // a get_nodal_coordinates function that wouldn't otherwise be needed
  using source_coordinate_type = typename decltype(get_nodal_coordinates(
    std::declval<SourceField>()))::value_type;
  using target_coordinate_type = typename decltype(get_nodal_coordinates(
    std::declval<TargetField>()))::value_type;

  // constexpr bool source_field_has_coordinate_system =
  //   detail::HasCoordinateSystem<typename SourceField::value_type>::value;
  // constexpr bool target_field_has_coordinate_system =
  //   detail::HasCoordinateSystem<typename TargetField::value_type>::value;
  // constexpr bool needs_coordinate_transform =
  //   source_field_has_coordinate_system && target_field_has_coordinate_system;
  constexpr bool source_field_has_coordinate_system =
    detail::HasCoordinateSystem<source_coordinate_type>::value;
  constexpr bool target_field_has_coordinate_system =
    detail::HasCoordinateSystem<target_coordinate_type>::value;
  constexpr bool needs_coordinate_transform =
    source_field_has_coordinate_system && target_field_has_coordinate_system;
  // field value_types must either both have coordinate systems, or both not
  // have coordinate systems
  static_assert(source_field_has_coordinate_system ==
                  target_field_has_coordinate_system,
                "cannot mix data with and without coordinate systems");

  auto coordinates = get_nodal_coordinates(target_field);
  auto coordinates_view = make_array_view(coordinates);
  // TODO: filter operation if working on part of field?
  // auto node_handles = target_field.get_node_handles();

  if constexpr (needs_coordinate_transform) {
    using target_coordinate_system =
      typename target_coordinate_type::coordinate_system;
    auto transformed_coordinates =
      coordinate_transform<target_coordinate_system>(coordinates_view);
    const auto data =
      evaluate(source_field, method, make_array_view(transformed_coordinates));
    // coordinate transform data if necessary
    // set(target_field, node_handles, make_array_view(data));
    set(target_field, make_array_view(data));
  } else {
    const auto data = evaluate(source_field, method, coordinates);
    // coordinate transform data if necessary
    // set(target_field, node_handles, make_array_view(data));
    set(target_field, make_array_view(data));
  }
}
} // namespace wdmcpl

#endif // WDM_COUPLING_TRANSFER_FIELD_H
