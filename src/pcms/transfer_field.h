#ifndef PCMS_COUPLING_TRANSFER_FIELD_H
#define PCMS_COUPLING_TRANSFER_FIELD_H
#include <utility>
#include "pcms/arrays.h"
#include "pcms/field_evaluation_methods.h"
#include "pcms/field.h"
#include "pcms/profile.h"

namespace pcms
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
  PCMS_FUNCTION_TIMER;
  const auto source_data = get_nodal_data(source_field);
  set_nodal_data(target_field, make_array_view(source_data));
}

template <typename SourceField, typename TargetField,
          typename EvaluationMethod = Lagrange<1>>
void interpolate_field(const SourceField& source_field,
                       TargetField& target_field, EvaluationMethod method = {})
{
  PCMS_FUNCTION_TIMER;
  // TODO: same_topology
  // if (same_topology(source_field, target_field)) {
  //  copy_field(source_field,target_field);
  // return;
  //}

  // One downside of this option is that it requires the source has to implement
  // a get_nodal_coordinates function that wouldn't otherwise be needed
  using source_coordinate_type = typename decltype(get_nodal_coordinates(
    std::declval<SourceField>()))::value_type;
  using target_coordinate_type = typename decltype(get_nodal_coordinates(
    std::declval<TargetField>()))::value_type;

  static constexpr bool source_field_has_coordinate_system =
    detail::HasCoordinateSystem<source_coordinate_type>::value;
  static constexpr bool target_field_has_coordinate_system =
    detail::HasCoordinateSystem<target_coordinate_type>::value;
  static constexpr bool needs_coordinate_transform =
    (source_field_has_coordinate_system &&
     target_field_has_coordinate_system) &&
    !std::is_same_v<target_coordinate_type, source_coordinate_type>;
  // field value_types must either both have coordinate systems, or both not
  // have coordinate systems
  static_assert(source_field_has_coordinate_system ==
                  target_field_has_coordinate_system,
                "cannot mix data with and without coordinate systems");

  auto coordinates = get_nodal_coordinates(target_field);
  auto coordinates_view = make_const_array_view(coordinates);

  if constexpr (needs_coordinate_transform) {
    static_assert(!needs_coordinate_transform,
                  "coordinate transforms not finalized");
    // using target_coordinate_system =
    //   typename target_coordinate_type::coordinate_system;
    // auto transformed_coordinates =
    //   coordinate_transform<target_coordinate_system>(coordinates_view);
    // const auto data =
    //   evaluate(source_field, method,
    //   make_array_view(transformed_coordinates));
    // set_nodal_data(target_field, make_array_view(data));
  } else {
    const auto data = evaluate(source_field, method, coordinates_view);
    set_nodal_data(target_field, make_array_view(data));
  }
}
} // namespace pcms

#endif // PCMS_COUPLING_TRANSFER_FIELD_H
