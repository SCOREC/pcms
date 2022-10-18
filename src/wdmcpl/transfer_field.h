#ifndef WDM_COUPLING_TRANSFER_FIELD_H
#define WDM_COUPLING_TRANSFER_FIELD_H
#include <utility>
#include "wdmcpl/arrays.h"
#include "wdmcpl/field_evaluation_methods.h"
#include "wdmcpl/field.h"

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
  set_nodal_data(target_field, make_array_view(source_data));
}

template <typename SourceField, typename TargetField,
          typename EvaluationMethod = Lagrange<1>>
void interpolate_field(const SourceField& source_field,
                       TargetField& target_field, EvaluationMethod method = {})
{
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
  auto coordinates_view = make_array_view(coordinates);

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
template <typename SourceField, typename TargetField>
void transfer_field(const SourceField& source, TargetField& target,
                    FieldTransferMethod transfer_method,
                    FieldEvaluationMethod evaluation_method)
{
  switch (transfer_method) {
    case FieldTransferMethod::None: return;
    case FieldTransferMethod::Copy:
      if constexpr(std::is_same_v<SourceField, TargetField>) {
        copy_field(source, target);
      }
      else {
        std::cerr<<"Source field and destination field must have same type to copy!\n";
        std::abort();
      }
      return;
    case FieldTransferMethod::Interpolate:
      switch (evaluation_method) {
        case FieldEvaluationMethod::None:
          std::cerr << "Cannot interpolate field with no evaluation method!";
          std::terminate();
        case FieldEvaluationMethod::Lagrange1:
          interpolate_field(source, target, Lagrange<1>{});
          break;
        case FieldEvaluationMethod::NearestNeighbor:
          interpolate_field(source, target, NearestNeighbor{});
          break;
          // no default case for compiler error on missing cases
      }
      return;
      // no default case for compiler error on missing transfer method
  }
}

} // namespace wdmcpl

#endif // WDM_COUPLING_TRANSFER_FIELD_H
