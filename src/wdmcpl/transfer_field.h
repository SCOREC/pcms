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

template <int o>
struct Lagrange
{
  static constexpr int order{o};
};
// variable order lagrange interpolation
template <>
struct Lagrange<0>
{
  Lagrange(int order) : order{order} {};
  int order;
};

template <typename SourceField, typename TargetField,
          typename EvaluationMethod = Lagrange<1>>
void project_field(const SourceField& source_field, TargetField& target_field,
                    EvaluationMethod method = {})
{
  throw;
//  auto coordinates = target_field.get_nodal_coordinates();
//  auto transformed_coordinates =
//coordinate_transform<TargetField::coordinate_system>(std::move(coordinates));
//  auto node_handles = target_field.get_node_handles();
//  auto data = source_field.evaluate(method, transformed_coordinates);
//  // coordinate transform data if necessary
//  target_field.set(node_handles, data);
}
} // namespace wdmcpl

#endif // WDM_COUPLING_TRANSFER_FIELD_H
