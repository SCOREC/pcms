#ifndef PCMS_FIELD_POINT_ROWS_HPP
#define PCMS_FIELD_POINT_ROWS_HPP

// Private: the evaluation of a field at a set of points as sparse rows of DOF
// weights.

#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>

namespace pcms::detail
{

/// @brief The evaluation at a set of points as CSR rows: row p holds the
///        (DOF holder, weight) entries of point p, the same for every
///        component. A holder indexes the layout's DOF holders, ghosts
///        included, in the order of FieldLayout::GetOwnedHost.
struct PointRows
{
  Kokkos::View<const LO*, DeviceMemorySpace> offsets; // [rows + 1]
  Kokkos::View<const LO*, DeviceMemorySpace> holders;
  Kokkos::View<const Real*, DeviceMemorySpace> values;
};

/// @brief An evaluator whose evaluation is a PointRows.
class RowsSource
{
public:
  [[nodiscard]] virtual const PointRows& Rows() const = 0;

  /// @return the global id of each DOF holder the rows index; a ghost has
  ///         its owner's id
  [[nodiscard]] virtual GlobalIDView<DeviceMemorySpace> HolderGids() const = 0;

protected:
  ~RowsSource() = default;
};

/// @return `evaluator` as a RowsSource, or null if its evaluation has no rows
[[nodiscard]] inline const RowsSource* AsRowsSource(
  const PointEvaluator<Real>& evaluator)
{
  return dynamic_cast<const RowsSource*>(&evaluator);
}

} // namespace pcms::detail

#endif // PCMS_FIELD_POINT_ROWS_HPP
