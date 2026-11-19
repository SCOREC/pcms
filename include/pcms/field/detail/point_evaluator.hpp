#ifndef PCMS_FIELD_DETAIL_POINT_EVALUATOR_HPP
#define PCMS_FIELD_DETAIL_POINT_EVALUATOR_HPP

#include <Kokkos_Core.hpp>
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"

namespace pcms::detail
{

/// @brief Sets every column of the listed rows of `values` to `value`.
/// @param values output view
/// @param rows indices of the rows to set
/// @param value value to write
template <typename T, typename LayoutPolicy>
void FillRows(Rank2View<T, DeviceMemorySpace, LayoutPolicy> values,
              Kokkos::View<const LO*, DeviceMemorySpace> rows, T value)
{
  const auto num_cols = static_cast<int>(values.extent(1));
  const auto fill = KOKKOS_LAMBDA(const LO k)
  {
    for (int c = 0; c < num_cols; ++c) {
      values(rows(k), c) = value;
    }
  };
  detail::ForEach("PointEvaluator::FillRows", static_cast<LO>(rows.extent(0)),
                  fill);
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_POINT_EVALUATOR_HPP
