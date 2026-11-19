#ifndef PCMS_UTILITY_DETAIL_UNIFORM_GRID_HPP
#define PCMS_UTILITY_DETAIL_UNIFORM_GRID_HPP

#include "pcms/utility/uniform_grid.hpp"

namespace pcms::detail
{

/**
 * \brief Divisions per axis giving roughly one element per grid cell.
 *
 */
template <unsigned dim>
[[nodiscard]] LO DivisionsPerAxisForMesh(LO nelems, LO min_divisions = 10)
{
  const auto per_axis = static_cast<LO>(std::ceil(
    std::pow(static_cast<double>(nelems), 1.0 / static_cast<double>(dim))));
  return std::max(per_axis, min_divisions);
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_UNIFORM_GRID_HPP
