#ifndef PCMS_BACKENDS_UNIFORM_GRID_GEOMETRY_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <memory>

namespace pcms::uniform_grid
{

/// @brief The geometry of a uniform grid
/// @param grid the grid
/// @param system the grid's coordinate system, or null for Cartesian
/// @throws pcms_error if the grid is not inside the fundamental interval of
///         a periodic axis of `system`
template <unsigned Dim>
[[nodiscard]] std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<Dim> grid, std::shared_ptr<const CoordinateSystem> system);

} // namespace pcms::uniform_grid

#endif // PCMS_BACKENDS_UNIFORM_GRID_GEOMETRY_HPP
