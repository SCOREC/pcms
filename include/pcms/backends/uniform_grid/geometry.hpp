#ifndef PCMS_BACKENDS_UNIFORM_GRID_GEOMETRY_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <memory>

namespace pcms::uniform_grid
{

/// The geometry of a 2D uniform grid: an affine map of the unit cells, located
/// in closed form (no search to build).
/// @param grid the grid
/// @param system the grid's coordinate system, or null for Cartesian
[[nodiscard]] std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<2> grid, std::shared_ptr<const CoordinateSystem> system);

/// The geometry of a 3D uniform grid.
/// @param grid the grid
/// @param system the grid's coordinate system, or null for Cartesian
[[nodiscard]] std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<3> grid, std::shared_ptr<const CoordinateSystem> system);

} // namespace pcms::uniform_grid

#endif // PCMS_BACKENDS_UNIFORM_GRID_GEOMETRY_HPP
