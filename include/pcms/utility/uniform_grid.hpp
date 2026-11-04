#ifndef PCMS_UTILITY_UNIFORM_GRID_HPP
#define PCMS_UTILITY_UNIFORM_GRID_HPP
#include "pcms/utility/bounding_box.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <array>
#include <cmath>
namespace pcms
{

template <unsigned dim>
struct UniformGrid
{
  // Make private?
  Kokkos::Array<Real, dim> edge_length;
  Kokkos::Array<Real, dim> bot_left;
  Kokkos::Array<LO, dim> divisions;

public:
  [[nodiscard]] LO GetNumCells() const
  {
    LO total = 1;
    for (std::size_t i = 0; i < dim; ++i)
      total *= divisions[i];
    return total;
  }
  /// return the grid cell ID that the input point is inside or closest to if
  /// the point lies outside
  /// `point` is any type indexable with point[i] for i < dim (e.g.
  /// Kokkos::Array<Real, dim>, a subview, Omega_h::Vector<dim>).
  template <typename Point>
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO ClosestCellID(const Point& point) const
  {
    Kokkos::Array<LO, dim> indexes;

    for (size_t i = 0; i < dim; ++i) {
      indexes[i] = AxisCellIndex(i, point[i]);
    }
    return CellIndexFromAxisIndices(indexes);
  }

  /// Index of the cell containing `coord` along `axis`, clamped to the grid.
  /// The single definition of the grid's binning rule: keeping the point query
  /// and the cell-range query on this one formula is what guarantees a point
  /// always lands in a cell that its own element was binned into.
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO AxisCellIndex(size_t axis,
                                                        Real coord) const
  {
    const auto index = static_cast<LO>(std::floor(
      (coord - bot_left[axis]) * divisions[axis] / edge_length[axis]));
    return Kokkos::clamp(index, 0, divisions[axis] - 1);
  }

  /// Linear cell ID from per-axis indices given in coordinate order (x,y,z).
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  CellIndexFromAxisIndices(Kokkos::Array<LO, dim> indexes) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(indexes);
    return GetCellIndex(indexes);
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<dim> GetCellBBOX(LO idx) const
  {
    auto index = GetDimensionedIndex(idx);
    reverse(index);

    Kokkos::Array<Real, dim> half_width, center;

    for (size_t i = 0; i < dim; ++i) {
      half_width[i] = edge_length[i] / divisions[i] / 2;
    }

    for (size_t i = 0; i < dim; ++i) {
      center[i] = (2.0 * index[i] + 1.0) * half_width[i] + bot_left[i];
    }

    return {.center = center, .half_width = half_width};
  }

  /// Check if a point is within the grid bounds. `point` is any type
  /// indexable with point[i] for i < dim.
  template <typename Point>
  [[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsPointInBounds(
    const Point& point) const
  {
    for (size_t i = 0; i < dim; ++i) {
      Real coord = point[i];
      Real grid_min = bot_left[i];
      Real grid_max = bot_left[i] + edge_length[i];
      if (coord < grid_min || coord > grid_max) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<LO, dim>
  GetDimensionedIndex(LO idx) const
  {
    LO stride = 1;
    for (std::size_t i = 0; i < dim - 1; ++i) {
      stride *= divisions[i];
    }
    Kokkos::Array<LO, dim> result;

    for (size_t i = 0; i < dim; ++i) {
      result[i] = idx / stride;
      idx -= result[i] * stride;
      stride /= divisions[i];
    }

    return result;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  GetCellIndex(Kokkos::Array<LO, dim> dimensionedIndex) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(dimensionedIndex);

    LO idx = 0;
    LO stride = 1;

    for (size_t i = 0; i < dim; ++i) {
      idx += dimensionedIndex[i] * stride;
      stride *= divisions[i];
    }

    return idx;
  }

private:
  template <typename T, std::size_t N>
  KOKKOS_INLINE_FUNCTION static void reverse(Kokkos::Array<T, N>& arr)
  {
    for (size_t i = 0, j = N - 1; i < j; ++i, --j) {
      auto temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }
};

using Uniform2DGrid = UniformGrid<2>;
using Uniform3DGrid = UniformGrid<3>;

/**
 * \brief Divisions per axis giving roughly one element per grid cell.
 *
 */
template <unsigned dim>
LO DivisionsPerAxisForMesh(LO nelems, LO min_divisions = 10)
{
  const auto per_axis = static_cast<LO>(std::ceil(
    std::pow(static_cast<double>(nelems), 1.0 / static_cast<double>(dim))));
  return std::max(per_axis, min_divisions);
}

} // namespace pcms

#endif // PCMS_UTILITY_UNIFORM_GRID_HPP
