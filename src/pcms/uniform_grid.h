#ifndef PCMS_COUPLING_UNIFORM_GRID_H
#define PCMS_COUPLING_UNIFORM_GRID_H
#include "pcms/bounding_box.h"
#include "Omega_h_vector.hpp"
#include <numeric>
namespace pcms
{

template <unsigned dim = 2>
struct UniformGrid
{
  // Make private?
  std::array<Real, dim> edge_length;
  std::array<Real, dim> bot_left;
  std::array<LO, dim> divisions;

public:
  [[nodiscard]] LO GetNumCells() const {
    return std::accumulate(divisions.begin(), divisions.end(), 1, std::multiplies<LO>{});
  }
  /// return the grid cell ID that the input point is inside or closest to if
  /// the point lies outside
  // take the view as a template because it might be a subview type
  //template <typename T>
  //[[nodiscard]] KOKKOS_INLINE_FUNCTION LO ClosestCellID(const T& point) const
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO ClosestCellID(const Omega_h::Vector<dim>& point) const
  {
    std::array<Real, dim> distance_within_grid;

    for (int i = 0; i < dim; ++i) {
      distance_within_grid[i] = point[i] - bot_left[i];
    }

    std::array<LO, dim> indexes;

    for (auto& index : indexes) {
      index = -1;
    }

    for (int i = 0; i < dim; ++i) {
      auto index = static_cast<LO>(std::floor(distance_within_grid[i] * divisions[i]/edge_length[i]));
      indexes[i] = std::clamp(index, 0, divisions[i] - 1);
    }
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(indexes);
    return GetCellIndex(indexes);
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<dim> GetCellBBOX(LO idx) const
  {
    auto index = GetDimensionedIndex(idx);
    reverse(index);

    std::array<Real, dim> half_width, center;

    for (int i = 0; i < dim; ++i) {
      half_width[i] = edge_length[i] / divisions[i] / 2;
    }

    for (int i = 0; i < dim; ++i) {
      center[i] = (2.0 * index[i] + 1.0) * half_width[i] + bot_left[i];
    }

    return { .center = center, .half_width = half_width };
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION std::array<LO, dim> GetDimensionedIndex(LO idx) const
  {
    LO stride = 1;
    for (std::size_t i = 0; i < divisions.size() - 1; ++i) {
      stride *= divisions[i];
    }
    std::array<LO, dim> result;

    for (int i = 0; i < dim; ++i) {
      result[i] = idx / stride;
      idx -= result[i] * stride;
      stride /= divisions[i];
    }

    return result;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO GetCellIndex(std::array<LO, dim> dimensionedIndex) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(dimensionedIndex);

    LO idx = 0;
    LO stride = 1;

    for (int i = 0; i < dim; ++i) {
      idx += dimensionedIndex[i] * stride;
      stride *= divisions[i];
    }

    return idx;
  }

private:
  template <typename T, std::size_t N>
  KOKKOS_INLINE_FUNCTION static void reverse(std::array<T, N>& arr)
  {
    for (size_t i = 0, j = arr.size() - 1; i < j; ++i, --j) {
      auto temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }
};

using Uniform2DGrid = UniformGrid<>;

} // namespace pcms

#endif // PCMS_COUPLING_UNIFORM_GRID_H
