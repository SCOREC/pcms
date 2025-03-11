#ifndef PCMS_COUPLING_UNIFORM_GRID_H
#define PCMS_COUPLING_UNIFORM_GRID_H
#include "pcms/bounding_box.h"
#include "Omega_h_vector.hpp"
#include <iostream>
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
    std::transform(point.begin(), point.end(), bot_left.begin(),
                   distance_within_grid.begin(), std::minus<>());
    
    std::array<LO, dim> indexes;
    indexes.fill(-1);

    for (int i = 0; i < dim; ++i) {
      auto index = static_cast<LO>(std::floor(distance_within_grid[i] * divisions[i]/edge_length[i]));
      indexes[i] = std::clamp(index, 0, divisions[i] - 1);
    }
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    std::reverse(indexes.begin(), indexes.end());
    return GetCellIndex(indexes);
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<dim> GetCellBBOX(LO idx) const
  {
    auto [i, j] = GetDimensionedIndex(idx);
    std::array<Real, dim> half_width = {edge_length[0] / (2.0 * divisions[0]),
                                        edge_length[1] / (2.0 * divisions[1])};
    AABBox<dim> bbox{.center = {(2.0 * j + 1.0) * half_width[0]+bot_left[0],
                                (2.0 * i + 1.0) * half_width[1]+bot_left[1]},
                     .half_width = half_width};
    return bbox;
  }
  [[nodiscard]] KOKKOS_INLINE_FUNCTION std::array<LO, dim> GetDimensionedIndex(LO idx) const
  {
    return {idx / divisions[0], idx % divisions[0]};
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO GetCellIndex(std::array<LO, dim> dimensionedIndex) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    std::reverse(dimensionedIndex.begin(), dimensionedIndex.end());

    LO idx = 0;
    LO stride = 1;

    for (int i = 0; i < dim; ++i) {
      idx += dimensionedIndex[i] * stride;
      stride *= divisions[i];
    }

    return idx;
  }
};

using Uniform2DGrid = UniformGrid<>;

} // namespace pcms

#endif // PCMS_COUPLING_UNIFORM_GRID_H
