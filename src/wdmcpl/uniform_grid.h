#ifndef WDM_COUPLING_UNIFORM_GRID_H
#define WDM_COUPLING_UNIFORM_GRID_H
#include "wdmcpl/bounding_box.h"
#include "Omega_h_vector.hpp"
#include <iostream>
#include <numeric>
namespace wdmcpl
{
struct UniformGrid
{
  // Make private?
  static constexpr int dim = 2;
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
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO ClosestCellID(const Omega_h::Vector<2>& point) const
  {
    std::array<Real, dim> distance_within_grid{point[0] - bot_left[0],
                                               point[1] - bot_left[1]};
    std::array<LO, dim> indexes{-1, -1};
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    for (int i = 0; i < dim; ++i) {
      if (distance_within_grid[i] <= 0) {
        indexes[dim - (i + 1)] = 0;
      } else if (distance_within_grid[i] >= edge_length[i]) {
        indexes[dim - (i + 1)] = divisions[i] - 1;
      } else {
        indexes[dim - (i + 1)] = 
          static_cast<LO>(std::floor(distance_within_grid[i] * divisions[i]/edge_length[i]));
      }
    }
    return GetCellIndex(indexes[0], indexes[1]);
  }
  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<dim> GetCellBBOX(LO idx) const
  {
    auto [i, j] = GetTwoDCellIndex(idx);
    std::array<Real, dim> half_width = {edge_length[0] / (2.0 * divisions[0]),
                                        edge_length[1] / (2.0 * divisions[1])};
    AABBox<dim> bbox{.center = {(2.0 * j + 1.0) * half_width[0]+bot_left[0],
                                (2.0 * i + 1.0) * half_width[1]+bot_left[0]},
                     .half_width = half_width};
    return bbox;
  }
  [[nodiscard]] KOKKOS_INLINE_FUNCTION std::array<LO, 2> GetTwoDCellIndex(LO idx) const
  {
    return {idx / divisions[0], idx % divisions[0]};
  }
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO GetCellIndex(LO i, LO j) const
  {
    OMEGA_H_CHECK(i >= 0 && j >= 0 && i < divisions[1] && j < divisions[0]);
    return i * divisions[0] + j;
  }
};
} // namespace wdmcpl

#endif // WDM_COUPLING_UNIFORM_GRID_H
