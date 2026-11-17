#ifndef PCMS_LOCALIZATION_UNIFORM_GRID_SEARCH_HPP
#define PCMS_LOCALIZATION_UNIFORM_GRID_SEARCH_HPP

#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <Kokkos_Core.hpp>
#include <Omega_h_vector.hpp>

namespace pcms
{

/// @brief The cell of a uniform grid containing each point, and the point's
///        position in that cell. Closed form: there is nothing to build.
/// @tparam Dim grid dimension
template <unsigned Dim>
class UniformGridPointSearch
{
public:
  struct Result
  {
    /// the cell containing the point, or OutsideId(nearest cell) for a point
    /// outside the grid
    LO cell;
    /// the position in the cell, each coordinate in [0, 1] for a point
    /// inside; for a point outside it is relative to the nearest cell and
    /// lies outside [0, 1]
    Omega_h::Vector<Dim> parametric;

    /// `cell` of a point outside the grid whose nearest cell is `id`.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION static constexpr LO OutsideId(LO id)
    {
      return -id - 1;
    }

    /// True when the point is outside the grid.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION constexpr bool IsOutside() const
    {
      return cell < 0;
    }

    /// The cell containing the point, or the nearest one for a point outside
    /// the grid.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION constexpr LO CellId() const
    {
      return cell < 0 ? -cell - 1 : cell;
    }
  };

  explicit UniformGridPointSearch(UniformGrid<Dim> grid) : grid_(grid) {}

  /// Per point: the cell containing it and its position in the cell, each
  /// coordinate in [0, 1]; a point outside the grid gets OutsideId(nearest
  /// cell) and the position relative to that cell (outside [0, 1]).
  /// @param points one row per point, `Dim` columns
  [[nodiscard]] Kokkos::View<Result*> operator()(
    Kokkos::View<const Real**, DeviceMemorySpace> points) const
  {
    const auto n = static_cast<LO>(points.extent(0));
    Kokkos::View<Result*> results(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "uniform_grid_search_results"),
      n);
    // the lambda must not capture this
    const auto grid = grid_;
    pcms::detail::ForEach(
      "UniformGridPointSearch", n, KOKKOS_LAMBDA(const LO i) {
        Omega_h::Vector<Dim> point;
        for (unsigned d = 0; d < Dim; ++d) {
          point[d] = points(i, d);
        }
        const LO cell = grid.ClosestCellID(point);
        const auto bbox = grid.GetCellBBOX(cell);
        Result result;
        result.cell =
          grid.IsPointInBounds(point) ? cell : Result::OutsideId(cell);
        for (unsigned d = 0; d < Dim; ++d) {
          const Real cell_min = bbox.center[d] - bbox.half_width[d];
          const Real cell_max = bbox.center[d] + bbox.half_width[d];
          result.parametric[d] = (point[d] - cell_min) / (cell_max - cell_min);
        }
        results(i) = result;
      });
    return results;
  }

private:
  UniformGrid<Dim> grid_;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_UNIFORM_GRID_SEARCH_HPP
