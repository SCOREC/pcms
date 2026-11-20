#ifndef PCMS_GEOMETRY_DETAIL_BOUND_POINTS_HPP
#define PCMS_GEOMETRY_DETAIL_BOUND_POINTS_HPP

#include <Kokkos_Core.hpp>
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>

namespace pcms::detail
{

/// @brief The geometric factors at bound points (BoundPoints::ComputeFactors):
///        the Jacobian J of the map from the geometry's parameter coordinates
///        xi to its coordinate system, and the metric g and volume factor in
///        the parameter coordinates; one row per query point.
class GeometricFactors
{
public:
  using JacobianView =
    View<3, const Real, DeviceMemorySpace, Kokkos::layout_stride>;

  /// @param num_points number of points
  /// @param jacobian [points or 1][system][parameter]: d x_j / d xi_i
  /// @param metric the metric and volume factor in the parameter coordinates
  ///        at the points
  GeometricFactors(LO num_points,
                   Kokkos::View<const Real***, DeviceMemorySpace> jacobian,
                   std::shared_ptr<const BoundCoordinateSystem> metric)
    : num_points_(num_points),
      jacobian_(std::move(jacobian)),
      metric_(std::move(metric))
  {
  }

  /// @return J(p, j, i) = d x_j / d xi_i
  [[nodiscard]] JacobianView Jacobian() const noexcept
  {
    using Extents = Kokkos::dextents<LO, 3>;
    const LO stride0 =
      jacobian_.extent(0) == 1 ? 0 : static_cast<LO>(jacobian_.stride(0));
    return JacobianView(
      jacobian_.data(),
      Kokkos::layout_stride::mapping<Extents>(
        Extents(num_points_, jacobian_.extent(1), jacobian_.extent(2)),
        std::array<LO, 3>{stride0, static_cast<LO>(jacobian_.stride(1)),
                          static_cast<LO>(jacobian_.stride(2))}));
  }

  /// @return g(p, a, b), the metric in the parameter coordinates
  [[nodiscard]] BoundCoordinateSystem::MetricView Metric() const noexcept
  {
    return metric_->Metric();
  }

  /// @return the physical volume per unit of parameter volume at each point
  [[nodiscard]] BoundCoordinateSystem::FactorView VolumeFactor() const noexcept
  {
    return metric_->VolumeFactor();
  }

  [[nodiscard]] LO NumPoints() const noexcept { return num_points_; }

private:
  LO num_points_;
  Kokkos::View<const Real***, DeviceMemorySpace> jacobian_;
  std::shared_ptr<const BoundCoordinateSystem> metric_;
};

/// @brief The geometric factors of a geometry's coordinate system at points
/// located in
///        it: J = I and the coordinate system's closed-form metric and volume
///        factor.
/// @throws pcms_error if the coordinate system has no metric
[[nodiscard]] GeometricFactors SystemFactors(
  const CoordinateSystem& system,
  Kokkos::View<const Real**, DeviceMemorySpace> points);

/// @brief Makes `storage` hold at least `rows` entries, reallocating it,
///        uninitialized, when it does not; the first `rows` are then the ones
///        in use.
/// @param label label of a reallocated view
template <typename T>
void ReserveRows(Kokkos::View<T*, DeviceMemorySpace>& storage,
                 const std::string& label, std::size_t rows)
{
  if (storage.extent(0) < rows) {
    storage = Kokkos::View<T*, DeviceMemorySpace>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, label), rows);
  }
}

/// @brief Makes `storage` hold at least `rows` rows of `cols` columns, as the
///        rank-1 overload.
///
/// The first rows of a taller view are a padded view on a LayoutLeft memory
/// space (the default on GPUs), which the conversions to mdspan do not
/// honour, so there a view is reused only for as many rows as it has.
template <typename T>
void ReserveRows(Kokkos::View<T**, DeviceMemorySpace>& storage,
                 const std::string& label, std::size_t rows, std::size_t cols)
{
  constexpr bool rows_contiguous =
    std::is_same_v<typename Kokkos::View<T**, DeviceMemorySpace>::array_layout,
                   Kokkos::LayoutRight>;
  if (storage.extent(0) < rows || storage.extent(1) != cols ||
      (!rows_contiguous && storage.extent(0) != rows)) {
    storage = Kokkos::View<T**, DeviceMemorySpace>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, label), rows, cols);
  }
}

} // namespace pcms::detail

#endif // PCMS_GEOMETRY_DETAIL_BOUND_POINTS_HPP
