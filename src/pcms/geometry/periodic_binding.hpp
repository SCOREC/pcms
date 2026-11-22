#ifndef PCMS_GEOMETRY_PERIODIC_BINDING_HPP
#define PCMS_GEOMETRY_PERIODIC_BINDING_HPP

// Private: the periodic geometry class and its binding.

#include "pcms/geometry/geometry.hpp"
#include "pcms/geometry/periodic_geometry.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms::detail
{

/// Locates points on the periodic interval: the phase (the wrapped position as
/// a fraction of the period, in [0, 1)), the cell and the position in the
/// cell. (A free function: CUDA does not allow device lambdas in members.)
inline void LocatePeriodic(Kokkos::View<const Real**, DeviceMemorySpace> points,
                           LO num_cells, Real origin, Real period,
                           Kokkos::View<Real*, DeviceMemorySpace> phase,
                           Kokkos::View<LO*, DeviceMemorySpace> cells,
                           Kokkos::View<Real*, DeviceMemorySpace> parametric)
{
  ForEach(
    "LocatePeriodic", static_cast<LO>(points.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      Real t = (points(p, 0) - origin) / period;
      t -= Kokkos::floor(t);
      if (t >= 1) {
        t = 0; // t was a rounding error below a whole period
      }
      const Real s = t * num_cells;
      LO c = static_cast<LO>(Kokkos::floor(s));
      c = c >= num_cells ? num_cells - 1 : c;
      phase(p) = t;
      cells(p) = c;
      parametric(p) = s - c;
    });
}

/// Periodic BoundPoints: the phase, cell and position in the cell of each
/// point. Every point is found.
class PeriodicBinding final : public BoundPoints
{
public:
  PeriodicBinding(std::shared_ptr<const Geometry> geometry,
                  const BindRequest& request, LO num_cells, Real origin,
                  Real period)
    : BoundPoints(std::move(geometry), request),
      num_cells_(num_cells),
      origin_(origin),
      period_(period)
  {
    Locate();
  }

  /// the wrapped position as a fraction of the period, in [0, 1)
  [[nodiscard]] Kokkos::View<const Real*, DeviceMemorySpace> Phase()
    const noexcept
  {
    return phase_;
  }
  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> CellIndices()
    const noexcept
  {
    return cells_;
  }
  /// the position in the cell, in [0, 1)
  [[nodiscard]] Kokkos::View<const Real*, DeviceMemorySpace> Parametric()
    const noexcept
  {
    return parametric_;
  }

protected:
  void RebindImpl() override { Locate(); }

private:
  void Locate()
  {
    const auto n = static_cast<size_t>(NumPoints());
    phase_ = Kokkos::View<Real*, DeviceMemorySpace>("periodic_phase", n);
    cells_ = Kokkos::View<LO*, DeviceMemorySpace>("periodic_cells", n);
    parametric_ =
      Kokkos::View<Real*, DeviceMemorySpace>("periodic_parametric", n);
    LocatePeriodic(LocatedPointValues(), num_cells_, origin_, period_, phase_,
                   cells_, parametric_);
  }

  LO num_cells_;
  Real origin_;
  Real period_;
  Kokkos::View<Real*, DeviceMemorySpace> phase_;
  Kokkos::View<LO*, DeviceMemorySpace> cells_;
  Kokkos::View<Real*, DeviceMemorySpace> parametric_;
};

} // namespace pcms::detail

namespace pcms
{

class PeriodicGeometry final : public Geometry
{
public:
  PeriodicGeometry(std::shared_ptr<const CoordinateSystem> system, LO num_cells,
                   Real period, Real origin);

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return system_;
  }
  [[nodiscard]] LO NumCells() const noexcept { return num_cells_; }
  [[nodiscard]] Real Period() const noexcept { return period_; }
  [[nodiscard]] Real Origin() const noexcept { return origin_; }

protected:
  [[nodiscard]] std::shared_ptr<BoundPoints> BindImpl(
    const BindRequest& request) const override;

  /// @brief Gauss-Legendre per cell, or the trapezoid rule of the period.
  [[nodiscard]] Quadrature CreateQuadrature(
    const pcms::detail::QuadratureRule& rule) const override;

private:
  std::shared_ptr<const CoordinateSystem> system_;
  LO num_cells_;
  Real period_;
  Real origin_;
};

} // namespace pcms

#endif // PCMS_GEOMETRY_PERIODIC_BINDING_HPP
