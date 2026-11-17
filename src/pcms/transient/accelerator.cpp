#include "pcms/transient/accelerator.hpp"

#include "pcms/utility/assert.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace pcms
{

// The two sums that Aitken relaxation needs from the residual history.
struct AitkenSums
{
  Real numerator = 0.0;   // dot(r_{k-1}, r_k - r_{k-1})
  Real denominator = 0.0; // dot(r_k - r_{k-1}, r_k - r_{k-1})

  KOKKOS_INLINE_FUNCTION AitkenSums& operator+=(const AitkenSums& other)
  {
    numerator += other.numerator;
    denominator += other.denominator;
    return *this;
  }
};

} // namespace pcms

template <>
struct Kokkos::reduction_identity<pcms::AitkenSums>
{
  KOKKOS_FORCEINLINE_FUNCTION static pcms::AitkenSums sum()
  {
    return pcms::AitkenSums{};
  }
};

namespace pcms
{

// Sums of the Aitken update.
AitkenSums ReduceAitkenSums(
  Rank1View<const Real, DeviceMemorySpace> x_in,
  Rank1View<const Real, DeviceMemorySpace> x_out,
  Kokkos::View<const Real*, DeviceMemorySpace> previous)
{
  AitkenSums sums;
  Kokkos::parallel_reduce(
    "AitkenRelaxation::ReduceAitkenSums",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, x_in.extent(0)),
    KOKKOS_LAMBDA(const LO i, AitkenSums& local) {
      const Real delta = (x_out(i) - x_in(i)) - previous(i);
      local.numerator += previous(i) * delta;
      local.denominator += delta * delta;
    },
    Kokkos::Sum<AitkenSums, Kokkos::HostSpace>(sums));
  return sums;
}

// Writes the next state, and returns the squared norm of the residual.
Real ApplyAitkenUpdate(Rank1View<const Real, DeviceMemorySpace> x_in,
                       Rank1View<const Real, DeviceMemorySpace> x_out,
                       Rank1View<Real, DeviceMemorySpace> x_next,
                       Kokkos::View<Real*, DeviceMemorySpace> residual_out,
                       Real omega)
{
  Real residual_squared = 0.0;
  Kokkos::parallel_reduce(
    "AitkenRelaxation::ApplyAitkenUpdate",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, x_in.extent(0)),
    KOKKOS_LAMBDA(const LO i, Real& local) {
      const Real in = x_in(i);
      const Real residual = x_out(i) - in;
      x_next(i) = in + omega * residual;
      residual_out(i) = residual;
      local += residual * residual;
    },
    residual_squared);
  return residual_squared;
}

AitkenRelaxation::AitkenRelaxation(Real omega0, Real omega_max)
  : omega0_(omega0), omega_(omega0), omega_max_(omega_max)
{
}

void AitkenRelaxation::BeginWindow()
{
  omega_ = omega0_;
  have_previous_residual_ = false;
  last_status_ = UpdateStatus::Initial;
}

Real AitkenRelaxation::Update(Rank1View<const Real, DeviceMemorySpace> x_in,
                              Rank1View<const Real, DeviceMemorySpace> x_out,
                              Rank1View<Real, DeviceMemorySpace> x_next)
{
  const std::size_t size = x_in.extent(0);
  if (x_out.extent(0) != size || x_next.extent(0) != size) {
    throw pcms_error(
      "AitkenRelaxation::Update: x_in, x_out and x_next must have the same "
      "size");
  }
  if (have_previous_residual_ && previous_residual_.extent(0) != size) {
    throw pcms_error(
      "AitkenRelaxation::Update: the size changed within a window, call "
      "BeginWindow() first");
  }

  last_status_ = UpdateStatus::Initial;
  if (have_previous_residual_) {
    const auto sums = ReduceAitkenSums(x_in, x_out, previous_residual_);
    last_status_ = UpdateStatus::Fallback;
    if (sums.denominator > 0.0) {
      omega_ = -omega_ * sums.numerator / sums.denominator;
      last_status_ = UpdateStatus::Updated;

      // A zero or non-finite omega is useless, so restart from omega0.
      if (!std::isfinite(omega_) || omega_ == 0.0) {
        omega_ = omega0_;
        last_status_ = UpdateStatus::Fallback;
      }

      // Limit |omega| to omega_max.
      const Real clamped = std::clamp(omega_, -omega_max_, omega_max_);
      if (clamped != omega_ && last_status_ == UpdateStatus::Updated)
        last_status_ = UpdateStatus::Limited;
      omega_ = clamped;
    }
  } else if (previous_residual_.extent(0) != size) {
    // First iteration of the window: omega is omega0.
    previous_residual_ = Kokkos::View<Real*, DeviceMemorySpace>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "AitkenRelaxation::previous_residual"),
      size);
  }
  const Real residual_squared =
    ApplyAitkenUpdate(x_in, x_out, x_next, previous_residual_, omega_);
  have_previous_residual_ = true;

  return std::sqrt(residual_squared);
}
} // namespace pcms
