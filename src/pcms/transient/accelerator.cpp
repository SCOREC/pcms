#include "pcms/transient/accelerator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace pcms::transient
{

Real InterfaceAccelerator::Norm(const std::vector<Real>& residual)
{
  Real squared_norm = 0.0;
  for (Real value : residual)
    squared_norm += value * value;
  return std::sqrt(squared_norm);
}

AitkenRelaxation::AitkenRelaxation(Real omega0, Real omega_max)
  : omega0_(omega0), omega_(omega0), omega_max_(omega_max)
{
}

void AitkenRelaxation::BeginWindow()
{
  omega_ = omega0_;
  have_previous_residual_ = false;
}

Real AitkenRelaxation::Update(const std::vector<Real>& x_in,
                              const std::vector<Real>& x_out,
                              std::vector<Real>& x_next)
{
  const std::size_t size = x_in.size();
  std::vector<Real> residual(size);
  for (std::size_t i = 0; i < size; ++i)
    residual[i] = x_out[i] - x_in[i];

  const Real residual_norm = Norm(residual);

  if (have_previous_residual_) {
    Real numerator = 0.0;
    Real denominator = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
      const Real residual_delta = residual[i] - previous_residual_[i];
      numerator += previous_residual_[i] * residual_delta;
      denominator += residual_delta * residual_delta;
    }

    if (denominator > 0.0) {
      omega_ = -omega_ * numerator / denominator;
      if (!std::isfinite(omega_) || omega_ == 0.0)
        omega_ = omega0_;
      omega_ = std::clamp(omega_, -omega_max_, omega_max_);
    }
  }

  x_next.resize(size);
  for (std::size_t i = 0; i < size; ++i)
    x_next[i] = x_in[i] + omega_ * residual[i];

  previous_residual_ = std::move(residual);
  have_previous_residual_ = true;
  return residual_norm;
}

} // namespace pcms::transient
