#include "pcms/transient/accelerator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace pcms::transient
{

Real InterfaceAccelerator::Norm(const std::vector<Real>& r)
{
  Real sum = 0.0;
  for (Real value : r)
    sum += value * value;
  return std::sqrt(sum);
}

AitkenRelaxation::AitkenRelaxation(Real omega0, Real omega_max)
  : omega0_(omega0), omega_(omega0), omega_max_(omega_max)
{
}

void AitkenRelaxation::BeginWindow()
{
  omega_ = omega0_;
  have_prev_ = false;
}

Real AitkenRelaxation::Update(const std::vector<Real>& x_in,
                              const std::vector<Real>& x_out,
                              std::vector<Real>& x_next)
{
  const std::size_t size = x_in.size();
  std::vector<Real> residual_values(size);
  for (std::size_t i = 0; i < size; ++i)
    residual_values[i] = x_out[i] - x_in[i];
  const Real residual = Norm(residual_values);

  if (have_prev_) {
    Real numerator = 0.0;
    Real denominator = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
      const Real residual_delta = residual_values[i] - r_prev_[i];
      numerator += r_prev_[i] * residual_delta;
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
    x_next[i] = x_in[i] + omega_ * residual_values[i];

  r_prev_ = std::move(residual_values);
  have_prev_ = true;
  return residual;
}

} // namespace pcms::transient
