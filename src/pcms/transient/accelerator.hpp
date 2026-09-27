#ifndef PCMS_TRANSIENT_ACCELERATOR_HPP
#define PCMS_TRANSIENT_ACCELERATOR_HPP

#include "pcms/utility/types.h"
#include <vector>

namespace pcms::transient
{

// Coupling interface values are updated here during a fixed-point iteration.
class InterfaceAccelerator
{
public:
  // Reset iteration history before solving a new coupling window.
  virtual void BeginWindow() = 0;

  // Compute the next interface state and return the residual norm.
  virtual Real Update(const std::vector<Real>& x_in,
                      const std::vector<Real>& x_out,
                      std::vector<Real>& x_next) = 0;

  virtual ~InterfaceAccelerator() = default;

protected:
  // Compute the L2 norm of the residual vector.
  static Real Norm(const std::vector<Real>& residual);
};

// Applies Aitken relaxation using residuals from consecutive iterations.
class AitkenRelaxation : public InterfaceAccelerator
{
public:
  // omega0 is the initial factor, omega_max is the maximum limit for omega.
  explicit AitkenRelaxation(Real omega0 = 0.5, Real omega_max = 1e6);

  void BeginWindow() override;

  Real Update(const std::vector<Real>& x_in, const std::vector<Real>& x_out,
              std::vector<Real>& x_next) override;

private:
  Real omega0_;
  Real omega_;
  Real omega_max_;
  std::vector<Real> previous_residual_;
  bool have_previous_residual_ = false;
};

} // namespace pcms::transient

#endif // PCMS_TRANSIENT_ACCELERATOR_HPP
