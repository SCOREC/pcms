#ifndef PCMS_TRANSIENT_ACCELERATOR_HPP
#define PCMS_TRANSIENT_ACCELERATOR_HPP

#include "pcms/utility/types.h"
#include <vector>

namespace pcms::transient
{

// Updates interface values during a fixed-point iteration.
class InterfaceAccelerator
{
public:
  // Reset iteration history before solving a new coupling window.
  virtual void BeginWindow() = 0;

  // Compute the next interface iterate from equally sized input and output
  // vectors. Return the residual norm used to test convergence.
  virtual Real Update(const std::vector<Real>& x_in,
                      const std::vector<Real>& x_out,
                      std::vector<Real>& x_next) = 0;

  virtual ~InterfaceAccelerator() = default;

protected:
  // Compute the local L2 norm of a residual vector. Distributed
  // implementations must reduce the squared sum across ranks.
  static Real Norm(const std::vector<Real>& r);
};

// Aitken relaxation with a factor updated from consecutive residuals.
class AitkenRelaxation : public InterfaceAccelerator
{
public:
  // omega0 is the first-iteration relaxation factor. omega_max limits the
  // magnitude of factors computed from later residuals.
  explicit AitkenRelaxation(Real omega0 = 0.5, Real omega_max = 1e6);

  // Restore the initial factor and discard the previous residual.
  void BeginWindow() override;

  // Form r = x_out - x_in, update the relaxation factor when history exists,
  // and write x_next = x_in + omega*r.
  Real Update(const std::vector<Real>& x_in, const std::vector<Real>& x_out,
              std::vector<Real>& x_next) override;

private:
  Real omega0_;
  Real omega_;
  Real omega_max_;
  std::vector<Real> r_prev_;
  bool have_prev_ = false;
};

} // namespace pcms::transient

#endif // PCMS_TRANSIENT_ACCELERATOR_HPP
