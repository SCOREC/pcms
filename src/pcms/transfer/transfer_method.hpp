#ifndef PCMS_TRANSFER_METHOD_H
#define PCMS_TRANSFER_METHOD_H

#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/transfer/copy.h"
#include "pcms/transfer/interpolator.h"
#include "pcms/transfer/omega_h_conservative_projection.hpp"
#include "pcms/transfer/omega_h_control_variate_projection.hpp"
#include "pcms/transfer/omega_h_mc_rhs_integrator.hpp"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/utility/types.h"
#include <cstdint>
#include <memory>

namespace pcms
{

class FunctionSpace;

// A "transfer method" is a parameter-carrying recipe: it owns the
// method-specific parameters, and its
//
//   std::unique_ptr<TransferOperator<T>> Build(const FunctionSpace& source,
//                                               const FunctionSpace& target)
//                                               const
//
// member turns a space pair into a transfer operator. Selecting or comparing
// methods is then a one-line change at the AddTransfer call site, e.g.
//   AddTransfer(src, tgt,
//   method::ConservativeMonteCarlo{.samples_per_element = 64});
//
// Methods are aggregates (no polymorphic base) so their parameters can be set
// with designated initializers; they are matched by AddTransfer as a
// duck-typed concept rather than a base class. The value type T is enforced
// through Build's return type: a Real-only method will not compile against
// FunctionHandle<float>, and templated methods work for any supported T.
//
// The recipes live in `pcms::method` so their short names (e.g. `Copy`) do not
// collide with the transfer-operator types they wrap.
namespace method
{

// Identity copy; source and target must share a layout. Any supported T.
template <typename T>
struct Copy
{
  std::unique_ptr<TransferOperator<T>> Build(const FunctionSpace& source,
                                             const FunctionSpace& target) const
  {
    return std::make_unique<pcms::Copy<T>>(source, target);
  }
};

// Point interpolation: evaluate the source field at the target DOF sites. Any
// supported T.
template <typename T>
struct Interpolation
{
  OutOfBoundsPolicy policy{};

  std::unique_ptr<TransferOperator<T>> Build(const FunctionSpace& source,
                                             const FunctionSpace& target) const
  {
    return std::make_unique<Interpolator<T>>(source, target, policy);
  }
};

// Conservative Galerkin (L2) projection via mesh intersection. Real only.
struct ConservativeIntersection
{
  std::unique_ptr<TransferOperator<Real>> Build(
    const FunctionSpace& source, const FunctionSpace& target) const
  {
    return std::make_unique<OmegaHConservativeProjection>(source, target);
  }
};

// Conservative projection with a Monte-Carlo / control-variate RHS. Real only.
struct ConservativeMonteCarlo
{
  int samples_per_element = 32;
  MonteCarloSampling sampling = MonteCarloSampling::UniformRandom;
  std::uint64_t seed = 8675309;

  std::unique_ptr<TransferOperator<Real>> Build(
    const FunctionSpace& source, const FunctionSpace& target) const
  {
    return std::make_unique<OmegaHControlVariateProjection>(
      source, target, samples_per_element, sampling, seed);
  }
};

} // namespace method

} // namespace pcms

#endif // PCMS_TRANSFER_METHOD_H
