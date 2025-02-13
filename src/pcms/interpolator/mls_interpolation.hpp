#ifndef MLS_RBF_OPTIONS_HPP
#define MLS_RBF_OPTIONS_HPP

#include "mls_interpolation_impl.hpp"

using namespace Omega_h;

namespace pcms
{
enum class RadialBasisFunction : LO
{
  RBF_GAUSSIAN = 0,
  RBF_C4,
  RBF_CONST

};

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, RadialBasisFunction bf);
} // namespace pcms
#endif
