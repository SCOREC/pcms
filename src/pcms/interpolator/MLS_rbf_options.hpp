#ifndef MLS_RBF_OPTIONS_HPP
#define MLS_RBF_OPTIONS_HPP

#include "adj_search.hpp"

using namespace Omega_h;

enum class RadialBasisFunction : LO
{
  RBF_GAUSSIAN = 0,
  RBF_C4,
  RBF_CONSTANT

};

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, Write<Real> radii2,
                              RadialBasisFunction bf);

#endif
