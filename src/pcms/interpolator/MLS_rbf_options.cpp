#include <cmath>
#include "MLS_rbf_options.hpp"
#include "MLSInterpolation.hpp"
#include <Omega_h_fail.hpp>

// TODO add PCMS namespace to all interpolation code

// RBF_GAUSSIAN Functor
struct RBF_GAUSSIAN
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    OMEGA_H_CHECK_PRINTF(rho_sq > 0,
                         "ERROR: square of cutoff distance should always be "
                         "positive but the value is %.16f\n",
                         rho_sq);

    OMEGA_H_CHECK_PRINTF(r_sq > 0,
                         "ERROR: square of  distance should always be positive "
                         "but the value is %.16f\n",
                         r_sq);

    // 'a' is a spreading factor/decay factor
    // the value of 'a' is higher if the data is localized
    // the value of 'a' is smaller if the data is farther
    int a = 20;

    double r = sqrt(r_sq);
    double rho = sqrt(rho_sq);
    double ratio = r / rho;
    double limit = 1 - ratio;

    if (limit < 0) {
      phi = 0;

    } else {
      phi = exp(-a * a * r * r);
    }

    return phi;
  }
};

// RBF_C4 Functor
struct RBF_C4
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    double r = sqrt(r_sq);
    OMEGA_H_CHECK_PRINTF(
      rho_sq > 0, "ERROR: rho_sq in rbf has to be positive, but got %.16f\n",
      rho_sq);
    double rho = sqrt(rho_sq);
    double ratio = r / rho;
    double limit = 1 - ratio;
    if (limit < 0) {
      phi = 0;

    } else {
      phi = 5 * pow(ratio, 5) + 30 * pow(ratio, 4) + 72 * pow(ratio, 3) +
            82 * pow(ratio, 2) + 36 * ratio + 6;
      phi = phi * pow(limit, 6);
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);
    return phi;
  }
};

// RBF_const Functor
//
struct RBF_CONST
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    double r = sqrt(r_sq);
    OMEGA_H_CHECK_PRINTF(
      rho_sq > 0, "ERROR: rho_sq in rbf has to be positive, but got %.16f\n",
      rho_sq);
    double rho = sqrt(rho_sq);
    double ratio = r / rho;
    double limit = 1 - ratio;
    if (limit < 0) {
      phi = 0;

    } else {
      phi = 1.0;
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);
    return phi;
  }
};

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, Write<Real> radii2,
                              RadialBasisFunction bf)
{

  const auto nvertices_target = target_coordinates.size() / dim;

  Write<Real> interpolated_values(nvertices_target, 0,
                                  "approximated target values");
  switch (bf) {
    case RadialBasisFunction::RBF_GAUSSIAN:
      interpolated_values =
        mls_interpolation(source_values, source_coordinates, target_coordinates,
                          support, dim, degree, radii2, RBF_GAUSSIAN{});
      break;

    case RadialBasisFunction::RBF_C4:
      interpolated_values =
        mls_interpolation(source_values, source_coordinates, target_coordinates,
                          support, dim, degree, radii2, RBF_C4{});
      break;

    case RadialBasisFunction::RBF_CONST:
      interpolated_values =
        mls_interpolation(source_values, source_coordinates, target_coordinates,
                          support, dim, degree, radii2, RBF_CONST{});
      break;
  }

  return interpolated_values;
}
