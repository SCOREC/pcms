#include <pcms/interpolator/mls_interpolation.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <cmath>

namespace pcms
{

// RBF_GAUSSIAN Functor
struct RBF_GAUSSIAN
{
  // 'a' is a spreading factor/decay factor
  // the value of 'a' is higher if the data is localized
  // the value of 'a' is smaller if the data is farther

  double a;

  RBF_GAUSSIAN(double a_val) : a(a_val) {}

  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    OMEGA_H_CHECK_PRINTF(rho_sq >= 0,
                         "ERROR: square of cutoff distance should always be "
                         "positive but the value is %.16f\n",
                         rho_sq);

    OMEGA_H_CHECK_PRINTF(r_sq >= 0,
                         "ERROR: square of  distance should always be positive "
                         "but the value is %.16f\n",
                         r_sq);

    double r = Kokkos::sqrt(r_sq);
    if (rho_sq < r_sq) {
      phi = 0;

    } else {
      phi = Kokkos::exp(-a * a * r * r);
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
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    double limit = 1 - ratio;

    OMEGA_H_CHECK_PRINTF(
      rho_sq > 0, "ERROR: rho_sq in rbf has to be positive, but got %.16f\n",
      rho_sq);
    if (rho_sq < r_sq) {
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
    OMEGA_H_CHECK_PRINTF(
      rho_sq > 0, "ERROR: rho_sq in rbf has to be positive, but got %.16f\n",
      rho_sq);

    if (rho_sq < r_sq) {
      phi = 0.0;

    } else {
      phi = 1.0;
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);
    return phi;
  }
};

struct NoOp
{
  OMEGA_H_INLINE
  double operator()(double, double) const { return 1.0; }
};

struct RBF_MULTIQUADRIC
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    if (rho_sq < r_sq) {
      phi = 0.0;
    } else {
      phi = Kokkos::sqrt(1.0 + ratio * ratio);
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);

    return phi;
  }
};

struct RBF_INVMULTIQUADRIC
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho if (rho_sq < r_sq)
    {
      phi = 0.0;
    }
    else
    {
      phi = 1.0 / Kokkos::sqrt(1.0 + ratio * ratio);
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);

    return phi;
  }
};

struct RBF_THINPLATESPLINE
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho if (rho_sq < r_sq)
    {
      phi = 0.0;
    }
    else
    {
      phi = r * r * Kokkos::log(ratio);
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);

    return phi;
  }
};

struct RBF_CUBIC
{
  OMEGA_H_INLINE
  double operator()(double r_sq, double rho_sq) const
  {
    double phi;
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho double limit = 1 - ratio if (rho_sq < r_sq)
    {
      phi = 0.0;
    }
    else
    {
      phi = limit * limit * limit;
    }

    OMEGA_H_CHECK_PRINTF(!std::isnan(phi),
                         "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                         r_sq, rho_sq);

    return phi;
  }
};

Omega_h::Write<Omega_h::Real> mls_interpolation(
  const Omega_h::Reals source_values, const Omega_h::Reals source_coordinates,
  const Omega_h::Reals target_coordinates, const SupportResults& support,
  const Omega_h::LO& dim, const Omega_h::LO& degree, RadialBasisFunction bf,
  double lambda, double tol, double decay_factor)
{

  const auto nvertices_target = target_coordinates.size() / dim;

  Omega_h::Write<Omega_h::Real> interpolated_values(
    nvertices_target, 0, "approximated target values");
  switch (bf) {
    case RadialBasisFunction::RBF_GAUSSIAN:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, RBF_GAUSSIAN{decay_factor}, lambda, tol);
      break;

    case RadialBasisFunction::RBF_C4:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, RBF_C4{}, lambda, tol);
      break;

    case RadialBasisFunction::RBF_CONST:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, RBF_CONST{}, lambda, tol);
      break;

    case RadialBasisFunction::NO_OP:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, NoOp{}, lambda, tol);
      break;
  }

  return interpolated_values;
}

namespace detail
{

void calculate_basis_slice_lengths(IntHostMatView& array)
{
  int degree = array.extent(0);
  int dim = array.extent(1);

  for (int j = 0; j < dim; ++j) {
    array(0, j) = 1;
  }

  for (int i = 0; i < degree; ++i) {
    array(i, 0) = 1;
  }

  for (int i = 1; i < degree; ++i) {
    for (int j = 1; j < dim; ++j) {
      array(i, j) = array(i, j - 1) + array(i - 1, j);
    }
  }
}

int calculate_basis_vector_size(const IntHostMatView& array)
{
  int sum = 1;
  int degree = array.extent(0);
  int dim = array.extent(1);

  for (int i = 0; i < degree; ++i) {
    for (int j = 0; j < dim; ++j) {
      sum += array(i, j);
    }
  }

  return sum;
}

int calculate_scratch_shared_size(const SupportResults& support,
                                  const int nvertices_target, int basis_size,
                                  int dim)
{

  IntDeviceVecView shmem_each_team("stores the size required for each team",
                                   nvertices_target);
  Kokkos::parallel_for(
    "calculate the size required for scratch for each team", nvertices_target,
    KOKKOS_LAMBDA(const int i) {
      int start_ptr = support.supports_ptr[i];
      int end_ptr = support.supports_ptr[i + 1];
      int nsupports = end_ptr - start_ptr;

      size_t total_shared_size = 0;
      total_shared_size +=
        ScratchMatView::shmem_size(basis_size, basis_size); // Vt
      total_shared_size +=
        ScratchMatView::shmem_size(basis_size, nsupports) * 2; // temp_matrix
      total_shared_size +=
        ScratchMatView::shmem_size(nsupports, basis_size); // vandermonde matrix
      total_shared_size +=
        ScratchVecView::shmem_size(basis_size) *
        3; // sigma, target_basis_vector, solution_coefficients
      total_shared_size += ScratchVecView::shmem_size(nsupports) *
                           3; // work, phi_vector, support_values
      total_shared_size +=
        ScratchMatView::shmem_size(nsupports, dim); // local_source_points
      total_shared_size +=
        ScratchMatView::shmem_size(nsupports, nsupports) * 2; // U, Ut
      shmem_each_team(i) = total_shared_size;
    });

  int shared_size = 0;
  Kokkos::parallel_reduce(
    "find_max", nvertices_target,
    KOKKOS_LAMBDA(const int i, int& max_val_temp) {
      if (shmem_each_team(i) > max_val_temp) {
        max_val_temp = shmem_each_team(i);
      }
    },
    Kokkos::Max<int>(shared_size));

  return shared_size;
}

} // namespace detail
} // namespace pcms
