#include <cmath>
#include "mls_interpolation.hpp"

namespace pcms
{

// RBF_GAUSSIAN Functor
struct RBF_GAUSSIAN
{
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

    // 'a' is a spreading factor/decay factor
    // the value of 'a' is higher if the data is localized
    // the value of 'a' is smaller if the data is farther
    int a = 5;

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

struct NoOp
{
  OMEGA_H_INLINE
  double operator()(double, double) const { return 1.0; }
};

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, RadialBasisFunction bf)
{

  const auto nvertices_target = target_coordinates.size() / dim;

  Write<Real> interpolated_values(nvertices_target, 0,
                                  "approximated target values");
  switch (bf) {
    case RadialBasisFunction::RBF_GAUSSIAN:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, RBF_GAUSSIAN{});
      break;

    case RadialBasisFunction::RBF_C4:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, RBF_C4{});
      break;

    case RadialBasisFunction::RBF_CONST:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, RBF_CONST{});
      break;

    case RadialBasisFunction::NO_OP:
      interpolated_values = detail::mls_interpolation(
        source_values, source_coordinates, target_coordinates, support, dim,
        degree, NoOp{});
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

      int max_size;
      int min_size;
      if (nsupports > basis_size) {
        max_size = nsupports;
        min_size = basis_size;
      } else {
        max_size = basis_size;
        min_size = nsupports;
      }
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
