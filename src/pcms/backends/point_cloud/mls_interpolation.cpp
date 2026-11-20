#include <pcms/backends/point_cloud/mls_interpolation.hpp>
#include <pcms/backends/point_cloud/mls_interpolation_impl.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <cmath>

namespace pcms::field
{

// 'a' is a shape parameter
// the value of 'a' is higher if the data is localized
// the value of 'a' is smaller if the data is farther

// gaussian
struct RBF_GAUSSIAN
{
  double a;
  RBF_GAUSSIAN(double a_val) : a(a_val) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    PCMS_CHECK_PRINTF(rho_sq >= 0,
                      "ERROR: square of cutoff distance should always be "
                      "positive but the value is %.16f\n",
                      rho_sq);
    PCMS_CHECK_PRINTF(r_sq >= 0,
                      "ERROR: square of  distance should always be positive "
                      "but the value is %.16f\n",
                      r_sq);
    double r = Kokkos::sqrt(r_sq);
    if (rho_sq < r_sq) {
      return 0;
    } else {
      return Kokkos::exp(-a * a * r * r);
    }
  }
};

struct RBF_C4
{
  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    PCMS_CHECK_PRINTF(
      rho_sq > 0, "ERROR: rho_sq in rbf has to be positive, but got %.16f\n",
      rho_sq);
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    double limit = 1 - ratio;
    double phi;
    if (rho_sq < r_sq) {
      phi = 0;
    } else {
      phi = 5 * pow(ratio, 5) + 30 * pow(ratio, 4) + 72 * pow(ratio, 3) +
            82 * pow(ratio, 2) + 36 * ratio + 6;
      phi = phi * pow(limit, 6);
    }
    PCMS_CHECK_PRINTF(!std::isnan(phi),
                      "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                      r_sq, rho_sq);
    return phi;
  }
};

struct RBF_CONST
{
  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    PCMS_CHECK_PRINTF(
      rho_sq > 0, "ERROR: rho_sq in rbf has to be positive, but got %.16f\n",
      rho_sq);
    double phi = (rho_sq < r_sq) ? 0.0 : 1.0;
    PCMS_CHECK_PRINTF(!std::isnan(phi),
                      "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                      r_sq, rho_sq);
    return phi;
  }
};

struct NoOp
{
  KOKKOS_INLINE_FUNCTION
  double operator()(double, double) const { return 1.0; }
};

struct RBF_MULTIQUADRIC
{
  double a;
  RBF_MULTIQUADRIC(double a_val) : a(a_val) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    double phi =
      (rho_sq < r_sq) ? 0.0 : Kokkos::sqrt(1.0 + a * a * ratio * ratio);
    PCMS_CHECK_PRINTF(!std::isnan(phi),
                      "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                      r_sq, rho_sq);
    return phi;
  }
};

// inverse multiquadric
struct RBF_INVMULTIQUADRIC
{
  double a;
  RBF_INVMULTIQUADRIC(double a_val) : a(a_val) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    double phi =
      (rho_sq < r_sq) ? 0.0 : 1.0 / Kokkos::sqrt(1.0 + a * a * ratio * ratio);
    PCMS_CHECK_PRINTF(!std::isnan(phi),
                      "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                      r_sq, rho_sq);
    return phi;
  }
};

struct RBF_THINPLATESPLINE
{
  double a;
  RBF_THINPLATESPLINE(double a_val) : a(a_val) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    double phi =
      (rho_sq < r_sq) ? 0.0 : a * a * ratio * ratio * Kokkos::log(a * ratio);
    PCMS_CHECK_PRINTF(!std::isnan(phi),
                      "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                      r_sq, rho_sq);
    return phi;
  }
};

struct RBF_CUBIC
{
  double a;
  RBF_CUBIC(double a_val) : a(a_val) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(double r_sq, double rho_sq) const
  {
    double r = Kokkos::sqrt(r_sq);
    double rho = Kokkos::sqrt(rho_sq);
    double ratio = r / rho;
    double phi = (rho_sq < r_sq) ? 0.0 : a * a * a * ratio * ratio * ratio;
    PCMS_CHECK_PRINTF(!std::isnan(phi),
                      "ERROR: phi in rbf is NaN. r_sq, rho_sq = (%f, %f)\n",
                      r_sq, rho_sq);
    return phi;
  }
};

} // namespace pcms::field

namespace pcms::detail
{

namespace
{

// Calls `f` with the radial basis function object of `bf`.
template <typename F>
auto WithRadialBasis(RadialBasisFunction bf, double decay_factor, F&& f)
{
  switch (bf) {
    case RadialBasisFunction::RBF_GAUSSIAN:
      return f(field::RBF_GAUSSIAN{decay_factor});
    case RadialBasisFunction::RBF_C4: return f(field::RBF_C4{});
    case RadialBasisFunction::RBF_CONST: return f(field::RBF_CONST{});
    case RadialBasisFunction::NO_OP: return f(field::NoOp{});
    case RadialBasisFunction::RBF_MULTIQUADRIC:
      return f(field::RBF_MULTIQUADRIC{decay_factor});
    case RadialBasisFunction::RBF_INVMULTIQUADRIC:
      return f(field::RBF_INVMULTIQUADRIC{decay_factor});
    case RadialBasisFunction::RBF_THINPLATESPLINE:
      return f(field::RBF_THINPLATESPLINE{decay_factor});
    case RadialBasisFunction::RBF_CUBIC:
      return f(field::RBF_CUBIC{decay_factor});
  }
  return f(field::NoOp{});
}

} // namespace

Kokkos::View<Real*, DeviceMemorySpace> mls_interpolation(
  Kokkos::View<const Real*, DeviceMemorySpace> source_values,
  Kokkos::View<const Real*, DeviceMemorySpace> source_coordinates,
  Kokkos::View<const Real*, DeviceMemorySpace> target_coordinates,
  const SupportResults& support, LO dim, LO degree, RadialBasisFunction bf,
  double lambda, double tol, double decay_factor)
{
  return WithRadialBasis(bf, decay_factor, [&](auto rbf) {
    return field::detail::mls_interpolation(source_values, source_coordinates,
                                            target_coordinates, support, dim,
                                            degree, rbf, lambda, tol);
  });
}

Kokkos::View<Real*, DeviceMemorySpace> mls_coefficients(
  Kokkos::View<const Real*, DeviceMemorySpace> source_coordinates,
  Kokkos::View<const Real*, DeviceMemorySpace> target_coordinates,
  const SupportResults& support, LO dim, LO degree, RadialBasisFunction bf,
  double lambda, double tol, double decay_factor)
{
  return WithRadialBasis(bf, decay_factor, [&](auto rbf) {
    return field::detail::mls_coefficients(source_coordinates,
                                           target_coordinates, support, dim,
                                           degree, rbf, lambda, tol);
  });
}

} // namespace pcms::detail

namespace pcms::field
{

namespace detail
{

void calculate_basis_slice_lengths(IntHostMatView& array)
{
  int degree = array.extent(0);
  int dim = array.extent(1);

  if (degree == 0)
    return;

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
} // namespace pcms::field
