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

std::vector<int> get_partitions(int total_size, int num_partitions)
{
  std::vector<int> partitions(num_partitions + 1);
  int partition_size = total_size / num_partitions;
  for (int i = 0; i < num_partitions; ++i) {
    partitions[i] = i * partition_size;
  }
  partitions[num_partitions] = total_size;
  return partitions;
}

// TODO move this to Omega_h
#ifdef OMEGA_H_USE_KOKKOS
Reals get_subview(Reals const& input_read, std::size_t min, std::size_t max) {
  // if min max contains full range, return the original view
  if (min == 0 && max == input_read.size()) {
        return input_read;
  }
  if (min > max || max > input_read.size()) {
    throw std::out_of_range("Invalid range for subview: min or max out of bounds");
  }

  auto subview = Kokkos::subview(input_read.view(), Kokkos::make_pair(min, max));

  // for now, copy the subview to a new Reals object
  // TODO need a constructor for Reals that takes a Kokkos::Subview
  Write<Real> subview_copy(subview.size());
  Kokkos::parallel_for(
    "CopySubView", subview.size(),
    KOKKOS_LAMBDA(const int i) { subview_copy[i] = subview[i]; });

    Kokkos::fence();
    return subview_copy;
}
#endif

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, Write<Real> radii2,
                              RadialBasisFunction bf, bool partitioned, int num_partitions)
{
  assert(!(partitioned==true && num_partitions<=0));
  if (!partitioned) { num_partitions = 1; }
  const auto nvertices_target = target_coordinates.size() / dim;
  auto partition_indices = get_partitions(nvertices_target, num_partitions);

  Write<Real> interpolated_values(nvertices_target, 0,
                                  "approximated target values");
  for(int part=0; part < num_partitions; ++part) {

    auto min = partition_indices[part];
    auto max = partition_indices[part + 1];
    printf("!!Partition %d: min = %d, max = %d\n", part, min, max);
    auto partitioned_target_coordinates =
      get_subview(target_coordinates, min*dim, max*dim);
    Write<Real> partitioned_interpolated_values(max - min, 0,
                                        "partitioned approximated target values");


    switch (bf) {
      case RadialBasisFunction::RBF_GAUSSIAN:
        partitioned_interpolated_values = mls_interpolation(
          source_values, source_coordinates, partitioned_target_coordinates, support, dim,
          degree, radii2, RBF_GAUSSIAN{});
        break;

      case RadialBasisFunction::RBF_C4:
        partitioned_interpolated_values = mls_interpolation(
          source_values, source_coordinates, partitioned_target_coordinates, support, dim,
          degree, radii2, RBF_C4{});
        break;

      case RadialBasisFunction::RBF_CONST:
        partitioned_interpolated_values = mls_interpolation(
          source_values, source_coordinates, partitioned_target_coordinates, support, dim,
          degree, radii2, RBF_CONST{});
        break;
    }
    // copy the partitioned interpolated values to the full array
    Omega_h::parallel_for(
      "copy_partitioned_interpolated_values", max - min,
      OMEGA_H_LAMBDA(const int i) {
        interpolated_values[min + i] = partitioned_interpolated_values[i];
      });
  }

  return interpolated_values;
}
