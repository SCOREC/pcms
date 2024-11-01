#ifndef MLS_COEFFICIENTS_HPP
#define MLS_COEFFICIENTS_HPP

#include <cmath>
#include <Omega_h_fail.hpp>

#include "points.hpp"

#define PI_M 3.14159265358979323846

static constexpr int MAX_DIM = 3;

KOKKOS_INLINE_FUNCTION
double func(Coord& p)
{
  auto x = (p.x - 0.5) * PI_M * 2;
  auto y = (p.y - 0.5) * PI_M * 2;
  double Z = sin(x) * sin(y) + 2;
  return Z;
}

// computes the slice lengths of the of the polynomial basis
KOKKOS_INLINE_FUNCTION
void basisSliceLengths(Kokkos::View<int**, Kokkos::HostSpace>& array)
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

// finds the size of the polynomial basis vector
KOKKOS_INLINE_FUNCTION
int basisSize(const Kokkos::View<int**, Kokkos::HostSpace>& array)
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

// evaluates the polynomial basis
KOKKOS_INLINE_FUNCTION
void BasisPoly(ScratchVecView basis_monomial, const MatViewType& slice_length,
               Coord& p)
{
  basis_monomial(0) = 1;
  int dim = slice_length.extent(1);
  int degree = slice_length.extent(0);

  int prev_col = 0;
  int curr_col = 1;

  double point[MAX_DIM];
  point[0] = p.x;
  point[1] = p.y;

  if (dim == 3) {
    point[2] = p.z;
  }

  for (int i = 0; i < degree; ++i) {
    int offset = curr_col;
    for (int j = 0; j < dim; ++j) {
      for (int k = 0; k < slice_length(i, j); ++k) {
        basis_monomial(offset + k) = basis_monomial(prev_col + k) * point[j];
      }

      offset += slice_length(i, j);
    }

    prev_col = curr_col;
    curr_col = offset;
  }
}

KOKKOS_INLINE_FUNCTION
double rbf(double r_sq, double rho_sq)
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

// create vandermondeMatrix
KOKKOS_INLINE_FUNCTION
void VandermondeMatrix(ScratchMatView V,
                       const ScratchMatView local_source_points, int j,
                       const MatViewType slice_length)
{
  int N = local_source_points.extent(0);
  int dim = local_source_points.extent(1);

  Coord source_point;
  source_point.x = local_source_points(j, 0);
  source_point.y = local_source_points(j, 1);
  if (dim == 3) {
    source_point.z = local_source_points(j, 2);
  }
  ScratchVecView basis_monomial = Kokkos::subview(V, j, Kokkos::ALL());
  BasisPoly(basis_monomial, slice_length, source_point);
}

// moment matrix
KOKKOS_INLINE_FUNCTION
void PTphiMatrix(ScratchMatView pt_phi, ScratchMatView V, ScratchVecView Phi,
                 int j)
{
  int M = V.extent(0);
  int N = V.extent(1);

  ScratchVecView vandermonde_mat_row = Kokkos::subview(V, j, Kokkos::ALL());
  for (int k = 0; k < N; k++) {
    OMEGA_H_CHECK_PRINTF(!std::isnan(vandermonde_mat_row(k)),
                         "ERROR: vandermonde_mat_row is NaN for k = %d\n", k);
    OMEGA_H_CHECK_PRINTF(!std::isnan(Phi(j)),
                         "ERROR: Phi(j) in PTphiMatrix is NaN for j = %d\n", j);
    pt_phi(k, j) = vandermonde_mat_row(k) * Phi(j);
    OMEGA_H_CHECK_PRINTF(
      !std::isnan(pt_phi(k, j)),
      "ERROR: pt_phi in PTphiMatrix is NaN for k = %d, j = %d\n", k, j);
  }
}

// radial basis function vector
KOKKOS_INLINE_FUNCTION
void PhiVector(ScratchVecView Phi, const Coord target_point,
               const ScratchMatView local_source_points, int j,
               double cuttoff_dis_sq)
{
  int N = local_source_points.extent(0);
  double dx = target_point.x - local_source_points(j, 0);
  double dy = target_point.y - local_source_points(j, 1);
  double ds_sq = dx * dx + dy * dy;
  Phi(j) = rbf(ds_sq, cuttoff_dis_sq);
  OMEGA_H_CHECK_PRINTF(!std::isnan(Phi(j)),
                       "ERROR: Phi(j) in PhiVector is NaN for j = %d "
                       "ds_sq=%.16f, cuttoff_dis_sq=%.16f",
                       j, ds_sq, cuttoff_dis_sq);
}

// matrix matrix multiplication
KOKKOS_INLINE_FUNCTION
void MatMatMul(member_type team, ScratchMatView& moment_matrix,
               const ScratchMatView& pt_phi, const ScratchMatView& vandermonde)
{
  int M = pt_phi.extent(0);
  int N = vandermonde.extent(1);
  int K = pt_phi.extent(1);

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, M), [=](const int i) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int j) {
      double sum = 0.0;
      Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(team, K),
        [=](const int k, double& lsum) {
          OMEGA_H_CHECK_PRINTF(!std::isnan(pt_phi(i, k)),
                               "ERROR: pt_phi is NaN for i = %d\n", i);
          OMEGA_H_CHECK_PRINTF(!std::isnan(vandermonde(k, j)),
                               "ERROR: vandermonde is NaN for k = %d\n", k);
          lsum += pt_phi(i, k) * vandermonde(k, j);
        },
        sum);
      OMEGA_H_CHECK_PRINTF(!std::isnan(sum),
                           "ERROR: sum is NaN for i = %d, j = %d\n", i, j);
      moment_matrix(i, j) = sum;
    });
  });
}

// Matrix vector multiplication
KOKKOS_INLINE_FUNCTION
void MatVecMul(member_type team, const ScratchVecView& vector,
               const ScratchMatView& matrix, ScratchVecView& result)
{
  int M = matrix.extent(0);
  int N = matrix.extent(1);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int i) {
    double sum = 0;
    Kokkos::parallel_reduce(
      Kokkos::ThreadVectorRange(team, M),
      [=](const int j, double& lsum) { lsum += vector(j) * matrix(j, i); },
      sum);
    result(i) = sum;
    OMEGA_H_CHECK_PRINTF(!std::isnan(result(i)),
                         "ERROR: sum is NaN for i = %d\n", i);
  });
  // team.team_barrier();
}

// dot product
KOKKOS_INLINE_FUNCTION
void dot_product(member_type team, const ScratchVecView& result_sub,
                 const ScratchVecView& SupportValues_sub, double& target_value)
{
  int N = result_sub.extent(0);
  for (int j = 0; j < N; ++j) {
    target_value += result_sub(j) * SupportValues_sub[j];
  }
  OMEGA_H_CHECK_PRINTF(!std::isnan(target_value),
                       "ERROR: NaN found in dot_product: N: %d\n", N);
}

// moment matrix
KOKKOS_INLINE_FUNCTION
void PtphiPMatrix(ScratchMatView& moment_matrix, member_type team,
                  const ScratchMatView& pt_phi,
                  const ScratchMatView& vandermonde)
{
  int M = pt_phi.extent(0);
  int N = vandermonde.extent(1);
  int K = pt_phi.extent(1);

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, M), [=](const int i) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int j) {
      double sum = 0.0;
      Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(team, K),
        [=](const int k, double& lsum) {
          lsum += pt_phi(i, k) * vandermonde(k, j);
        },
        sum);
      moment_matrix(i, j) = sum;
    });
  });
}

// inverse matrix
KOKKOS_INLINE_FUNCTION
void inverse_matrix(member_type team, const ScratchMatView& matrix,
                    ScratchMatView& lower, ScratchMatView& forward_matrix,
                    ScratchMatView& solution)
{
  int N = matrix.extent(0);

  for (int j = 0; j < N; ++j) {
    Kokkos::single(Kokkos::PerTeam(team), [=]() {
      double sum = 0;
      for (int k = 0; k < j; ++k) {
        sum += lower(j, k) * lower(j, k);
      }
      OMEGA_H_CHECK_PRINTF(!std::isnan(matrix(j, j)),
                           "ERROR: matrix(j,j) is NaN: j = %d\n", j);
      OMEGA_H_CHECK_PRINTF(
        matrix(j, j) - sum >= -1e-10, // TODO: how to check this reliably?
        "ERROR: (matrix(j,j) - sum) is negative: mat(jj)=%.16f, sum = %.16f \n",
        matrix(j, j), sum);
      lower(j, j) = sqrt(matrix(j, j) - sum);
      OMEGA_H_CHECK_PRINTF(!std::isnan(lower(j, j)),
                           "Lower in inverse_matrix is NaN (j,j) = (%d,%d)\n",
                           j, j);
    });

    team.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, j + 1, N), [=](int i) {
      double inner_sum = 0;
      for (int k = 0; k < j; ++k) {
        inner_sum += lower(i, k) * lower(j, k);
      }
      lower(i, j) = (matrix(i, j) - inner_sum) / lower(j, j);
      lower(j, i) = lower(i, j);
      OMEGA_H_CHECK_PRINTF(!std::isnan(lower(i, j)),
                           "Lower in inverse_matrix is NaN (i,j) = (%d,%d)\n",
                           i, j);
      OMEGA_H_CHECK_PRINTF(!std::isnan(lower(j, i)),
                           "Lower in inverse_matrix is NaN (j,i) = (%d,%d)\n",
                           j, i);
    });

    team.team_barrier();
  }

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int i) {
    forward_matrix(i, i) = 1.0 / lower(i, i);
    for (int j = i + 1; j < N; ++j) {
      forward_matrix(j, i) = 0.0; // Initialize to zero
      for (int k = 0; k < j; ++k) {
        forward_matrix(j, i) -= lower(j, k) * forward_matrix(k, i);
      }
      forward_matrix(j, i) /= lower(j, j);
      OMEGA_H_CHECK_PRINTF(!std::isnan(forward_matrix(j, i)),
                           "Forward in inverse_matrix is NaN (j,i) = (%d,%d)\n",
                           j, i);
    }
  });

  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int i) {
    solution(N - 1, i) = forward_matrix(N - 1, i) / lower(N - 1, N - 1);
    for (int j = N - 2; j >= 0; --j) {
      solution(j, i) = forward_matrix(j, i);
      for (int k = j + 1; k < N; ++k) {
        solution(j, i) -= lower(j, k) * solution(k, i);
      }
      solution(j, i) /= lower(j, j);
      OMEGA_H_CHECK_PRINTF(
        !std::isnan(solution(j, i)),
        "Solution in inverse_matrix is NaN (j,i) = (%d,%d)\n", j, i);
    }
  });
}

#endif
