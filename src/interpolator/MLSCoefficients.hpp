#ifndef MLS_COEFFICIENTS_HPP
#define MLS_COEFFICIENTS_HPP

#include <cmath>

#include "points.hpp"

#define PI_M 3.14159265358979323846

KOKKOS_INLINE_FUNCTION
double func(Coord& p) {
  auto x = (p.x - 0.5) * PI_M * 2;
  auto y = (p.y - 0.5) * PI_M * 2;
  double Z = sin(x) * sin(y);
  return Z;
}

// polynomial basis vector
KOKKOS_INLINE_FUNCTION
void BasisPoly(ScratchVecView basis_monomial, Coord& p1) {
  basis_monomial(0) = 1.0;
  basis_monomial(1) = p1.x;
  basis_monomial(2) = p1.y;
  basis_monomial(3) = p1.x * p1.x;
  basis_monomial(4) = p1.x * p1.y;
  basis_monomial(5) = p1.y * p1.y;
}

// radial basis function
KOKKOS_INLINE_FUNCTION
double rbf(double r_sq, double rho_sq) {
  double phi;
  double r = sqrt(r_sq);
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

  return phi;
}

// create vandermondeMatrix
KOKKOS_INLINE_FUNCTION
void VandermondeMatrix(ScratchMatView V, ScratchMatView local_source_points,
                       int j) {
  int N = local_source_points.extent(0);
  Coord source_point;
  source_point.x = local_source_points(j, 0);
  source_point.y = local_source_points(j, 1);
  ScratchVecView basis_monomial = Kokkos::subview(V, j, Kokkos::ALL());
  BasisPoly(basis_monomial, source_point);
}

// moment matrix
KOKKOS_INLINE_FUNCTION
void PTphiMatrix(ScratchMatView pt_phi, ScratchMatView V, ScratchVecView Phi,
                 int j) {
  int M = V.extent(0);
  int N = V.extent(1);

  ScratchVecView vandermonde_mat_row = Kokkos::subview(V, j, Kokkos::ALL());
  for (int k = 0; k < N; k++) {
    pt_phi(k, j) = vandermonde_mat_row(k) * Phi(j);
  }
}

// radial basis function vector
KOKKOS_INLINE_FUNCTION
void PhiVector(ScratchVecView Phi, Coord target_point,
               ScratchMatView local_source_points, int j,
               double cuttoff_dis_sq) {
  int N = local_source_points.extent(0);
  double dx = target_point.x - local_source_points(j, 0);
  double dy = target_point.y - local_source_points(j, 1);
  double ds_sq = dx * dx + dy * dy;
  Phi(j) = rbf(ds_sq, cuttoff_dis_sq);
}

// matrix matrix multiplication
KOKKOS_INLINE_FUNCTION
void MatMatMul(member_type team, ScratchMatView moment_matrix,
               ScratchMatView pt_phi, ScratchMatView vandermonde) {
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

// Matrix vector multiplication
KOKKOS_INLINE_FUNCTION
void MatVecMul(member_type team, ScratchVecView vector, ScratchMatView matrix,
               ScratchVecView result) {
  int M = matrix.extent(0);
  int N = matrix.extent(1);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int i) {
    double sum = 0;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(team, M),
        [=](const int j, double& lsum) { lsum += vector(j) * matrix(j, i); },
        sum);
    result(i) = sum;
  });
  // team.team_barrier();
}

// dot product
KOKKOS_INLINE_FUNCTION
void dot_product(member_type team, ScratchVecView result_sub,
                 ScratchVecView exactSupportValues_sub, double& target_value) {
  int N = result_sub.extent(0);
  for (int j = 0; j < N; ++j) {
    target_value += result_sub(j) * exactSupportValues_sub[j];
  }
}

// moment matrix
KOKKOS_INLINE_FUNCTION
void PtphiPMatrix(ScratchMatView moment_matrix, member_type team,
                  ScratchMatView pt_phi, ScratchMatView vandermonde) {
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
void inverse_matrix(member_type team, ScratchMatView matrix,
                    ScratchMatView lower, ScratchMatView forward_matrix,
                    ScratchMatView solution) {
  int N = matrix.extent(0);

  for (int j = 0; j < N; ++j) {
    Kokkos::single(Kokkos::PerTeam(team), [=]() {
      double sum = 0;
      for (int k = 0; k < j; ++k) {
        sum += lower(j, k) * lower(j, k);
      }
      lower(j, j) = sqrt(matrix(j, j) - sum);
    });

    team.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, j + 1, N), [=](int i) {
      double inner_sum = 0;
      for (int k = 0; k < j; ++k) {
        inner_sum += lower(i, k) * lower(j, k);
      }
      lower(i, j) = (matrix(i, j) - inner_sum) / lower(j, j);
      lower(j, i) = lower(i, j);
    });

    team.team_barrier();
  }

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [=](const int i) {
    forward_matrix(i, i) = 1.0 / lower(i, i);
    for (int j = i + 1; j < N; ++j) {
      forward_matrix(j, i) = 0.0;  // Initialize to zero
      for (int k = 0; k < j; ++k) {
        forward_matrix(j, i) -= lower(j, k) * forward_matrix(k, i);
      }
      forward_matrix(j, i) /= lower(j, j);
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
    }
  });
}

#endif
