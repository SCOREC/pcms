#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/mls_interpolation_impl.hpp>
#include <pcms/interpolator/pcms_interpolator_aliases.hpp>
#include <pcms/interpolator/pcms_interpolator_view_utils.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "KokkosBatched_SVD_Decl.hpp"
#include "KokkosBatched_SVD_Serial_Impl.hpp"
#include <KokkosBlas2_serial_gemv_impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <vector>
#include <iostream>

using namespace std;
using namespace Omega_h;
using namespace pcms;

TEST_CASE("test_serial_svd")
{

  constexpr int row = 3;
  constexpr int column = 3;

  double tolerance = 1e-4;
  double A_data[row][column] = {
    {4.0, 2.0, 3.0}, {3.0, 5.0, 1.0}, {2.0, 3.0, 6.0}};

  double rhs_data[row] = {1.28571, 1.42857, 3.1428};

  double expected_solution[column] = {-0.142857, 0.285714,
                                      0.428571}; // Approximate values

  SECTION("test_svd_factorization")
  {
    Kokkos::View<double**> result("result", row, column);
    team_policy tp(1, Kokkos::AUTO);
    Kokkos::parallel_for(
      "Solve SVD", tp.set_scratch_size(1, Kokkos::PerTeam(500)),
      KOKKOS_LAMBDA(const member_type& team) {
        ScratchMatView A(team.team_scratch(1), row, column);
        ScratchMatView U(team.team_scratch(1), row, row);
        ScratchMatView Vt(team.team_scratch(1), column, column);
        ScratchVecView sigma(team.team_scratch(1), column);
        ScratchVecView work(team.team_scratch(1), row);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            A(i, j) = A_data[i][j];
          }
        });

        if (team.league_rank() == 0) {
          KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), A, U,
                                           sigma, Vt, work);
        }
        team.team_barrier();
        ScratchMatView sigma_mat(team.team_scratch(1), row, column);
        ScratchMatView Usigma(team.team_scratch(1), row, column);
        detail::fill(0.0, team, Usigma);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            sigma_mat(i, j) = (i == j) ? sigma(i) : 0.0;
          }
        });
        team.team_barrier();

        ScratchMatView reconstructedA(team.team_scratch(1), row, column);
        if (team.league_rank() == 0) {
          KokkosBatched::SerialGemm<
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, U, sigma_mat,
                                                          0.0, Usigma);

          KokkosBatched::SerialGemm<
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, Usigma, Vt, 0.0,
                                                          reconstructedA);
        }

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            result(i, j) = reconstructedA(i, j);
          }
        });
      });
    auto host_result = Kokkos::create_mirror_view(result);
    Kokkos::deep_copy(host_result, result);

    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < column; ++j) {
        CHECK_THAT(host_result(i, j),
                   Catch::Matchers::WithinAbs(A_data[i][j], tolerance));
      }
    }
  } // end section

  SECTION("test_transpose_function")
  {

    Kokkos::View<double**> result("result", row, row);
    Kokkos::View<double**> transpose_expected("result", row, row);
    Kokkos::deep_copy(result, 0.0);
    Kokkos::deep_copy(transpose_expected, 0.0);
    team_policy tp(1, Kokkos::AUTO);
    Kokkos::parallel_for(
      "Solve SVD", tp.set_scratch_size(1, Kokkos::PerTeam(500)),
      KOKKOS_LAMBDA(const member_type& team) {
        ScratchMatView A(team.team_scratch(1), row, column);
        ScratchMatView U(team.team_scratch(1), row, row);
        ScratchMatView Vt(team.team_scratch(1), column, column);
        ScratchVecView sigma(team.team_scratch(1), column);
        ScratchVecView work(team.team_scratch(1), row);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            A(i, j) = A_data[i][j];
          }
        });

        if (team.league_rank() == 0) {
          KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), A, U,
                                           sigma, Vt, work);
        }
        team.team_barrier();
        // detail::fill(0.0, team, sigma_inv_expected);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < row; ++j) {
            transpose_expected(i, j) = U(j, i);
          }
        });
        auto transposedU = detail::find_transpose(team, U);
        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < row; ++j) {
            result(i, j) = transposedU(i, j);
          }
        });
      });
    auto host_result = Kokkos::create_mirror_view(result);
    auto host_expected = Kokkos::create_mirror_view(transpose_expected);
    Kokkos::deep_copy(host_result, result);
    Kokkos::deep_copy(host_expected, transpose_expected);

    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < row; ++j) {
        CHECK_THAT(host_result(i, j),
                   Catch::Matchers::WithinAbs(host_expected(i, j), tolerance));
      }
    }

  } // end section

  SECTION("test_row_scaling_function")
  {
    double diagonal_entries[row] = {2.0, 0.5, -1.0};
    double expected_result[row][column] = {
      {8.0, 4.0, 6.0}, {1.5, 2.5, 0.5}, {-2.0, -3.0, -6.0}};

    Kokkos::View<double**> result("result", row, column);
    Kokkos::deep_copy(result, 0.0);
    team_policy tp(1, Kokkos::AUTO);
    Kokkos::parallel_for(
      "Solve SVD", tp.set_scratch_size(1, Kokkos::PerTeam(500)),
      KOKKOS_LAMBDA(const member_type& team) {
        ScratchMatView A(team.team_scratch(1), row, column);
        ScratchVecView weight(team.team_scratch(1), row);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            A(i, j) = A_data[i][j];
          }
          weight(i) = diagonal_entries[i];
        });

        team.team_barrier();

        detail::eval_row_scaling(team, weight, A);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            result(i, j) = A(i, j);
          }
        });
      });
    auto host_result = Kokkos::create_mirror_view(result);
    Kokkos::deep_copy(host_result, result);

    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < row; ++j) {
        CHECK_THAT(host_result(i, j), Catch::Matchers::WithinAbs(
                                        expected_result[i][j], tolerance));
      }
    }

  } // end section

  SECTION("test_scale_adjust_function")
  {
    constexpr int rowA = 6; // Matrix to scale (6x6)
    constexpr int colA = 6;
    constexpr int rowB = 3; // Adjusted matrix (3x6)
    constexpr int colB = 6;

    double A_data[rowA][colA] = {
      {4.0, 2.0, 3.0, 1.0, 0.0, 5.0}, {3.0, 5.0, 1.0, 4.0, 2.0, 6.0},
      {2.0, 3.0, 6.0, 5.0, 1.0, 7.0}, {7.0, 8.0, 9.0, 3.0, 2.0, 1.0},
      {1.0, 3.0, 2.0, 6.0, 4.0, 5.0}, {5.0, 7.0, 3.0, 8.0, 9.0, 6.0}};

    double diagonal_entries_data[rowB] = {2.0, 0.5, -1.0};

    double expected_result[rowB][colB] = {{8.0, 4.0, 6.0, 2.0, 0.0, 10.0},
                                          {1.5, 2.5, 0.5, 2.0, 1.0, 3.0},
                                          {-2.0, -3.0, -6.0, -5.0, -1.0, -7.0}};

    Kokkos::View<double**> result("result", rowB, colB);
    Kokkos::deep_copy(result, 0.0);
    team_policy tp(1, Kokkos::AUTO);
    Kokkos::parallel_for(
      tp.set_scratch_size(1, Kokkos::PerTeam(500)),
      KOKKOS_LAMBDA(const member_type& team) {
        ScratchMatView A(team.team_scratch(1), rowA, colA);
        ScratchMatView result_matrix(team.team_scratch(1), rowB, colB);
        ScratchVecView weight(team.team_scratch(1), rowB);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowA), [=](int i) {
          for (int j = 0; j < colA; ++j) {
            A(i, j) = A_data[i][j];
          }
        });

        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowB), [=](int i) {
          weight(i) = diagonal_entries_data[i];
        });

        team.team_barrier();

        detail::scale_and_adjust(team, weight, A, result_matrix);
        team.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowB), [=](int i) {
          for (int j = 0; j < colB; ++j) {
            result(i, j) = result_matrix(i, j);
          }
        });
      });
    auto host_result = Kokkos::create_mirror_view(result);
    Kokkos::deep_copy(host_result, result);

    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < row; ++j) {
        CHECK_THAT(host_result(i, j), Catch::Matchers::WithinAbs(
                                        expected_result[i][j], tolerance));
      }
    }

  } // end section

  SECTION("test_svd_solver")
  {
    Kokkos::View<double*> result("result", column);
    team_policy tp(1, Kokkos::AUTO);
    Kokkos::parallel_for(
      "Solve SVD", tp.set_scratch_size(1, Kokkos::PerTeam(800)),
      KOKKOS_LAMBDA(const member_type& team) {
        ScratchMatView A(team.team_scratch(1), row, column);
        ScratchVecView rhs(team.team_scratch(1), row);
        ScratchVecView x(team.team_scratch(1), column);
        detail::fill(0.0, team, x);

        ScratchVecView weight(team.team_scratch(1), row);
        detail::fill(1.0, team, weight);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            A(i, j) = A_data[i][j];
          }
          rhs(i) = rhs_data[i];
        });

        detail::solve_matrix_svd(team, weight, rhs, A, x);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, column),
                             [=](int i) { result(i) = x(i); });
      });
    auto host_result = Kokkos::create_mirror_view(result);
    Kokkos::deep_copy(host_result, result);

    for (int i = 0; i < column; ++i) {
      CHECK_THAT(host_result(i),
                 Catch::Matchers::WithinAbs(expected_solution[i], tolerance));
    }
  } // end section
}
