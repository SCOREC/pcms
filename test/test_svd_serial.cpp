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

  Kokkos::View<double**> A_data("Device A data", row, column);
  auto host_A_data = Kokkos::create_mirror_view(A_data);

  host_A_data(0, 0) = 4.0;
  host_A_data(0, 1) = 2.0;
  host_A_data(0, 2) = 3.0;
  host_A_data(1, 0) = 3.0;
  host_A_data(1, 1) = 5.0;
  host_A_data(1, 2) = 1.0;
  host_A_data(2, 0) = 2.0;
  host_A_data(2, 1) = 3.0;
  host_A_data(2, 2) = 6.0;

  Kokkos::deep_copy(A_data, host_A_data);

  Kokkos::View<double*> rhs_data("Device rhs data", row, column);
  auto host_rhs_data = Kokkos::create_mirror_view(rhs_data);

  host_rhs_data(0) = 1.28571;
  host_rhs_data(1) = 1.42857;
  host_rhs_data(2) = 3.1428;

  Kokkos::deep_copy(rhs_data, host_rhs_data);

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
        ScratchMatView reconstructedA(team.team_scratch(1), row, column);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            A(i, j) = A_data(i, j);
          }
        });

        ScratchMatView sigma_mat(team.team_scratch(1), row, column);
        ScratchMatView Usigma(team.team_scratch(1), row, column);
        detail::fill(0.0, team, Usigma);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            sigma_mat(i, j) = (i == j) ? sigma(i) : 0.0;
          }
        });

        if (team.team_rank() == 0) {
          printf("I am above svd solver function");
          printf("A before SVD:\n");
          for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
              printf("%f ", A(i, j));
            }
            printf("\n");
          }

          KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), A, U,
                                           sigma, Vt, work, 1e-6);

          printf("=== Sigma Values ===\n");
          for (int i = 0; i < column; ++i) {
            printf("sigma[%d] = %f\n", i, sigma(i));
          }

          printf("\n=== First few entries of U ===\n");
          for (int i = 0; i < row; ++i) {
            for (int j = 0; j < row; ++j) {
              printf("%f ", U(i, j));
            }
            printf("\n");
          }

          printf("\n=== First few entries of Vt ===\n");
          for (int i = 0; i < column; ++i) {
            for (int j = 0; j < column; ++j) {
              printf("%f ", Vt(i, j));
            }
            printf("\n");
          }

          printf("=============================\n");

          printf("I am below svd solver function");
          KokkosBatched::TeamGemm<
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, U, sigma_mat,
                                                          0.0, Usigma);

          printf("I am below serialgemm");
          KokkosBatched::TeamGemm<
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Algo::Gemm::Unblocked>::invoke(1.0, Usigma, Vt, 0.0,
                                                          reconstructedA);
          printf("I am below another serialgemm");
        }
        team.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
          for (int j = 0; j < column; ++j) {
            result(i, j) = reconstructedA(i, j);
            printf("Copied reconstructedA(%d,%d) = %f\n", i, j,
                   reconstructedA(i, j));
          }
        });
      });

    Kokkos::fence();
    auto host_result = Kokkos::create_mirror_view(result);
    Kokkos::deep_copy(host_result, result);

    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < column; ++j) {
        CHECK_THAT(host_result(i, j),
                   Catch::Matchers::WithinAbs(host_A_data(i, j), tolerance));
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
            A(i, j) = A_data(i, j);
          }
        });

        if (team.team_rank() == 0) {
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
    Kokkos::View<double*> diagonal_entries("Device diagonal enetries", row);
    auto host_diagonal_entries = Kokkos::create_mirror_view(diagonal_entries);

    host_diagonal_entries(0) = 2.0;
    host_diagonal_entries(1) = 0.5;
    host_diagonal_entries(2) = -1.0;

    Kokkos::deep_copy(diagonal_entries, host_diagonal_entries);

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
            A(i, j) = A_data(i, j);
          }
          weight(i) = diagonal_entries(i);
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

    Kokkos::View<double**> A_data("Device A data", row, column);
    auto host_A_data = Kokkos::create_mirror_view(A_data);

    host_A_data(0, 0) = 4.0;
    host_A_data(0, 1) = 2.0;
    host_A_data(0, 2) = 3.0;
    host_A_data(0, 3) = 1.0;
    host_A_data(0, 4) = 0.0;
    host_A_data(0, 5) = 5.0;
    host_A_data(1, 0) = 3.0;
    host_A_data(1, 1) = 5.0;
    host_A_data(1, 2) = 1.0;
    host_A_data(1, 3) = 4.0;
    host_A_data(1, 4) = 2.0;
    host_A_data(1, 5) = 6.0;
    host_A_data(2, 0) = 2.0;
    host_A_data(2, 1) = 3.0;
    host_A_data(2, 2) = 6.0;
    host_A_data(2, 3) = 5.0;
    host_A_data(2, 4) = 1.0;
    host_A_data(2, 5) = 7.0;
    host_A_data(3, 0) = 7.0;
    host_A_data(3, 1) = 8.0;
    host_A_data(3, 2) = 9.0;
    host_A_data(3, 3) = 3.0;
    host_A_data(3, 4) = 2.0;
    host_A_data(3, 5) = 1.0;
    host_A_data(4, 0) = 1.0;
    host_A_data(4, 1) = 3.0;
    host_A_data(4, 2) = 2.0;
    host_A_data(4, 3) = 6.0;
    host_A_data(4, 4) = 4.0;
    host_A_data(4, 5) = 5.0;
    host_A_data(5, 0) = 5.0;
    host_A_data(5, 1) = 7.0;
    host_A_data(5, 2) = 3.0;
    host_A_data(5, 3) = 8.0;
    host_A_data(5, 4) = 9.0;
    host_A_data(5, 5) = 6.0;

    Kokkos::deep_copy(A_data, host_A_data);

    Kokkos::View<double*> diagonal_entries_data("Device rhs data");
    auto host_diagonal_entries_data =
      Kokkos::create_mirror_view(diagonal_entries_data);

    host_diagonal_entries_data(0) = 2.0;
    host_diagonal_entries_data(1) = 0.5;
    host_diagonal_entries_data(2) = -1.0;

    Kokkos::deep_copy(diagonal_entries_data, host_diagonal_entries_data);

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
            A(i, j) = A_data(i, j);
          }
        });

        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowB), [=](int i) {
          weight(i) = diagonal_entries_data(i);
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
            A(i, j) = A_data(i, j);
          }
          rhs(i) = rhs_data(i);
        });

        detail::solve_matrix_svd(team, weight, rhs, A, x, 0);
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
