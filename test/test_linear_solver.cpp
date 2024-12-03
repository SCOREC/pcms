#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/MLSCoefficients.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_fail.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("solver test")
{
  // This test computes A^TQA and A^TQ b for the normal equation A^TQA x = A^TQ
  // b; where A(n,m) (n >= m) is a rectangular matrix, Q(m,m) is a diagonal
  // matrix with diagonal elements stored as vector and b(n,1) is a vector
  SECTION("check convert normal equation ")
  {

    int nvertices_target = 2;

    int nsupports = 3;

    int size = 2;

    Kokkos::View<double**> Q("phi vector", nvertices_target, nsupports);

    Kokkos::View<double**> b("rhs vector", nvertices_target, nsupports);

    Kokkos::View<double***> A("vandermonde matrix", nvertices_target, nsupports,
                              size);

    Kokkos::View<double***> TransAQA("moment matrix", nvertices_target, size,
                                     size);

    Kokkos::View<double***> TransAQ("moment matrix", nvertices_target, size,
                                    nsupports);

    Kokkos::View<double**> TransAQb("transformed rhs", nvertices_target, size);

    Kokkos::View<double**> solution("transformed rhs", nvertices_target, size);

    auto host_Q = Kokkos::create_mirror_view(Q);
    auto host_A = Kokkos::create_mirror_view(A);
    auto host_b = Kokkos::create_mirror_view(b);

    host_Q(0, 0) = 1.0;
    host_Q(0, 1) = 1.0;
    host_Q(0, 2) = 1.0;

    host_Q(1, 0) = 3.0;
    host_Q(1, 1) = 2.0;
    host_Q(1, 2) = 1.0;

    Kokkos::deep_copy(Q, host_Q);

    host_A(0, 0, 0) = 1.0;
    host_A(0, 0, 1) = 2.0;
    host_A(0, 1, 0) = 3.0;
    host_A(0, 1, 1) = 4.0;
    host_A(0, 2, 0) = 5.0;
    host_A(0, 2, 1) = 6.0;

    host_A(1, 0, 0) = 2.0;
    host_A(1, 0, 1) = 3.0;
    host_A(1, 1, 0) = 1.0;
    host_A(1, 1, 1) = -1.0;
    host_A(1, 2, 0) = 4.0;
    host_A(1, 2, 1) = 1.0;

    Kokkos::deep_copy(A, host_A);

    host_b(0, 0) = 7.0;
    host_b(0, 1) = 8.0;
    host_b(0, 2) = 9.0;

    host_b(1, 0) = 5.0;
    host_b(1, 1) = 6.0;
    host_b(1, 2) = 7.0;

    Kokkos::deep_copy(b, host_b);

    team_policy tp(nvertices_target, Kokkos::AUTO);
    Kokkos::parallel_for(
      "inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
      KOKKOS_LAMBDA(const member_type& team) {
        int i = team.league_rank();
        ScratchMatView vandermonde_matrix(team.team_scratch(0), nsupports,
                                          size);
        ScratchVecView phi(team.team_scratch(0), nsupports);
        ScratchVecView support_values(team.team_scratch(0), nsupports);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports),
                             [=](int j) {
                               phi(j) = Q(i, j);
                               support_values(j) = b(i, j);
                               for (int k = 0; k < size; ++k) {
                                 vandermonde_matrix(j, k) = A(i, j, k);
                               }
                             });

        team.team_barrier();

        auto result =
          ConvertNormalEq(vandermonde_matrix, phi, support_values, team);

        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size), [=](int j) {
          TransAQb(i, j) = result.transformed_rhs(j);
          for (int k = 0; k < size; ++k) {
            TransAQA(i, j, k) = result.square_matrix(j, k);
          }
        });

        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports),
                             [=](int j) {
                               for (int k = 0; k < size; ++k) {
                                 TransAQ(i, k, j) = result.scaled_matrix(k, j);
                               }
                             });

        team.team_barrier();

        SolveMatrix(result.square_matrix, result.transformed_rhs, team);

        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size), [=](int j) {
          solution(i, j) = result.transformed_rhs(j);
        });

        team.team_barrier();
      });

    auto host_result_lhs = Kokkos::create_mirror_view(TransAQA);
    auto host_result_rhs = Kokkos::create_mirror_view(TransAQb);
    auto host_result_scaled = Kokkos::create_mirror_view(TransAQ);
    auto host_result_solution = Kokkos::create_mirror_view(solution);

    Kokkos::deep_copy(host_result_lhs, TransAQA);
    Kokkos::deep_copy(host_result_rhs, TransAQb);
    Kokkos::deep_copy(host_result_scaled, TransAQ);
    Kokkos::deep_copy(host_result_solution, solution);

    Kokkos::View<double***, Kokkos::HostSpace> expected_lhs_matrix(
      "expected lhs matrix", nvertices_target, size, size);

    Kokkos::View<double**, Kokkos::HostSpace> expected_rhs_vector(
      "expected rhs vector", nvertices_target, size);

    Kokkos::View<double***, Kokkos::HostSpace> expected_scaled_matrix(
      "expected scaled matrix", nvertices_target, size, nsupports);

    Kokkos::View<double**, Kokkos::HostSpace> expected_solution(
      "expected solution", nvertices_target, size);

    expected_lhs_matrix(0, 0, 0) = 35.0;
    expected_lhs_matrix(0, 0, 1) = 44.0;
    expected_lhs_matrix(0, 1, 0) = 44.0;
    expected_lhs_matrix(0, 1, 1) = 56.0;

    expected_rhs_vector(0, 0) = 76.0;
    expected_rhs_vector(0, 1) = 100.0;

    expected_scaled_matrix(0, 0, 0) = 1.0;
    expected_scaled_matrix(0, 0, 1) = 3.0;
    expected_scaled_matrix(0, 0, 2) = 5.0;
    expected_scaled_matrix(0, 1, 0) = 2.0;
    expected_scaled_matrix(0, 1, 1) = 4.0;
    expected_scaled_matrix(0, 1, 2) = 6.0;

    expected_solution(0, 0) = -6.0;
    expected_solution(0, 1) = 6.5;

    expected_lhs_matrix(1, 0, 0) = 30.0;
    expected_lhs_matrix(1, 0, 1) = 20.0;
    expected_lhs_matrix(1, 1, 0) = 20.0;
    expected_lhs_matrix(1, 1, 1) = 30.0;

    expected_rhs_vector(1, 0) = 70.0;
    expected_rhs_vector(1, 1) = 40.0;

    expected_scaled_matrix(1, 0, 0) = 6.0;
    expected_scaled_matrix(1, 0, 1) = 2.0;
    expected_scaled_matrix(1, 0, 2) = 4.0;
    expected_scaled_matrix(1, 1, 0) = 9.0;
    expected_scaled_matrix(1, 1, 1) = -2.0;
    expected_scaled_matrix(1, 1, 2) = 1.0;

    expected_solution(1, 0) = 2.6;
    expected_solution(1, 1) = -0.4;

    for (int i = 0; i < nvertices_target; ++i) {
      for (int j = 0; j < size; ++j) {
        for (int l = 0; l < nsupports; ++l) {
          printf("scaled_matrix (A^T Q) (%f, %f) at (%d,%d,%d)\n",
                 expected_scaled_matrix(i, j, l), host_result_scaled(i, j, l),
                 i, j, l);
          REQUIRE_THAT(
            expected_scaled_matrix(i, j, l),
            Catch::Matchers::WithinAbs(host_result_scaled(i, j, l), 1E-10));
        }
      }
      for (int j = 0; j < size; ++j) {
        for (int k = 0; k < size; ++k) {
          printf("A^T Q A (%f, %f) at (%d,%d,%d)\n",
                 expected_lhs_matrix(i, j, k), host_result_lhs(i, j, k), i, j,
                 k);
          REQUIRE_THAT(
            expected_lhs_matrix(i, j, k),
            Catch::Matchers::WithinAbs(host_result_lhs(i, j, k), 1E-10));
        }
      }
      for (int j = 0; j < size; ++j) {
        printf("A^T Q b: (%f,%f) at (%d,%d)\n", expected_rhs_vector(i, j),
               host_result_rhs(i, j), i, j);
        REQUIRE_THAT(expected_rhs_vector(i, j),
                     Catch::Matchers::WithinAbs(host_result_rhs(i, j), 1E-10));
      }

      for (int j = 0; j < size; ++j) {
        printf("A^T Q b: (%f,%f) at (%d,%d)\n", expected_solution(i, j),
               host_result_solution(i, j), i, j);
        REQUIRE_THAT(
          expected_solution(i, j),
          Catch::Matchers::WithinAbs(host_result_solution(i, j), 1E-10));
      }
    }
  }
}
