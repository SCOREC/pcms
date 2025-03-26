#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/mls_interpolation.hpp>
#include <pcms/interpolator/mls_interpolation_impl.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

using namespace pcms;

TEST_CASE("basis test")
{
  SECTION("check basis slice lengths, size of the basis vector and basis "
          "vector for degree = 3 and dim = 3")
  {

    const int dim = 3;
    const int degree = 3;
    Kokkos::View<int**, Kokkos::HostSpace> array("array", degree, dim);
    Kokkos::deep_copy(array, 0);
    detail::calculate_basis_slice_lengths(array);
    auto size = detail::calculate_basis_vector_size(array);

    int expected[degree][dim] = {{1, 1, 1}, {1, 2, 3}, {1, 3, 6}};

    for (int i = 0; i < degree; ++i) {
      for (int j = 0; j < dim; ++j) {
        REQUIRE(array(i, j) == expected[i][j]);
      }
    }

    REQUIRE(size == 20);

    IntDeviceMatView d_array("array in device", degree, dim);
    auto array_hd = Kokkos::create_mirror_view(d_array);
    Kokkos::deep_copy(array_hd, array);
    Kokkos::deep_copy(d_array, array_hd);

    int nvertices_target = 2;
    Omega_h::Matrix<3, 2> coords{{2.0, 3.0, 4.0}, {0.4, 0.5, 0.2}};
    Kokkos::View<double**> results("stores results", nvertices_target, size);

    auto host_results = Kokkos::create_mirror_view(results);

    team_policy tp(nvertices_target, Kokkos::AUTO);
    Kokkos::parallel_for(
      "inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
      KOKKOS_LAMBDA(const member_type& team) {
        int i = team.league_rank();
        ScratchVecView basis_vector(team.team_scratch(0), size);
        pcms::detail::fill(0, team, basis_vector);

        double target_point[MAX_DIM] = {};

        for (int j = 0; j < dim; ++j) {
          target_point[j] = coords[i][j];
        }

        detail::eval_basis_vector(d_array, target_point, basis_vector);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size),
                             [=](int j) { results(i, j) = basis_vector(j); });
      });

    Kokkos::deep_copy(host_results, results);
    Kokkos::View<double*, Kokkos::HostSpace> expected_basis("expected vector",
                                                            size);
    std::vector<std::string> names = {"1",   "x",    "y",    "z",    "x^2",
                                      "xy",  "y^2",  "xz",   "yz",   "z^2",
                                      "x^3", "x^2y", "xy^2", "Y^3",  "x^2z",
                                      "xyz", "y^2z", "xz^2", "yz^2", "z^3"};

    for (int i = 0; i < nvertices_target; ++i) {
      double x = coords[i][0];
      double y = coords[i][1];
      double z = coords[i][2];
      expected_basis(0) = 1;
      expected_basis(1) = x;          // (x)
      expected_basis(2) = y;          // (y)
      expected_basis(3) = z;          // (z)
      expected_basis(4) = x * x;      // (x^2)
      expected_basis(5) = x * y;      // (xy)
      expected_basis(6) = y * y;      // (y^2)
      expected_basis(7) = x * z;      // (xz)
      expected_basis(8) = y * z;      // (yz)
      expected_basis(9) = z * z;      // (z^2)
      expected_basis(10) = x * x * x; // (x^3)
      expected_basis(11) = x * x * y; // (x^2y)
      expected_basis(12) = x * y * y; // (xy^2)
      expected_basis(13) = y * y * y; // (y^3)
      expected_basis(14) = x * x * z; // (x^2z)
      expected_basis(15) = x * y * z; // (xyz)
      expected_basis(16) = y * y * z; // (y^2z)
      expected_basis(17) = x * z * z; // (xz^2)
      expected_basis(18) = y * z * z; // (yz^2)
      expected_basis(19) = z * z * z; // (z^3)
      for (int j = 0; j < size; ++j) {
        std::cout << names[j] << " " << expected_basis(j) << " "
                  << host_results(i, j) << "\n";
        REQUIRE(expected_basis(j) == host_results(i, j));
      }
    }
  }

  SECTION("check basis slice lengths, size of the basis vector and basis "
          "vector for degree = 1 and dim = 3")
  {

    const int dim = 3;
    const int degree = 1;
    Kokkos::View<int**, Kokkos::HostSpace> array("array", degree, dim);
    Kokkos::deep_copy(array, 0);
    detail::calculate_basis_slice_lengths(array);

    auto size = detail::calculate_basis_vector_size(array);

    int expected[degree][dim] = {
      {1, 1, 1},
    };

    for (int i = 0; i < degree; ++i) {
      for (int j = 0; j < dim; ++j) {
        REQUIRE(array(i, j) == expected[i][j]);
      }
    }

    REQUIRE(size == 4);

    IntDeviceMatView d_array("array in device", degree, dim);
    auto array_hd = Kokkos::create_mirror_view(d_array);
    Kokkos::deep_copy(array_hd, array);
    Kokkos::deep_copy(d_array, array_hd);

    int nvertices_target = 2;
    Omega_h::Matrix<3, 2> coords{{2.0, 3.0, 4.0}, {0.4, 0.5, 0.2}};
    Kokkos::View<double**> results("stores results", nvertices_target, size);

    auto host_results = Kokkos::create_mirror_view(results);

    team_policy tp(nvertices_target, Kokkos::AUTO);
    Kokkos::parallel_for(
      "inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
      KOKKOS_LAMBDA(const member_type& team) {
        int i = team.league_rank();
        ScratchVecView basis_vector(team.team_scratch(0), size);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size),
                             [=](int j) { basis_vector(j) = 0; });

        double target_point[MAX_DIM] = {};

        for (int j = 0; j < dim; ++j) {
          target_point[j] = coords[i][j];
        }

        detail::eval_basis_vector(d_array, target_point, basis_vector);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size),
                             [=](int j) { results(i, j) = basis_vector(j); });
      });

    Kokkos::deep_copy(host_results, results);
    Kokkos::View<double*, Kokkos::HostSpace> expected_basis("expected vector",
                                                            size);
    std::vector<std::string> names = {"1", "x", "y", "z"};

    for (int i = 0; i < nvertices_target; ++i) {
      double x = coords[i][0];
      double y = coords[i][1];
      double z = coords[i][2];
      expected_basis(0) = 1;
      expected_basis(1) = x; // (x)
      expected_basis(2) = y; // (y)
      expected_basis(3) = z; // (z)
      for (int j = 0; j < size; ++j) {
        std::cout << names[j] << " " << expected_basis(j) << " "
                  << host_results(i, j) << "\n";
        REQUIRE(expected_basis(j) == host_results(i, j));
      }
    }
  }

  SECTION("check vandermonde matrix")
  {

    const int dim = 1;
    const int degree = 3;

    Kokkos::View<int**, Kokkos::HostSpace> host_slice_length("array", degree,
                                                             dim);
    Kokkos::deep_copy(host_slice_length, 0);
    detail::calculate_basis_slice_lengths(host_slice_length);

    IntDeviceMatView slice_length("slice array in device", degree, dim);
    auto slice_length_hd = Kokkos::create_mirror_view(slice_length);
    Kokkos::deep_copy(slice_length_hd, host_slice_length);
    Kokkos::deep_copy(slice_length, slice_length_hd);
    auto size = detail::calculate_basis_vector_size(host_slice_length);
    int nvertices_target = 2;
    int nsupports = 5;
    Kokkos::View<double**> points("vandermonde matrix", nvertices_target,
                                  nsupports);

    Kokkos::View<double***> vandermonde_matrix_combined(
      "stores vandermonde matrix", nvertices_target, nsupports, size);

    auto host_points = Kokkos::create_mirror_view(points);
    auto host_vandermonde_matrix_combined =
      Kokkos::create_mirror_view(vandermonde_matrix_combined);

    host_points(0, 1) = 0.5;
    host_points(0, 2) = 0.1;
    host_points(0, 3) = 0.3;
    host_points(0, 4) = 0.2;
    host_points(0, 5) = 0.6;

    host_points(1, 1) = 1.1;
    host_points(1, 2) = 2.6;
    host_points(1, 3) = 0.8;
    host_points(1, 4) = 0.4;
    host_points(1, 5) = 1.7;

    Kokkos::deep_copy(points, host_points);

    team_policy tp(nvertices_target, Kokkos::AUTO);
    Kokkos::parallel_for(
      "inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
      KOKKOS_LAMBDA(const member_type& team) {
        int i = team.league_rank();
        ScratchMatView vandermonde_matrix(team.team_scratch(0), nsupports,
                                          size);

        ScratchMatView local_source_points(team.team_scratch(0), nsupports, 1);
        for (int j = 0; j < nsupports; ++j) {
          local_source_points(j, 0) = points(i, j);
        }

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports),
                             [=](int j) {
                               for (int k = 0; k < size; ++k) {
                                 vandermonde_matrix(j, k) = 0.0;
                               }
                             });

        team.team_barrier();

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
            detail::create_vandermonde_matrix(local_source_points, j,
                                              slice_length, vandermonde_matrix);
          });

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
            for (int k = 0; k < size; ++k) {
              vandermonde_matrix_combined(i, j, k) = vandermonde_matrix(j, k);
            }
          });
      });

    Kokkos::deep_copy(host_vandermonde_matrix_combined,
                      vandermonde_matrix_combined);
    Kokkos::View<double***, Kokkos::HostSpace> expected_vandermonde_matrix(
      "expected vector", nvertices_target, nsupports, size);

    for (int i = 0; i < nvertices_target; ++i) {
      for (int j = 0; j < nsupports; ++j) {
        double x = host_points(i, j);
        expected_vandermonde_matrix(i, j, 0) = 1;
        expected_vandermonde_matrix(i, j, 1) = x;         // (x)
        expected_vandermonde_matrix(i, j, 2) = x * x;     // (x^2)
        expected_vandermonde_matrix(i, j, 3) = x * x * x; // (x^3)
      }
    }

    for (int i = 0; i < nvertices_target; ++i) {
      for (int j = 0; j < nsupports; ++j) {
        for (int k = 0; k < size; ++k) {
          REQUIRE(expected_vandermonde_matrix(i, j, k) ==
                  host_vandermonde_matrix_combined(i, j, k));
        }
      }
    }
  }
}
