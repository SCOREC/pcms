#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>

#include <pcms/transfer/mass_matrix_integrator.hpp>

TEST_CASE("Local mass matrix on reference triangle for linear elements",
          "[mass_matrix][local]")
{
  Omega_h::Library lib;

  // Reference triangle:
  // v0 = (0,0), v1 = (1,0), v2 = (0,1)
  Omega_h::Reals coords({
    0.0, 0.0, // v0
    1.0, 0.0, // v1
    0.0, 1.0  // v2
  });

  Omega_h::LOs ev2v({0, 1, 2});

  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, 2, ev2v, coords);

  REQUIRE(mesh.dim() == 2);
  REQUIRE(mesh.nelems() == 1);
  REQUIRE(mesh.nverts() == 3);

  MeshField::OmegahMeshField<DefaultExecutionSpace, 2,
                             MeshField::KokkosController>
    omf(mesh);

  constexpr int ShapeOrder = 1;
  auto coordField = omf.getCoordField();
  const auto [shp, map] =
    MeshField::Omegah::getTriangleElement<ShapeOrder>(mesh);
  MeshField::FieldElement coordFe(mesh.nelems(), coordField, shp, map);

  auto elemMass = buildMassMatrix(mesh, coordFe);

  auto elemMass_h = Kokkos::create_mirror_view(elemMass);
  Kokkos::deep_copy(elemMass_h, elemMass);

  // Assuming one element and 9 flattened entries per element:
  // [M00, M01, M02, M10, M11, M12, M20, M21, M22]
  REQUIRE(elemMass_h.size() == 9);

  // Reference-triangle P1 mass matrix:
  // (1/24) * [ [2,1,1], [1,2,1], [1,1,2] ]
  const double expected[9] = {
    1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0,
    1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0,
    1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0
  };

  const double tol = 1e-12;

  SECTION("Check local mass matrix entries")
  {
    for (int k = 0; k < 9; ++k) {
      CAPTURE(k, expected[k], elemMass_h(0, k));
      CHECK_THAT(elemMass_h(0, k),
                 Catch::Matchers::WithinAbs(expected[k], tol));
    }
  }

  SECTION("Check symmetry")
  {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const int a = 3 * i + j;
        const int b = 3 * j + i;
        CHECK_THAT(elemMass_h(0, a),
                   Catch::Matchers::WithinAbs(elemMass_h(0, b), tol));
      }
    }
  }

  SECTION("Check row sums")
  {
    // Each row sum = integral of basis function over ref triangle = area/3 = 1/6
    const double expected_row_sum = 1.0 / 6.0;

    for (int i = 0; i < 3; ++i) {
      double row_sum = 0.0;
      for (int j = 0; j < 3; ++j) {
        row_sum += elemMass_h(0, 3 * i + j);
      }
      CAPTURE(i, row_sum);
      CHECK_THAT(row_sum,
                 Catch::Matchers::WithinAbs(expected_row_sum, tol));
    }
  }
}
