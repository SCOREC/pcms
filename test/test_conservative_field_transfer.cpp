#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Omega_h_mesh.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_build.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <pcms/field_transfer/galerkin_projection_solver.hpp>
#include <pcms/field_transfer/mesh_intersection.hpp>

TEST_CASE("Galerkin projection mesh-to-mesh", "[projection]")
{

  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto rank = lib.world()->rank();

  // ----------------------------------------------
  // Build internal meshes (coarse target, fine source)
  // ----------------------------------------------
  int tgt_n = 4; // 4x4 = 32 triangles
  int src_n = 8; // 8x8 = 128 triangles (finer)

  Omega_h::Mesh tgt_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX,
                       /*Lx=*/1.0, /*Ly=*/1.0, /*Lz=*/0.0,
                       /*nx=*/tgt_n, /*ny=*/tgt_n, /*nz=*/0,
                       /*symmetric=*/false);

  Omega_h::Mesh src_mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0,
                                              0.0, src_n, src_n, 0, false);

  REQUIRE(src_mesh.dim() == 2);
  REQUIRE(tgt_mesh.dim() == 2);

  // Build intersections
  auto intersection = intersectTargets(src_mesh, tgt_mesh);

  auto src_coords = src_mesh.coords();
  auto tgt_coords = tgt_mesh.coords();

  // TEST 1: CONSTANT FIELD (Galerkin projection must match exactly)
  SECTION("Galerkin projection of constant field")
  {

    double constant_value = 1.0;

    // Source field is constant everywhere
    Omega_h::Write<Omega_h::Real> src_values(src_mesh.nverts(), constant_value);

    // Projection
    auto tgt_values = pcms::solveGalerkinProjection(tgt_mesh, src_mesh,
                                                    intersection, src_values);

    REQUIRE(tgt_values.size() == tgt_mesh.nverts());

    auto tgt_vals = Omega_h::HostRead(tgt_values);

    for (int v = 0; v < tgt_mesh.nverts(); v++) {
      CHECK_THAT(tgt_vals[v], Catch::Matchers::WithinAbs(constant_value, 5e-5));
    }
  }

  // TEST 2: LINEAR FIELD  f(x,y) = x + y  (Galerkin projection must reproduce
  // exactly)
  SECTION("Galerkin projection of linear field f(x,y)=x+y")
  {

    Omega_h::Write<Omega_h::Real> src_values(src_mesh.nverts());

    // Assign linear field on source
    Omega_h::parallel_for(
      src_mesh.nverts(), OMEGA_H_LAMBDA(int v) {
        double x = src_coords[2 * v + 0];
        double y = src_coords[2 * v + 1];
        src_values[v] = x + y;
      });

    // Projection
    auto tgt_values = pcms::solveGalerkinProjection(tgt_mesh, src_mesh,
                                                    intersection, src_values);

    auto tgt_vals = Omega_h::HostRead(tgt_values);
    auto tgt_coords = Omega_h::HostRead(tgt_mesh.coords());

    double tol = 5e-5;

    for (int v = 0; v < tgt_mesh.nverts(); v++) {
      double x = tgt_coords[2 * v + 0];
      double y = tgt_coords[2 * v + 1];
      double exact = x + y;

      CHECK_THAT(tgt_vals[v], Catch::Matchers::WithinAbs(exact, tol));
    }
  }
}
