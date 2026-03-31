
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/conservative_projection_solver.hpp>
#include <pcms/transfer/mesh_intersection.hpp>
#include <pcms/transfer/load_vector_integrator.hpp>

TEST_CASE("mesh intersection error evaluator",
          "[transfer][mesh_intersection][error_eval]")
{
  Omega_h::Library lib;

  Omega_h::Reals coords({
    0.0, 0.0, // v0
    1.0, 0.0, // v1
    1.0, 1.0, // v2
    0.0, 1.0  // v3
  });

  // Source mesh: diagonal (v0-v2)
  Omega_h::LOs ev2v_source({0, 1, 2, 0, 2, 3});
  Omega_h::Mesh source_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&source_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_source, coords);

  // Target mesh: opposite diagonal (v1-v3)
  Omega_h::LOs ev2v_target({0, 1, 3, 1, 2, 3});
  Omega_h::Mesh target_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&target_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_target, coords);

  const auto src_coords = source_mesh.coords();
  const auto tgt_coords = target_mesh.coords();

  auto intersections = pcms::intersectTargets(source_mesh, target_mesh);

  SECTION("constant field gives zero projection and conservation error")
  {
    const double c = 2.0;

    Omega_h::Write<Omega_h::Real> source_const(source_mesh.nverts());
    Omega_h::Write<Omega_h::Real> target_const(target_mesh.nverts());

    Omega_h::parallel_for(
      source_mesh.nverts(), OMEGA_H_LAMBDA(int i) { source_const[i] = c; });

    Omega_h::parallel_for(
      target_mesh.nverts(), OMEGA_H_LAMBDA(int i) { target_const[i] = c; });

    auto errs = pcms::evaluate_proj_and_cons_errors(
      target_mesh, source_mesh, intersections, target_const, source_const);

    REQUIRE(errs.proj_err == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(errs.cons_err == Catch::Approx(0.0).margin(1e-12));
  }

  SECTION("linear field gives zero projection and conservation error")
  {
    Omega_h::Write<Omega_h::Real> source_linear(source_mesh.nverts());
    Omega_h::Write<Omega_h::Real> target_linear(target_mesh.nverts());

    // u(x,y) = x + y
    Omega_h::parallel_for(
      source_mesh.nverts(), OMEGA_H_LAMBDA(int i) {
        const double x = src_coords[2 * i + 0];
        const double y = src_coords[2 * i + 1];
        source_linear[i] = x + y;
      });

    Omega_h::parallel_for(
      target_mesh.nverts(), OMEGA_H_LAMBDA(int i) {
        const double x = tgt_coords[2 * i + 0];
        const double y = tgt_coords[2 * i + 1];
        target_linear[i] = x + y;
      });

    auto errs = pcms::evaluate_proj_and_cons_errors(
      target_mesh, source_mesh, intersections, target_linear, source_linear);

    REQUIRE(errs.proj_err == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(errs.cons_err == Catch::Approx(0.0).margin(1e-12));
  }

  SECTION("different linear fields give nonzero error")
  {
    Omega_h::Write<Omega_h::Real> source_field(source_mesh.nverts());
    Omega_h::Write<Omega_h::Real> target_field(target_mesh.nverts());

    // source: u_s = x + y
    Omega_h::parallel_for(
      source_mesh.nverts(), OMEGA_H_LAMBDA(int i) {
        const double x = src_coords[2 * i + 0];
        const double y = src_coords[2 * i + 1];
        source_field[i] = x + y;
      });

    // target: u_t = 2x - y
    Omega_h::parallel_for(
      target_mesh.nverts(), OMEGA_H_LAMBDA(int i) {
        const double x = tgt_coords[2 * i + 0];
        const double y = tgt_coords[2 * i + 1];
        target_field[i] = 2.0 * x - y;
      });

    auto errs = pcms::evaluate_proj_and_cons_errors(
      target_mesh, source_mesh, intersections, target_field, source_field);

    REQUIRE(errs.proj_err > 0.0);
    REQUIRE(errs.cons_err >= 0.0);

  }
}
