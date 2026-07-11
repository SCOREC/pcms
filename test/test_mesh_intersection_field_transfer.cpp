#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/conservative_projection_solver.hpp>
#include <pcms/transfer/mesh_intersection.hpp>
#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/field/function_space/lagrange.h>
#include "field_test_utils.h"

namespace
{

double integrate_linear_field(Omega_h::Mesh& mesh, const Omega_h::Reals& u)
{
  const auto elem_areas = Omega_h::measure_elements_real(&mesh);
  const auto elem_verts = mesh.ask_elem_verts();
  const auto nverts = mesh.nverts();

  REQUIRE(static_cast<Omega_h::LO>(u.size()) == nverts);

  const auto elem_areas_h = Omega_h::HostRead<Omega_h::Real>(elem_areas);
  const auto elem_verts_h = Omega_h::HostRead<Omega_h::LO>(elem_verts);
  const auto u_h = Omega_h::HostRead<Omega_h::Real>(u);

  double integral = 0.0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const Omega_h::LO v0 = elem_verts_h[3 * e + 0];
    const Omega_h::LO v1 = elem_verts_h[3 * e + 1];
    const Omega_h::LO v2 = elem_verts_h[3 * e + 2];

    const double area = elem_areas_h[e];
    const double avg = (u_h[v0] + u_h[v1] + u_h[v2]) / 3.0;
    integral += area * avg;
  }
  return integral;
}

Omega_h::Reals make_omega_h_reals(
  pcms::Rank1View<const pcms::Real, pcms::HostMemorySpace> values)
{
  Omega_h::HostWrite<Omega_h::Real> values_h(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    values_h[i] = values[i];
  }
  return Omega_h::Reals(values_h);
}

} // namespace

TEST_CASE("mesh intersection linear/constant conservation",
          "[transfer][mesh_intersection]")
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

  auto intersections = pcms::intersectTargets(source_mesh, target_mesh);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 2.0;
    Omega_h::Write<Omega_h::Real> source_const(source_mesh.nverts());
    Omega_h::parallel_for(
      source_const.size(), OMEGA_H_LAMBDA(int i) { source_const[i] = c; });

    auto projected = pcms::solveGalerkinProjection(target_mesh, source_mesh,
                                                   intersections, source_const);
    auto projected_h = Omega_h::HostRead<Omega_h::Real>(projected);

    REQUIRE(static_cast<Omega_h::LO>(projected.size()) == target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(projected_h[i] == Catch::Approx(c).margin(1e-10));
    }

    const double source_integral =
      integrate_linear_field(source_mesh, source_const);
    const double target_integral =
      integrate_linear_field(target_mesh, projected);
    REQUIRE(target_integral == Catch::Approx(source_integral).margin(1e-10));
  }

  SECTION("linear field is reproduced on target vertices")
  {
    Omega_h::Write<Omega_h::Real> source_linear(source_mesh.nverts());
    Omega_h::parallel_for(
      source_mesh.nverts(), OMEGA_H_LAMBDA(int i) {
        const double x = src_coords[2 * i + 0];
        const double y = src_coords[2 * i + 1];
        source_linear[i] = x + y;
      });

    auto projected = pcms::solveGalerkinProjection(
      target_mesh, source_mesh, intersections, source_linear);
    auto projected_h = Omega_h::HostRead<Omega_h::Real>(projected);
    const auto tgt_coords_h =
      Omega_h::HostRead<Omega_h::Real>(target_mesh.coords());

    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = tgt_coords_h[2 * i + 0] + tgt_coords_h[2 * i + 1];
      REQUIRE(projected_h[i] == Catch::Approx(expected).margin(1e-9));
    }
  }
}

TEST_CASE("OmegaHConservativeProjection matches conservative projection solver",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Reals coords({0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0});

  Omega_h::LOs ev2v_source({0, 1, 2, 0, 2, 3});
  Omega_h::Mesh source_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&source_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_source, coords);

  Omega_h::LOs ev2v_target({0, 1, 3, 1, 2, 3});
  Omega_h::Mesh target_mesh(&lib);
  Omega_h::build_from_elems_and_coords(&target_mesh, OMEGA_H_SIMPLEX, 2,
                                       ev2v_target, coords);

  // Add classification for all dimensions
  // Vertices (dim 0)
  source_mesh.add_tag<Omega_h::I8>(
    0, "class_dim", 1,
    Omega_h::Read<Omega_h::I8>(source_mesh.nverts(), Omega_h::I8(0)));
  source_mesh.add_tag<Omega_h::ClassId>(
    0, "class_id", 1,
    Omega_h::Read<Omega_h::ClassId>(source_mesh.nverts(), Omega_h::ClassId(0)));

  // Edges (dim 1)
  source_mesh.add_tag<Omega_h::I8>(
    1, "class_dim", 1,
    Omega_h::Read<Omega_h::I8>(source_mesh.nedges(), Omega_h::I8(1)));
  source_mesh.add_tag<Omega_h::ClassId>(
    1, "class_id", 1,
    Omega_h::Read<Omega_h::ClassId>(source_mesh.nedges(), Omega_h::ClassId(0)));

  // Faces (dim 2)
  source_mesh.add_tag<Omega_h::I8>(
    2, "class_dim", 1,
    Omega_h::Read<Omega_h::I8>(source_mesh.nelems(), Omega_h::I8(2)));
  source_mesh.add_tag<Omega_h::ClassId>(
    2, "class_id", 1,
    Omega_h::Read<Omega_h::ClassId>(source_mesh.nelems(), Omega_h::ClassId(0)));

  // Add classification for all dimensions
  // Vertices (dim 0)
  target_mesh.add_tag<Omega_h::I8>(
    0, "class_dim", 1,
    Omega_h::Read<Omega_h::I8>(target_mesh.nverts(), Omega_h::I8(0)));
  target_mesh.add_tag<Omega_h::ClassId>(
    0, "class_id", 1,
    Omega_h::Read<Omega_h::ClassId>(target_mesh.nverts(), Omega_h::ClassId(0)));

  // Edges (dim 1)
  target_mesh.add_tag<Omega_h::I8>(
    1, "class_dim", 1,
    Omega_h::Read<Omega_h::I8>(target_mesh.nedges(), Omega_h::I8(1)));
  target_mesh.add_tag<Omega_h::ClassId>(
    1, "class_id", 1,
    Omega_h::Read<Omega_h::ClassId>(target_mesh.nedges(), Omega_h::ClassId(0)));

  // Faces (dim 2)
  target_mesh.add_tag<Omega_h::I8>(
    2, "class_dim", 1,
    Omega_h::Read<Omega_h::I8>(target_mesh.nelems(), Omega_h::I8(2)));
  target_mesh.add_tag<Omega_h::ClassId>(
    2, "class_id", 1,
    Omega_h::Read<Omega_h::ClassId>(target_mesh.nelems(), Omega_h::ClassId(0)));

  auto source_space = pcms::LagrangeFunctionSpace::FromMesh(
    source_mesh, 1, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto target_space = pcms::LagrangeFunctionSpace::FromMesh(
    target_mesh, 1, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);

  auto source = source_space.CreateField<pcms::Real>();
  auto target = target_space.CreateField<pcms::Real>();

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return x * x + x * y + 0.5 * y * y;
    });

  const auto intersections = pcms::intersectTargets(source_mesh, target_mesh);
  const auto expected =
    pcms::solveGalerkinProjection(target_mesh, source_mesh, intersections,
                                  make_omega_h_reals(pcms::FlattenToRank1View(
                                    source.GetDOFHolderDataHost())));
  const auto expected_h = Omega_h::HostRead<Omega_h::Real>(expected);

  pcms::OmegaHConservativeProjection projection(source_space, target_space);
  projection.Apply(source, target);

  const auto target_values =
    pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
  REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
          target_mesh.nverts());
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    REQUIRE(target_values[i] == Catch::Approx(expected_h[i]).margin(1e-12));
  }
}
