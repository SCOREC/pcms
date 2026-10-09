#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/transfer/mesh_intersection.hpp>
#include <pcms/transfer/omega_h_intersection_quadrature.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/utility/mesh_geometry.h>
#include "field_test_utils.h"

#include <algorithm>
#include <cmath>

namespace
{

void AddReorderedVertexGlobalIds(Omega_h::Mesh& mesh)
{
  mesh.add_tag<Omega_h::GO>(
    0, "reordered_global", 1,
    Omega_h::Read<Omega_h::GO>({2, 0, 3, 1}, "reordered_global"));
}

void AddSparseVertexGlobalIds(Omega_h::Mesh& mesh)
{
  mesh.add_tag<Omega_h::GO>(
    0, "sparse_global", 1,
    Omega_h::Read<Omega_h::GO>({102, 7, 41, 19}, "sparse_global"));
}

} // namespace

// Fields that are in the target order-1 space (constants and affine
// functions) should be reproduced exactly by the L2 projection, and the
// integral must be conserved.
TEST_CASE("OmegaHConservativeProjection reproduces constant and linear fields",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  // Source and target triangulate the same square along opposite diagonals.
  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 2.0;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(target_values[i] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-10));
  }

  SECTION("linear field is reproduced on target vertices and conserved")
  {
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + y; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    const auto tgt_coords_h = pcms::test::CopyCoordinatesToHost(
      pcms::MakeConstRank2View(target_mesh.coords(), 2), target_mesh.nverts(),
      2);
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = tgt_coords_h(i, 0) + tgt_coords_h(i, 1);
      REQUIRE(target_values[i] == Catch::Approx(expected).margin(1e-9));
    }
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-9));
  }
}

TEST_CASE("OmegaHConservativeProjection conserves the integral",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  const auto mass_type =
    GENERATE(pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped);
  CAPTURE(static_cast<int>(mass_type));

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space,
                                                mass_type);

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return x * x + x * y + 0.5 * y * y;
    });
  projection.Apply(source, target);
  REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
          Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
            .margin(1e-12));

  // Second Apply with a different source field — verifies cached state is
  // reused correctly and is not stale.
  pcms::test::SetField(
    source,
    OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return 2.0 * x - y + 0.5; });
  projection.Apply(source, target);
  REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
          Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
            .margin(1e-12));
}

// Row-sum lumping preserves constants and integral. It also bounds the values.
TEST_CASE("OmegaHConservativeProjection with lumped mass reproduces constants, "
          "not linears, and stays bounded",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space,
                                                pcms::MassMatrixType::Lumped);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 2.0;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(target_values[i] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-10));
  }

  SECTION("linear field conserves the integral but is not nodally exact")
  {
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + y; });

    projection.Apply(source, target);

    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-9));

    // On this mesh the lumped solution at the corner vertices is the
    // phi-weighted patch average of x + y, which is off by 0.5 — far above
    // solver tolerance. A consistent mass matrix reproduces the field to
    // 1e-9 (see the exact-reproduction test above), so a large error here
    // proves the lumped path was actually taken.
    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    const auto tgt_coords_h =
      Omega_h::HostRead<Omega_h::Real>(target_mesh.coords());
    double max_err = 0.0;
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = tgt_coords_h[2 * i + 0] + tgt_coords_h[2 * i + 1];
      max_err = std::max(max_err, std::abs(target_values[i] - expected));
    }
    CHECK(max_err > 1e-3);

    pcms::test::RequireBoundedBy(target, source);
  }

  SECTION("hat function stays inside source field range")
  {
    // A hat centered on v0: nodal values (1, 0, 0, 0).
    // This function can cause overshoot/undershoot with
    // L2 projection
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
        return (1.0 - x) * (1.0 - y);
      });

    projection.Apply(source, target);

    pcms::test::RequireBoundedBy(target, source);
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-9));
  }
}

TEST_CASE("OmegaHConservativeProjection writes reordered target GIDs in local "
          "DOF order",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));
  AddReorderedVertexGlobalIds(target_mesh);

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh, "reordered_global");

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);
  projection.Apply(source, target);

  const auto target_values = target.GetDOFHolderDataHost();
  const auto target_coords =
    target_space->GetLayout()->GetDOFHolderCoordinates().GetValues();
  const auto target_coords_h = pcms::test::CopyCoordinatesToHost(
    target_coords, target_mesh.nverts(), target_mesh.dim());
  REQUIRE(static_cast<Omega_h::LO>(target_values.extent(0)) ==
          target_mesh.nverts());
  REQUIRE(target_values.extent(1) == 1);
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = target_coords_h(i, 0) + 2.0 * target_coords_h(i, 1);
    CAPTURE(i);
    REQUIRE(target_values(i, 0) == Catch::Approx(expected).margin(1e-9));
  }
}

TEST_CASE("OmegaHConservativeProjection maps sparse target GIDs to active "
          "PETSc rows",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));
  AddSparseVertexGlobalIds(target_mesh);

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh, "sparse_global");

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);
  projection.Apply(source, target);

  const auto target_values = target.GetDOFHolderDataHost();
  const auto target_coords =
    target_space->GetLayout()->GetDOFHolderCoordinates().GetValues();
  const auto target_coords_h = pcms::test::CopyCoordinatesToHost(
    target_coords, target_mesh.nverts(), target_mesh.dim());
  REQUIRE(static_cast<Omega_h::LO>(target_values.extent(0)) ==
          target_mesh.nverts());
  REQUIRE(target_values.extent(1) == 1);
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = target_coords_h(i, 0) + 2.0 * target_coords_h(i, 1);
    CAPTURE(i);
    REQUIRE(target_values(i, 0) == Catch::Approx(expected).margin(1e-9));
  }
}

// P0 (piecewise-constant) source -> P1 target. A constant lives in P1, so it is
// reproduced exactly; a non-constant P0 source is not reproduced but its
// integral must still be conserved.
TEST_CASE("OmegaHConservativeProjection P0 source to P1 target",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP0Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 3.5;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(target_values[i] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP0Field(source_mesh, source))
              .margin(1e-10));
  }

  SECTION("non-constant P0 source conserves the integral")
  {
    pcms::test::SetField(
      source,
      OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

    projection.Apply(source, target);

    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP0Field(source_mesh, source))
              .margin(1e-10));
  }
}

// P1 source -> P0 (piecewise-constant) target. A constant is reproduced
// exactly; a linear source projects to its cell average, which for a linear
// function equals the value at each element centroid. The integral is conserved
// in every case.
TEST_CASE("OmegaHConservativeProjection P1 source to P0 target",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP0Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = -1.25;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nelems());
    for (Omega_h::LO e = 0; e < target_mesh.nelems(); ++e) {
      REQUIRE(target_values[e] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(pcms::test::IntegrateP0Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-10));
  }

  SECTION("linear field projects to cell average and conserves the integral")
  {
    const auto f = [](double x, double y) { return x + 2.0 * y; };
    pcms::test::SetField(
      source,
      OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    const auto centroids_h = Omega_h::HostRead<Omega_h::Real>(
      pcms::get_entity_centroids(target_mesh, target_mesh.dim()));
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nelems());
    for (Omega_h::LO e = 0; e < target_mesh.nelems(); ++e) {
      const double expected = f(centroids_h[2 * e + 0], centroids_h[2 * e + 1]);
      REQUIRE(target_values[e] == Catch::Approx(expected).margin(1e-9));
    }
    REQUIRE(pcms::test::IntegrateP0Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-9));
  }
}

TEST_CASE("OmegaHConservativeProjection reuses a shared mesh intersection",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_p1 = pcms::test::MakeP1Space(target_mesh);
  auto target_p0 = pcms::test::MakeP0Space(target_mesh);

  auto intersection = pcms::IntersectMeshes(*source_space, *target_p1);
  REQUIRE(intersection->GetTargetToSource().tgt2src_offsets.size() ==
          target_mesh.nelems() + 1);
  REQUIRE(intersection->GetSourceDiscretization()->SameEntities(
    *source_space->GetDiscretization()));

  auto source = source_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return x * x + x * y + 0.5 * y * y;
    });

  SECTION("shared and self-built intersections give identical results")
  {
    const auto mass_type =
      GENERATE(pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped);
    CAPTURE(static_cast<int>(mass_type));
    pcms::OmegaHConservativeProjection own(*source_space, *target_p1,
                                           mass_type);
    pcms::OmegaHConservativeProjection shared(*source_space, *target_p1,
                                              mass_type, intersection);
    auto expected = target_p1->CreateFunction<pcms::Real>();
    auto got = target_p1->CreateFunction<pcms::Real>();
    own.Apply(source, expected);
    shared.Apply(source, got);
    const auto e = pcms::FlattenToRank1View(expected.GetDOFHolderDataHost());
    const auto g = pcms::FlattenToRank1View(got.GetDOFHolderDataHost());
    REQUIRE(e.size() == g.size());
    for (pcms::LO i = 0; i < static_cast<pcms::LO>(e.size()); ++i) {
      CAPTURE(i);
      CHECK(g[i] == e[i]);
    }
  }

  SECTION("the same intersection serves a target space of another order")
  {
    pcms::OmegaHConservativeProjection projection(
      *source_space, *target_p0, pcms::MassMatrixType::Consistent,
      intersection);
    auto target = target_p0->CreateFunction<pcms::Real>();
    projection.Apply(source, target);
    REQUIRE(pcms::test::IntegrateP0Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-12));
  }

  SECTION("an intersection of other meshes is rejected")
  {
    auto reversed = pcms::IntersectMeshes(*target_p1, *source_space);
    CHECK_THROWS_AS(
      pcms::OmegaHConservativeProjection(
        *source_space, *target_p1, pcms::MassMatrixType::Consistent, reversed),
      pcms::pcms_error);
  }
}

TEST_CASE(
  "OmegaHConservativeProjection reuses a shared intersection quadrature",
  "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    pcms::test::BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_p1 = pcms::test::MakeP1Space(target_mesh);
  auto target_p0 = pcms::test::MakeP0Space(target_mesh);

  auto intersection = pcms::IntersectMeshes(*source_space, *target_p1);
  auto quadrature = std::make_shared<pcms::OmegaHIntersectionQuadrature>(
    *source_space, *target_p1, intersection);
  REQUIRE(quadrature->GetNumPoints() > 0);
  REQUIRE(quadrature->GetSourceEvaluator() != nullptr);
  REQUIRE(static_cast<pcms::LO>(quadrature->GetSourceElementIds().extent(0)) ==
          quadrature->GetNumPoints());

  auto source = source_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return x * x + x * y + 0.5 * y * y;
    });

  SECTION("shared quadrature gives the self-built result for both mass types")
  {
    const auto mass_type =
      GENERATE(pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped);
    CAPTURE(static_cast<int>(mass_type));
    pcms::OmegaHConservativeProjection own(*source_space, *target_p1,
                                           mass_type);
    pcms::OmegaHConservativeProjection shared(*source_space, *target_p1,
                                              mass_type, nullptr, quadrature);
    auto expected = target_p1->CreateFunction<pcms::Real>();
    auto got = target_p1->CreateFunction<pcms::Real>();
    own.Apply(source, expected);
    shared.Apply(source, got);
    const auto e = pcms::FlattenToRank1View(expected.GetDOFHolderDataHost());
    const auto g = pcms::FlattenToRank1View(got.GetDOFHolderDataHost());
    REQUIRE(e.size() == g.size());
    for (pcms::LO i = 0; i < static_cast<pcms::LO>(e.size()); ++i) {
      CAPTURE(i);
      CHECK(g[i] == Catch::Approx(e[i]).margin(1e-14));
    }
  }

  SECTION("the search-free evaluator matches the searched one")
  {
    auto searched = source_space->CreatePointEvaluator<pcms::Real>(
      pcms::EvaluationRequest::FromCoordinates(
        quadrature->GetIntegrationPoints()));
    const auto n = quadrature->GetNumPoints();
    Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> direct("direct", n, 1);
    Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> reference("reference",
                                                                  n, 1);
    quadrature->GetSourceEvaluator()->Evaluate(source,
                                               pcms::MakeRank2View(direct));
    searched->Evaluate(source, pcms::MakeRank2View(reference));
    auto direct_h =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, direct);
    auto reference_h =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, reference);
    for (pcms::LO i = 0; i < n; ++i) {
      CAPTURE(i);
      CHECK(direct_h(i, 0) == Catch::Approx(reference_h(i, 0)).margin(1e-13));
    }
  }

  SECTION("a quadrature for other spaces is rejected")
  {
    auto other = std::make_shared<pcms::OmegaHIntersectionQuadrature>(
      *source_space, *target_p0, intersection);
    CHECK_THROWS_AS(pcms::OmegaHConservativeProjection(
                      *source_space, *target_p1,
                      pcms::MassMatrixType::Consistent, nullptr, other),
                    pcms::pcms_error);
  }

  SECTION("element ids outside the mesh are rejected")
  {
    const auto n = quadrature->GetNumPoints();
    Kokkos::View<pcms::LO*, pcms::DeviceMemorySpace> bad("bad", n);
    Kokkos::deep_copy(bad, source_mesh.nelems());
    CHECK_THROWS_AS(source_space->CreatePointEvaluator<pcms::Real>(
                      pcms::EvaluationRequest::FromElements(
                        quadrature->GetIntegrationPoints(), bad)),
                    pcms::pcms_error);
  }
}
