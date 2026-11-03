#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/field/layout/omega_h_lagrange.h>
#include <pcms/transfer/omega_h_mass_integrator.hpp>
#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include "field_test_utils.h"
#include <petscksp.h>

// Reference corner tet:
//   v0=(0,0,0), v1=(1,0,0), v2=(0,1,0), v3=(0,0,1)
// Volume V = 1/6.
//
// P1 consistent mass on a tet:
//   M_ij = V/20 * (1 + delta_ij)
//        = 1/60 on diagonal, 1/120 off-diagonal.

TEST_CASE("OmegaHMassIntegrator (3D): P1 mass on reference tet matches "
          "analytic V/20*(1+delta)",
          "[mass_integrator][3d][analytic]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildReferenceTet(lib);
  REQUIRE(mesh.dim() == 3);
  REQUIRE(mesh.nverts() == 4);
  REQUIRE(mesh.nelems() == 1);

  constexpr pcms::Real V = 1.0 / 6.0;
  constexpr pcms::Real M_diag = V / 10.0; // 1/60
  constexpr pcms::Real M_off = V / 20.0;  // 1/120

  auto space = pcms::test::MakeP1Space(mesh);
  auto integrator = pcms::BuildOmegaHMassIntegrator(*space);
  Mat mat = integrator->GetMatrix();

  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space->GetLayout());
  REQUIRE(layout != nullptr);

  // PETSc rows are active indices = global_to_local(local_vertex).
  const auto perm = layout->GetGlobalToLocalPermutationHost();
  REQUIRE(static_cast<int>(perm.extent(0)) == 4);

  pcms::Real grand_total = 0.0;
  for (int i = 0; i < 4; ++i) {
    pcms::Real row_sum = 0.0;
    for (int j = 0; j < 4; ++j) {
      PetscInt r = static_cast<PetscInt>(perm(i));
      PetscInt c = static_cast<PetscInt>(perm(j));
      PetscScalar got = 0.0;
      MatGetValues(mat, 1, &r, 1, &c, &got);
      const pcms::Real expected = (i == j) ? M_diag : M_off;
      CAPTURE(i, j, r, c, expected, got);
      CHECK(static_cast<pcms::Real>(got) ==
            Catch::Approx(expected).epsilon(1e-12));
      row_sum += static_cast<pcms::Real>(got);
    }
    // Row sum =\int \lamda_i = V/4
    CHECK(row_sum == Catch::Approx(V / 4.0).epsilon(1e-12));
    grand_total += row_sum;
  }
  // Sum of row sums = volume
  CHECK(grand_total == Catch::Approx(V).epsilon(1e-12));
}

TEST_CASE("OmegaHMassIntegrator (3D): P0 mass on reference tet equals volume",
          "[mass_integrator][3d][analytic]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildReferenceTet(lib);

  constexpr pcms::Real V = 1.0 / 6.0;

  auto space = pcms::test::MakeP0Space(mesh);
  auto integrator = pcms::BuildOmegaHMassIntegrator(*space);
  Mat mat = integrator->GetMatrix();

  PetscInt r = 0;
  PetscInt c = 0;
  PetscScalar got = 0.0;
  MatGetValues(mat, 1, &r, 1, &c, &got);
  CHECK(static_cast<pcms::Real>(got) == Catch::Approx(V).epsilon(1e-12));
}

TEST_CASE("OmegaHConservativeProjection (3D): same-mesh reference tet is "
          "exact for constant and linear fields",
          "[transfer][mesh_intersection][3d][analytic]")
{
  // Purpose: source == target == one tet. Intersection is trivial; P1 must
  // reproduce constants/linears at vertices. If this fails, bug is in
  // RHS/mass/Apply/KSP — not multi-tet clipping.
  Omega_h::Library lib;
  Omega_h::Mesh mesh = pcms::test::BuildReferenceTet(lib);
  REQUIRE(mesh.nverts() == 4);
  REQUIRE(mesh.nelems() == 1);
  auto space = pcms::test::MakeP1Space(mesh);
  auto source = space->CreateFunction<pcms::Real>();
  auto target = space->CreateFunction<pcms::Real>();
  pcms::OmegaHConservativeProjection projection(*space, *space);
  SECTION("constant field")
  {
    const double c = 2.5;
    pcms::test::SetField(
      source, KOKKOS_LAMBDA(pcms::Real, pcms::Real, pcms::Real) { return c; });
    projection.Apply(source, target);
    const auto values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
    REQUIRE(static_cast<Omega_h::LO>(values.size()) == 4);
    for (Omega_h::LO i = 0; i < 4; ++i) {
      CAPTURE(i, values[i]);
      REQUIRE(values[i] == Catch::Approx(c).margin(1e-12));
    }
  }
  SECTION("linear field")
  {
    pcms::test::SetField(
      source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
        return 1.0 + x + 2.0 * y + 3.0 * z;
      });
    projection.Apply(source, target);
    const auto values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
    const auto coords_h = pcms::test::CopyCoordinatesToHost(
      pcms::MakeConstRank2View(mesh.coords(), 3));
    for (Omega_h::LO i = 0; i < 4; ++i) {
      const double expected =
        1.0 + coords_h(i, 0) + 2.0 * coords_h(i, 1) + 3.0 * coords_h(i, 2);
      CAPTURE(i, expected, values[i]);
      REQUIRE(values[i] == Catch::Approx(expected).margin(1e-12));
    }
  }
}
