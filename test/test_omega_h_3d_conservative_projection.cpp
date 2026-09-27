#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/transfer/omega_h_control_variate_projection.hpp>
#include <pcms/transfer/monte_carlo_sampling.hpp>
#include <pcms/transfer/copy.h>
#include <pcms/transfer/interpolator.h>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/utility/arrays.h>
#include "field_test_utils.h"

#include <vector>

TEST_CASE("OmegaHConservativeProjection (3D tets) reproduces constant and "
          "linear fields",
          "[transfer][mesh_intersection][3d]")
{
  Omega_h::Library lib;

  // Two independent tessellations of the same unit cube.
  Omega_h::Mesh source_mesh = pcms::test::BuildUnitCube(lib, 1);
  Omega_h::Mesh target_mesh = pcms::test::BuildUnitCube(lib, 2);

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 2.0;
    pcms::test::SetField(
      source, KOKKOS_LAMBDA(pcms::Real, pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(target_values[i] == Catch::Approx(c).margin(1e-9));
    }
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-9));
  }

  SECTION("linear field is reproduced on target vertices and conserved")
  {
    pcms::test::SetField(
      source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
        return 1.0 + x + 2.0 * y + 3.0 * z;
      });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
    const auto tgt_coords_h = pcms::test::CopyCoordinatesToHost(
      pcms::MakeConstRank2View(target_mesh.coords(), 3));
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = 1.0 + tgt_coords_h(i, 0) +
                              2.0 * tgt_coords_h(i, 1) +
                              3.0 * tgt_coords_h(i, 2);
      REQUIRE(target_values[i] == Catch::Approx(expected).margin(1e-8));
    }
    REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
              .margin(1e-8));
  }
}

TEST_CASE("OmegaHConservativeProjection (3D tets) conserves the integral for a "
          "P0 target",
          "[transfer][mesh_intersection][3d]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh = pcms::test::BuildUnitCube(lib, 2);
  Omega_h::Mesh target_mesh = pcms::test::BuildUnitCube(lib, 1);

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP0Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  pcms::test::SetField(
    source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
      return 1.0 + x + 2.0 * y + 3.0 * z;
    });

  projection.Apply(source, target);

  REQUIRE(pcms::test::IntegrateP0Field(target_mesh, target) ==
          Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
            .margin(1e-8));
}

TEST_CASE("Copy transfer (3D tets) reproduces the source field",
          "[transfer][copy][3d]")
{
  Omega_h::Library lib;
  Omega_h::Mesh mesh = pcms::test::BuildUnitCube(lib, 2);
  auto space = pcms::test::MakeP1Space(mesh);

  auto source = space->CreateFunction<pcms::Real>();
  auto target = space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
      return 1.0 + x + 2.0 * y + 3.0 * z;
    });

  pcms::Copy<pcms::Real> copy(*space, *space);
  copy.Apply(source, target);

  const auto sv =
    pcms::FlattenToRank1View(source.GetDOFHolderDataHost().GetValues());
  const auto tv =
    pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
  REQUIRE(tv.size() == sv.size());
  for (std::size_t i = 0; i < tv.size(); ++i) {
    REQUIRE(tv[i] == Catch::Approx(sv[i]));
  }
}

// Point interpolation evaluates the source field at each target DOF site. A
// linear field is reproduced exactly at the target vertices in 3D.
TEST_CASE("Interpolation transfer (3D tets) reproduces a linear field",
          "[transfer][interpolation][3d]")
{
  Omega_h::Library lib;
  Omega_h::Mesh source_mesh = pcms::test::BuildUnitCube(lib, 2);
  Omega_h::Mesh target_mesh = pcms::test::BuildUnitCube(lib, 3);

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
      return 1.0 + x + 2.0 * y + 3.0 * z;
    });

  pcms::Interpolator<pcms::Real> interp(*source_space, *target_space);
  interp.Apply(source, target);

  const auto target_values =
    pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
  const auto tgt_coords_h = pcms::test::CopyCoordinatesToHost(
    pcms::MakeConstRank2View(target_mesh.coords(), 3));
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = 1.0 + tgt_coords_h(i, 0) +
                            2.0 * tgt_coords_h(i, 1) + 3.0 * tgt_coords_h(i, 2);
    REQUIRE(target_values[i] == Catch::Approx(expected).margin(1e-8));
  }
}

// The Monte-Carlo/control-variate projection uses the source field interpolated
// onto the target space as a control variate, so a field already representable
// in the target P1 space (an affine function) is reproduced exactly and its
// integral conserved, even with very few stochastic samples.
TEST_CASE("OmegaHControlVariateProjection (3D tets) is exact for target-space "
          "fields",
          "[transfer][monte_carlo][3d]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh = pcms::test::BuildUnitCube(lib, 1);
  Omega_h::Mesh target_mesh = pcms::test::BuildUnitCube(lib, 2);

  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
      return 1.0 + x + 2.0 * y + 3.0 * z;
    });

  pcms::OmegaHControlVariateProjection projection(
    *source_space, *target_space, /*samples_per_element=*/8,
    pcms::MonteCarloSampling::UniformRandom, /*seed=*/12345);
  projection.Apply(source, target);

  const auto target_values =
    pcms::FlattenToRank1View(target.GetDOFHolderDataHost().GetValues());
  const auto tgt_coords_h = pcms::test::CopyCoordinatesToHost(
    pcms::MakeConstRank2View(target_mesh.coords(), 3));
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = 1.0 + tgt_coords_h(i, 0) +
                            2.0 * tgt_coords_h(i, 1) + 3.0 * tgt_coords_h(i, 2);
    REQUIRE(target_values[i] == Catch::Approx(expected).margin(1e-8));
  }
  REQUIRE(pcms::test::IntegrateP1Field(target_mesh, target) ==
          Catch::Approx(pcms::test::IntegrateP1Field(source_mesh, source))
            .margin(1e-8));
}
