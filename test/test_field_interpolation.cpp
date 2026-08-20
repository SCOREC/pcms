#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include "pcms/transfer/interpolator.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/assert.h"
#include "field_test_utils.h"
#include <Kokkos_Core.hpp>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"

using pcms::Real;

KOKKOS_INLINE_FUNCTION static Real interpolation_linear_f(Real x, Real y)
{
  return -0.3 * x + 0.5 * y;
}

TEST_CASE("interpolate linear 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred());
  auto field = factory->CreateFunction<Real>();
  auto interpolated = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field, OMEGA_H_LAMBDA(Real x, Real y) { return -0.3 * x + 0.5 * y; });

  pcms::Interpolator<Real> interp(*factory, *factory);
  interp.Apply(field, interpolated);
  auto interpolated_dof =
    pcms::FlattenToRank1View(interpolated.GetDOFHolderDataHost());
  auto original_dof =
    pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  REQUIRE(interpolated_dof.size() == original_dof.size());
  for (int i = 0; i < static_cast<int>(interpolated_dof.size()); ++i) {
    REQUIRE_THAT(interpolated_dof[i],
                 Catch::Matchers::WithinRel(original_dof[i], 0.001) ||
                   Catch::Matchers::WithinAbs(original_dof[i], 1E-10));
  }
}

#ifdef PCMS_ENABLE_MESHFIELDS
TEST_CASE("interpolate quadratic 2d meshfields_field")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto factory2 = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 2, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);
  auto field = factory2->CreateFunction<Real>();
  auto interpolated = factory2->CreateFunction<Real>();
  pcms::test::SetField(
    field,
    OMEGA_H_LAMBDA(Real x, Real y) { return interpolation_linear_f(x, y); });

  pcms::Interpolator<Real> interp(*factory2, *factory2);
  interp.Apply(field, interpolated);

  auto interpolated_dof =
    pcms::FlattenToRank1View(interpolated.GetDOFHolderDataHost());
  auto original_dof =
    pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  REQUIRE(interpolated_dof.size() == original_dof.size());
  for (int i = 0; i < static_cast<int>(interpolated_dof.size()); ++i) {
    REQUIRE_THAT(interpolated_dof[i],
                 Catch::Matchers::WithinRel(original_dof[i], 0.001) ||
                   Catch::Matchers::WithinAbs(original_dof[i], 1E-10));
  }
}
#endif

TEST_CASE("interpolate quadratic 2d omega_h_field throws")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);

  REQUIRE_THROWS_AS(pcms::LagrangeFunctionSpace::FromMesh(
                      mesh, 2, 1, pcms::csys::Cartesian::Deferred(), "global",
                      pcms::LagrangeFunctionSpace::Backend::OmegaH),
                    pcms::pcms_error);
}

// Interpolator<T> test: construct once (localize), apply twice with different
// source data, verify correct results each time. This demonstrates that the
// localization only happens once (at construction) while Apply can be called
// cheaply in a coupling loop.
TEST_CASE("Interpolator: construct once, apply twice with different data")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 100, 100, 0, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred());

  auto source = factory->CreateFunction<Real>();
  auto target = factory->CreateFunction<Real>();

  // Construct Interpolator once — localization happens here.
  // Because source and target share the same layout (same factory),
  // the target DOF holder coordinates are the source mesh node coordinates.
  pcms::Interpolator<pcms::Real> interp(*factory, *factory);

  // First application: linear function f(x,y) = -0.3*x + 0.5*y
  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(Real x, Real y) { return -0.3 * x + 0.5 * y; });
  interp.Apply(source, target);

  {
    auto src_dof =
      pcms::FlattenToRank1View(source.GetDOFHolderDataHost());
    auto tgt_dof =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(tgt_dof.size() == src_dof.size());
    for (int i = 0; i < static_cast<int>(tgt_dof.size()); ++i) {
      REQUIRE_THAT(tgt_dof[i], Catch::Matchers::WithinRel(src_dof[i], 0.001) ||
                                 Catch::Matchers::WithinAbs(src_dof[i], 1E-10));
    }
  }

  // Second application: constant function f(x,y) = 7.0
  // Apply is called without re-constructing the interpolator (no
  // re-localization).
  pcms::test::SetField(source, OMEGA_H_LAMBDA(Real, Real) { return 7.0; });
  interp.Apply(source, target);

  {
    auto tgt_dof =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    for (int i = 0; i < static_cast<int>(tgt_dof.size()); ++i) {
      REQUIRE_THAT(tgt_dof[i], Catch::Matchers::WithinAbs(7.0, 1E-10));
    }
  }
}
