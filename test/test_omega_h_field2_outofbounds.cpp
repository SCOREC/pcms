#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/function_space/spline.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/uniform_grid.h"
#include "field_test_utils.h"
#include <Kokkos_Core.hpp>
#include <vector>

using pcms::Real;

TEST_CASE("omega_h_field2 out of bounds FILL mode")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, 10, 10, 0, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory->CreateFunction<Real>();
  pcms::test::SetField(
    field.GetData(), *factory->GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + y; });

  Real fill_value = -999.0;

  // Test points - mix of inside and outside
  std::vector<Real> coords = {
    0.5,  0.5,  // inside - should evaluate normally
    1.5,  0.5,  // outside (x > 1) - should return fill_value
    0.5,  -0.1, // outside (y < 0) - should return fill_value
    0.3,  0.7,  // inside - should evaluate normally
    -0.1, 0.5,  // outside (x < 0) - should return fill_value
  };

  std::vector<bool> is_inside = {true, false, false, true, false};
  pcms::test::CheckEvaluationWithFill(
    factory, field, coords, is_inside,
    OMEGA_H_LAMBDA(Real x, Real y) { return x + y; }, fill_value, 1.0e-10);
}

TEST_CASE("uniform_grid_field out of bounds FILL mode")
{
  const int N = 10;
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {N, N};

  auto factory = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::CoordinateSystem::Cartesian, 1);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});
  pcms::test::SetField(
    field.GetData(), *factory.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + y; });

  Real fill_value = -999.0;

  // Test points - mix of inside and outside the [0,1]^2 grid
  std::vector<Real> coords = {
    0.5,  0.5,  // inside - should evaluate normally
    1.5,  0.5,  // outside (x > 1) - should return fill_value
    0.5,  -0.1, // outside (y < 0) - should return fill_value
    0.3,  0.7,  // inside - should evaluate normally
    -0.1, 0.5,  // outside (x < 0) - should return fill_value
  };

  std::vector<bool> is_inside = {true, false, false, true, false};
  pcms::test::CheckEvaluationWithFill(
    factory, field, coords, is_inside,
    OMEGA_H_LAMBDA(Real x, Real y) { return x + y; }, fill_value, 1.0e-8);
}

TEST_CASE("spline_field out of bounds FILL mode")
{
  const int N = 10;
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {N, N};

  auto factory = pcms::SplineFunctionSpace::FromUniformGrid(
    grid, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});
  pcms::test::SetField(
    field.GetData(), *factory.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + y; });

  Real fill_value = -999.0;

  // Test points - mix of inside and outside the [0,1]^2 grid
  std::vector<Real> coords = {
    0.5,  0.5,  // inside - should evaluate normally
    1.5,  0.5,  // outside (x > 1) - should return fill_value
    0.5,  -0.1, // outside (y < 0) - should return fill_value
    0.3,  0.7,  // inside - should evaluate normally
    -0.1, 0.5,  // outside (x < 0) - should return fill_value
  };

  std::vector<bool> is_inside = {true, false, false, true, false};
  pcms::test::CheckEvaluationWithFill(
    factory, field, coords, is_inside,
    OMEGA_H_LAMBDA(Real x, Real y) { return x + y; }, fill_value, 1.0e-8);
}

// Note: polynomial_reconstruction_field does not support out-of-bounds FILL
// mode. The MLS (Moving Least Squares) implementation always attempts
// interpolation/ extrapolation for all points and does not check the
// OutOfBoundsPolicy.
