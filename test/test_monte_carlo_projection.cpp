#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/configuration.hpp>
#include <pcms/field/field_metadata.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/backends/omega_h/lagrange_layout.hpp>
#include "field_test_utils.hpp"
#include <vector>
#include "pcms/backends/omega_h/omega_h.hpp"

// ---------------------------------------------------------------------------
// Control variate projection
// ---------------------------------------------------------------------------

TEST_CASE("CreateControlVariateProjection: exact for fields in the target "
          "space even with few samples",
          "[monte_carlo][control_variate]")
{
  // For a globally linear field the control variate (target-space interpolant
  // of the source field) equals the source field, so the sampled residual is
  // identically zero and the projection is exact regardless of sample count.
  // That holds for any invertible mass matrix (delta = M^-1 0 = 0), so the
  // lumped variant must be exact here too.
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitSquare(lib, 1);
  auto target_mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source_field = source_space->CreateFunction<pcms::Real>();
  auto target_field = target_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source_field,
    OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return 3.0 * x - y + 0.25; });

  const auto mass_type =
    GENERATE(pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped);
  CAPTURE(static_cast<int>(mass_type));

  auto projection = pcms::CreateControlVariateProjection(
    source_space, target_space,
    {.samples_per_element = 4, .mass_matrix_type = mass_type});
  projection->Apply(source_field, target_field);

  const auto values = pcms::detail::FlattenToRank1View(
    target_field.GetDOFHolderDataHost().GetValues());
  const auto coords_h = pcms::test::CopyCoordinatesToHost(
    pcms::omega_h::MakeConstRank2View(target_mesh.coords(), 2));
  REQUIRE(static_cast<Omega_h::LO>(values.size()) == target_mesh.nverts());
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = 3.0 * coords_h(i, 0) - coords_h(i, 1) + 0.25;
    CAPTURE(i, expected, values[i]);
    CHECK(values[i] == Catch::Approx(expected).margin(1e-9));
  }
}

TEST_CASE("CreateControlVariateProjection: reduces error vs plain Monte Carlo",
          "[monte_carlo][control_variate]")
{
  // For a non-linear field the control variate does not cancel f exactly, but
  // the residual is much smaller than f, so with identical sample counts and
  // seeds the control-variate projection must be closer to the exact
  // (intersection-quadrature) projection than the plain Monte Carlo one.
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitSquare(lib, 1);
  auto target_mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source_field = source_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source_field, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return x * x + x * y + 0.5 * y * y;
    });

  // Reference: exact mesh-intersection quadrature projection.
  auto reference_field = target_space->CreateFunction<pcms::Real>();
  auto reference_projection =
    pcms::CreateConservativeProjection(source_space, target_space);
  reference_projection->Apply(source_field, reference_field);
  const auto reference = pcms::detail::FlattenToRank1View(
    reference_field.GetDOFHolderDataHost().GetValues());

  const int samples_per_element = 64;
  const uint64_t seed = 20240611;

  // Plain Monte Carlo projection.
  pcms::ProjectionOptions mc_options;
  mc_options.quadrature = pcms::CreateMonteCarloQuadrature(
    target_space->GetGeometry(), samples_per_element, seed);
  auto mc_field = target_space->CreateFunction<pcms::Real>();
  pcms::CreateProjection(source_space, target_space, mc_options)
    ->Apply(source_field, mc_field);
  const auto mc_projected_h = pcms::detail::FlattenToRank1View(
    mc_field.GetDOFHolderDataHost().GetValues());

  // Control-variate projection with the same sampling configuration.
  auto cv_field = target_space->CreateFunction<pcms::Real>();
  auto cv_projection = pcms::CreateControlVariateProjection(
    source_space, target_space,
    {.samples_per_element = samples_per_element, .seed = seed});
  cv_projection->Apply(source_field, cv_field);
  const auto cv_values = pcms::detail::FlattenToRank1View(
    cv_field.GetDOFHolderDataHost().GetValues());

  double mc_error = 0.0;
  double cv_error = 0.0;
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    mc_error = std::max(mc_error, std::abs(mc_projected_h[i] - reference[i]));
    cv_error = std::max(cv_error, std::abs(cv_values[i] - reference[i]));
  }
  CAPTURE(mc_error, cv_error);
  CHECK(cv_error < mc_error);
}

// ---------------------------------------------------------------------------
// 3D
// ---------------------------------------------------------------------------

TEST_CASE(
  "CreateControlVariateProjection (3D): reduces error vs plain Monte Carlo",
  "[monte_carlo][control_variate][3d]")
{
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitCube(lib, 1);
  auto target_mesh = pcms::test::BuildUnitCube(lib, 2);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source, KOKKOS_LAMBDA(pcms::Real x, pcms::Real y, pcms::Real z) {
      return x * x + y * y + z * z;
    });

  auto reference = target_space->CreateFunction<pcms::Real>();
  auto ref_proj =
    pcms::CreateConservativeProjection(source_space, target_space);
  ref_proj->Apply(source, reference);
  const auto ref_vals = pcms::detail::FlattenToRank1View(
    reference.GetDOFHolderDataHost().GetValues());

  const int nsample = 64;
  const uint64_t seed = 20240611;

  pcms::ProjectionOptions mc_options;
  mc_options.quadrature = pcms::CreateMonteCarloQuadrature(
    target_space->GetGeometry(), nsample, seed);
  auto mc_field = target_space->CreateFunction<pcms::Real>();
  pcms::CreateProjection(source_space, target_space, mc_options)
    ->Apply(source, mc_field);
  const auto mc_h = pcms::detail::FlattenToRank1View(
    mc_field.GetDOFHolderDataHost().GetValues());

  auto cv_field = target_space->CreateFunction<pcms::Real>();
  auto cv = pcms::CreateControlVariateProjection(
    source_space, target_space, {.samples_per_element = nsample, .seed = seed});
  cv->Apply(source, cv_field);
  const auto cv_vals = pcms::detail::FlattenToRank1View(
    cv_field.GetDOFHolderDataHost().GetValues());

  double mc_err = 0.0, cv_err = 0.0;
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    mc_err = std::max(mc_err, std::abs(mc_h[i] - ref_vals[i]));
    cv_err = std::max(cv_err, std::abs(cv_vals[i] - ref_vals[i]));
  }
  CAPTURE(mc_err, cv_err);
  CHECK(cv_err < mc_err);
}

namespace
{

// (x, y) -> (y, x)
struct SwapXY
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<pcms::Real, 2>& x,
                                  Kokkos::Array<pcms::Real, 2>& out,
                                  Kokkos::Array<pcms::Real, 4>& jacobian) const
  {
    out[0] = x[1];
    out[1] = x[0];
    jacobian[0] = 0;
    jacobian[1] = 1;
    jacobian[2] = 1;
    jacobian[3] = 0;
  }
};

} // namespace

TEST_CASE("CreateControlVariateProjection through a coordinate map maps the "
          "sample points",
          "[monte_carlo][control_variate]")
{
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitSquare(lib, 1);
  auto target_mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source_field = source_space->CreateFunction<pcms::Real>();
  auto target_field = target_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source_field,
    OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return 3.0 * x - y + 0.25; });

  auto projection = pcms::CreateControlVariateProjection(
    source_space, target_space,
    {.samples_per_element = 4,
     .to_source = pcms::CreateCoordinateMap<2, 2>(
       pcms::csys::CreateCartesian(2), pcms::csys::CreateCartesian(2),
       SwapXY{})});
  projection->Apply(source_field, target_field);

  const auto values = pcms::detail::FlattenToRank1View(
    target_field.GetDOFHolderDataHost().GetValues());
  const auto coords_h = pcms::test::CopyCoordinatesToHost(
    pcms::omega_h::MakeConstRank2View(target_mesh.coords(), 2));
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = 3.0 * coords_h(i, 1) - coords_h(i, 0) + 0.25;
    CAPTURE(i, expected, values[i]);
    CHECK(values[i] == Catch::Approx(expected).margin(1e-9));
  }
}

TEST_CASE("the conservative and control variate projections check units",
          "[projection][metadata]")
{
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitSquare(lib, 1);
  auto target_mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);
  const auto conservative =
    pcms::CreateConservativeProjection(source_space, target_space);
  const auto control_variate =
    pcms::CreateControlVariateProjection(source_space, target_space);

  auto source = source_space->CreateFunction<pcms::Real>(
    "source", {.units = pcms::Units::Parse("m^-3")});
  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + y; });
  const auto target_with = [&](const char* units) {
    pcms::ValueBasis basis;
    if (units != nullptr) {
      basis.units = pcms::Units::Parse(units);
    }
    return target_space->CreateFunction<pcms::Real>("target", basis);
  };
  for (const auto* projection : {conservative.get(), control_variate.get()}) {
    auto same = target_with("1/m^3");
    REQUIRE_NOTHROW(projection->Apply(source, same));
    auto unknown = target_with(nullptr);
    REQUIRE_NOTHROW(projection->Apply(source, unknown));
#ifdef PCMS_ENABLE_UNITS
    auto scaled = target_with("cm^-3");
    REQUIRE_THROWS_WITH(
      projection->Apply(source, scaled),
      Catch::Matchers::ContainsSubstring("differ in scale or offset"));
    auto kelvin = target_with("K");
    REQUIRE_THROWS_WITH(
      projection->Apply(source, kelvin),
      Catch::Matchers::ContainsSubstring("different dimensions"));
#endif
  }
}
