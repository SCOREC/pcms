#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/field/evaluation_request.h>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/field/coordinate_systems/cartesian.hpp>
#include <pcms/field/layout/omega_h_lagrange.h>
#include <pcms/field/transformed_source.hpp>
#include <pcms/transfer/conservative_projection_solver.hpp>
#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/transfer/omega_h_control_variate_projection.hpp>
#include <pcms/transfer/omega_h_mass_integrator.hpp>
#include <pcms/transfer/omega_h_mc_rhs_integrator.hpp>
#include "field_test_utils.h"
#include <petscksp.h>
#include <vector>

// ---------------------------------------------------------------------------
// Monte Carlo RHS integrator
// ---------------------------------------------------------------------------

TEST_CASE("OmegaHMonteCarloRHSIntegrator: sample points lie inside the domain",
          "[mc_rhs_integrator]")
{
  Omega_h::Library lib;
  auto target_mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  const int samples_per_element = 16;
  for (const auto sampling : {pcms::MonteCarloSampling::UniformRandom}) {
    pcms::OmegaHMonteCarloRHSIntegrator integrator(
      *target_space, samples_per_element, sampling);
    const auto raw_coords = integrator.GetIntegrationPoints().GetValues();
    auto coords_h = pcms::test::CopyCoordinatesToHost(raw_coords);

    REQUIRE(
      coords_h.extent(0) ==
      static_cast<std::size_t>(target_mesh.nelems() * samples_per_element));
    for (std::size_t i = 0; i < coords_h.extent(0); ++i) {
      CAPTURE(i, coords_h(i, 0), coords_h(i, 1));
      CHECK(coords_h(i, 0) >= -1e-12);
      CHECK(coords_h(i, 0) <= 1.0 + 1e-12);
      CHECK(coords_h(i, 1) >= -1e-12);
      CHECK(coords_h(i, 1) <= 1.0 + 1e-12);
    }
  }
}

TEST_CASE("OmegaHMonteCarloRHSIntegrator: constant field integrates exactly",
          "[mc_rhs_integrator]")
{
  // For f = c the estimator sums to c * |domain| for any sample placement,
  // because the P1 basis functions partition unity at every sample point.
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitSquare(lib, 1);
  auto target_mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source_field = source_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source_field, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return 2.0; });

  for (const auto sampling : {pcms::MonteCarloSampling::UniformRandom}) {
    pcms::OmegaHMonteCarloRHSIntegrator integrator(*target_space, 8, sampling);
    pcms::test::EvaluateAndAssemble(integrator, source_space, source_field);

    PetscScalar sum = 0.0;
    VecSum(integrator.GetVector(), &sum);
    INFO("sampling=" << "UniformRandom");
    CHECK(static_cast<pcms::Real>(sum) == Catch::Approx(2.0).margin(1e-12));
  }
}

// ---------------------------------------------------------------------------
// Control variate projection
// ---------------------------------------------------------------------------

TEST_CASE("OmegaHControlVariateProjection: exact for fields in the target "
          "space even with few samples",
          "[mc_rhs_integrator][control_variate]")
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

  pcms::OmegaHControlVariateProjection projection(
    *source_space, *target_space, /*samples_per_element=*/4,
    pcms::MonteCarloSampling::UniformRandom, /*seed=*/8675309, mass_type);
  projection.Apply(source_field, target_field);

  const auto values =
    pcms::FlattenToRank1View(target_field.GetDOFHolderDataHost().GetValues());
  const auto coords_h = pcms::test::CopyCoordinatesToHost(
    pcms::MakeConstRank2View(target_mesh.coords(), 2));
  REQUIRE(static_cast<Omega_h::LO>(values.size()) == target_mesh.nverts());
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = 3.0 * coords_h(i, 0) - coords_h(i, 1) + 0.25;
    CAPTURE(i, expected, values[i]);
    CHECK(values[i] == Catch::Approx(expected).margin(1e-9));
  }
}

TEST_CASE("OmegaHControlVariateProjection: reduces error vs plain Monte Carlo",
          "[mc_rhs_integrator][control_variate]")
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
  pcms::OmegaHConservativeProjection reference_projection(*source_space,
                                                          *target_space);
  reference_projection.Apply(source_field, reference_field);
  const auto reference = pcms::FlattenToRank1View(
    reference_field.GetDOFHolderDataHost().GetValues());

  const int samples_per_element = 64;
  const uint64_t seed = 20240611;

  // Plain Monte Carlo projection.
  pcms::OmegaHMonteCarloRHSIntegrator mc_integrator(
    *target_space, samples_per_element, pcms::MonteCarloSampling::UniformRandom,
    seed);
  auto evaluator = source_space->CreatePointEvaluator<pcms::Real>(
    pcms::EvaluationRequest::FromCoordinates(
      mc_integrator.GetIntegrationPoints()));
  pcms::OmegaHMassIntegrator mass_integrator(*target_space);
  pcms::GalerkinProjectionSolver solver(mass_integrator, mc_integrator);
  const auto mc_projected = solver.Solve(*evaluator, source_field);
  const auto mc_projected_h = Omega_h::HostRead<Omega_h::Real>(mc_projected);

  // Control-variate projection with the same sampling configuration.
  auto cv_field = target_space->CreateFunction<pcms::Real>();
  pcms::OmegaHControlVariateProjection cv_projection(
    *source_space, *target_space, samples_per_element,
    pcms::MonteCarloSampling::UniformRandom, seed);
  cv_projection.Apply(source_field, cv_field);
  const auto cv_values =
    pcms::FlattenToRank1View(cv_field.GetDOFHolderDataHost().GetValues());

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
// 3D Monte Carlo RHS integrator
// ---------------------------------------------------------------------------

TEST_CASE("OmegaHMonteCarloRHSIntegrator (3D): sample points lie inside the "
          "unit cube",
          "[mc_rhs_integrator][3d]")
{
  // for  Dim=3 sample coordinates inside [0,1]^3 with correct count.
  Omega_h::Library lib;
  auto target_mesh = pcms::test::BuildUnitCube(lib, 1);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  const int samples_per_element = 16;
  pcms::OmegaHMonteCarloRHSIntegrator integrator(
    *target_space, samples_per_element,
    pcms::MonteCarloSampling::UniformRandom);

  const auto raw_coords = integrator.GetIntegrationPoints().GetValues();
  REQUIRE(raw_coords.extent(1) == 3);

  auto coords_h = pcms::test::CopyCoordinatesToHost(
    raw_coords, static_cast<int>(raw_coords.extent(0)),
    static_cast<int>(raw_coords.extent(1)));

  REQUIRE(coords_h.extent(0) ==
          static_cast<std::size_t>(target_mesh.nelems() * samples_per_element));

  for (std::size_t i = 0; i < coords_h.extent(0); ++i) {
    CAPTURE(i, coords_h(i, 0), coords_h(i, 1), coords_h(i, 2));
    CHECK(coords_h(i, 0) >= -1e-12);
    CHECK(coords_h(i, 0) <= 1.0 + 1e-12);
    CHECK(coords_h(i, 1) >= -1e-12);
    CHECK(coords_h(i, 1) <= 1.0 + 1e-12);
    CHECK(coords_h(i, 2) >= -1e-12);
    CHECK(coords_h(i, 2) <= 1.0 + 1e-12);
  }
}

TEST_CASE("OmegaHMonteCarloRHSIntegrator (3D): constant field integrates "
          "exactly",
          "[mc_rhs_integrator][3d]")
{
  // if f=c => VecSum(b) = c * volume = c * 1 (P1 partition of unity).
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildUnitCube(lib, 1);
  auto target_mesh = pcms::test::BuildUnitCube(lib, 1);
  auto source_space = pcms::test::MakeP1Space(source_mesh);
  auto target_space = pcms::test::MakeP1Space(target_mesh);

  auto source_field = source_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source_field,
    KOKKOS_LAMBDA(pcms::Real, pcms::Real, pcms::Real) { return 2.0; });

  pcms::OmegaHMonteCarloRHSIntegrator integrator(
    *target_space, /*samples_per_element=*/8,
    pcms::MonteCarloSampling::UniformRandom);

  pcms::test::EvaluateAndAssemble(integrator, source_space, source_field);

  PetscScalar sum = 0.0;
  VecSum(integrator.GetVector(), &sum);
  CHECK(static_cast<pcms::Real>(sum) == Catch::Approx(2.0).margin(1e-12));
}

TEST_CASE(
  "OmegaHControlVariateProjection (3D): reduces error vs plain Monte Carlo",
  "[mc_rhs_integrator][control_variate][3d]")
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
  pcms::OmegaHConservativeProjection ref_proj(*source_space, *target_space);
  ref_proj.Apply(source, reference);
  const auto ref_vals =
    pcms::FlattenToRank1View(reference.GetDOFHolderDataHost());

  const int nsample = 64;
  const uint64_t seed = 20240611;

  pcms::OmegaHMonteCarloRHSIntegrator mc(
    *target_space, nsample, pcms::MonteCarloSampling::UniformRandom, seed);
  auto ev = source_space->CreatePointEvaluator<pcms::Real>(
    pcms::EvaluationRequest::FromCoordinates(mc.GetIntegrationPoints()));
  pcms::OmegaHMassIntegrator mass(*target_space);
  pcms::GalerkinProjectionSolver solver(mass, mc);
  const auto mc_x = solver.Solve(*ev, source);
  const auto mc_h = Omega_h::HostRead<Omega_h::Real>(mc_x);

  auto cv_field = target_space->CreateFunction<pcms::Real>();
  pcms::OmegaHControlVariateProjection cv(
    *source_space, *target_space, nsample,
    pcms::MonteCarloSampling::UniformRandom, seed);
  cv.Apply(source, cv_field);
  const auto cv_vals =
    pcms::FlattenToRank1View(cv_field.GetDOFHolderDataHost());

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

class SwapXY final : public pcms::CoordinateMap
{
public:
  [[nodiscard]] std::shared_ptr<const pcms::CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept override
  {
    return pcms::csys::Cartesian::Create(2);
  }
  [[nodiscard]] std::shared_ptr<const pcms::CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept override
  {
    return pcms::csys::Cartesian::Create(2);
  }
  [[nodiscard]] std::unique_ptr<pcms::BoundCoordinateMap> Bind(
    const pcms::CoordinateView<pcms::DeviceMemorySpace>& query_points)
    const override
  {
    const auto in = query_points.GetValues();
    Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> out("swapped",
                                                            in.extent(0), 2);
    Kokkos::parallel_for(
      "swap_xy",
      Kokkos::RangePolicy<pcms::DeviceMemorySpace::execution_space>(
        0, static_cast<pcms::LO>(in.extent(0))),
      KOKKOS_LAMBDA(const pcms::LO i) {
        out(i, 0) = in(i, 1);
        out(i, 1) = in(i, 0);
      });
    Kokkos::fence();
    return std::make_unique<pcms::BoundCoordinateMap>(
      GetTargetCoordinateSystem(), out);
  }
};

} // namespace

TEST_CASE("OmegaHControlVariateProjection over a TransformedSource maps the "
          "sample points",
          "[mc_rhs_integrator][control_variate]")
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

  pcms::TransformedSource src_as_swapped(source_space,
                                         std::make_shared<SwapXY>());
  pcms::OmegaHControlVariateProjection projection(
    src_as_swapped, *target_space, /*samples_per_element=*/4,
    pcms::MonteCarloSampling::UniformRandom);
  projection.Apply(source_field, target_field);

  const auto values =
    pcms::FlattenToRank1View(target_field.GetDOFHolderDataHost().GetValues());
  const auto coords_h = pcms::test::CopyCoordinatesToHost(
    pcms::MakeConstRank2View(target_mesh.coords(), 2));
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = 3.0 * coords_h(i, 1) - coords_h(i, 0) + 0.25;
    CAPTURE(i, expected, values[i]);
    CHECK(values[i] == Catch::Approx(expected).margin(1e-9));
  }
}
