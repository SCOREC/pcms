#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/field/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_space/fusion_io.hpp>
#include <pcms/field/function_space/fusion_io_source.hpp>
#include <pcms/utility/eqdsk.h>
#include "field_test_utils.h"
#include <algorithm>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using pcms::DeviceMemorySpace;
using pcms::FusionIOFormat;
using pcms::FusionIOOptions;
using pcms::HostMemorySpace;
using pcms::Real;
using pcms::ValueBasis;

namespace
{

// The Fusion-IO example data (fusion-io/examples/data) in the test data.
const std::string GeqdskFile =
  std::string(PCMS_TEST_FUSIONIO_DIR) + "/geqdsk/g158115.04701";
const std::string M3DC1File =
  std::string(PCMS_TEST_FUSIONIO_DIR) + "/m3dc1/C1.h5";

// Fusion-IO's poloidal flux field is the total flux 2 pi R A_phi; geqdsk
// files and Fusion-IO's flux series are per radian.
constexpr Real TwoPi = 2.0 * std::numbers::pi;

// Values of `field` at `points`, (R, phi, Z) row-major, copied to the host.
auto EvaluateAt(const pcms::FunctionSpace& space,
                const pcms::Field<Real>& field, const std::vector<Real>& points)
{
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    points, space.GetCoordinateSystem(), 3);
  const auto evaluator = space.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(coords.coordinate_view));
  Kokkos::View<Real**, DeviceMemorySpace> out(
    "out", points.size() / 3, field.GetLayout().GetNumComponents());
  evaluator->Evaluate(field, pcms::test::TagLike(field, out));
  return Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out);
}

Real EvaluateSeries(fio_source& source, int type)
{
  fio_series* opened = nullptr;
  const int status = source.get_series(type, &opened);
  const std::unique_ptr<fio_series> series(opened);
  REQUIRE(status == FIO_SUCCESS);
  Real value = 0.0;
  REQUIRE(series->eval(0.0, &value) == FIO_SUCCESS);
  return value;
}

} // namespace

TEST_CASE("Fusion-IO geqdsk poloidal flux matches the grid values")
{
  auto source = pcms::OpenFusionIOSource(FusionIOFormat::GEQDSK, GeqdskFile);
  auto space = pcms::FusionIOFunctionSpace::Create(source, 1);
  auto psi = space->CreateFunction("psi", "poloidal flux");

  // Fusion-IO's interpolant reproduces the grid values (read here by PCMS's
  // own reader) at the grid nodes, and phi must not matter.
  const pcms::EQDSKData eqdsk(GeqdskFile);
  std::vector<Real> points;
  std::vector<Real> expected;
  for (int i : {eqdsk.GetNW() / 3, eqdsk.GetNW() / 2}) {
    for (int j : {eqdsk.GetNH() / 3, eqdsk.GetNH() / 2}) {
      for (Real phi : {0.0, 1.3}) {
        points.insert(points.end(), {eqdsk.GetR(i), phi, eqdsk.GetZ(j)});
        expected.push_back(TwoPi * eqdsk.GetPsi(i, j));
      }
    }
  }
  const auto result = EvaluateAt(*space, psi, points);
  for (size_t i = 0; i < expected.size(); ++i) {
    CAPTURE(i, points[3 * i], points[3 * i + 2]);
    REQUIRE_THAT(result(i, 0), WithinRel(expected[i], pcms::test::EvalTol));
  }
}

TEST_CASE("Fusion-IO M3D-C1 fields")
{
  auto source = pcms::OpenFusionIOSource(FusionIOFormat::M3DC1, M3DC1File);
  auto scalars = pcms::FusionIOFunctionSpace::Create(source, 1);
  auto vectors = pcms::FusionIOFunctionSpace::Create(source, 3);
  auto psi = scalars->CreateFunction("psi", "poloidal flux");
  auto a = vectors->CreateFunction("A", "vector potential");
  REQUIRE(pcms::SameValueBasis(
    a.GetData().GetValueBasis(),
    ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
               pcms::ComponentScaling::Physical, pcms::values::Vector}));

  // The evaluator stores each point's search hint as sizeof_search_hint()
  // zero-initialized bytes (an element guess for M3D-C1). Revisit it if this
  // changes.
  REQUIRE(source->Get().sizeof_search_hint() == sizeof(int));

  SECTION("the flux at the magnetic axis is the axis flux")
  {
    const Real r = EvaluateSeries(source->Get(), FIO_MAGAXIS_R);
    const Real z = EvaluateSeries(source->Get(), FIO_MAGAXIS_Z);
    const auto result = EvaluateAt(*scalars, psi, {r, 0.0, z});
    REQUIRE_THAT(
      result(0, 0),
      WithinRel(TwoPi * EvaluateSeries(source->Get(), FIO_MAGAXIS_PSI), 1e-3));
  }

  SECTION("vector components are in (R, phi, Z) order: psi = 2 pi R A_phi")
  {
    const std::vector<Real> points = {1.6, 0.1, 0.2, 2.0, 2.5, -0.4};
    const auto psi_values = EvaluateAt(*scalars, psi, points);
    const auto a_values = EvaluateAt(*vectors, a, points);
    for (size_t i = 0; i < psi_values.extent(0); ++i) {
      CAPTURE(i);
      REQUIRE_THAT(TwoPi * points[3 * i] * a_values(i, 1),
                   WithinRel(psi_values(i, 0), pcms::test::EvalTol));
    }
  }

  SECTION("(R, phi, Z) is right-handed: B_R = -dA_phi/dZ")
  {
    // The axisymmetric part has no phi derivative, so this holds exactly for
    // a right-handed system and with the opposite sign otherwise.
    const FusionIOOptions equilibrium{.part =
                                        FusionIOOptions::Part::Equilibrium};
    auto a_eq = vectors->CreateFunction("A", "vector potential", equilibrium);
    auto b_eq = vectors->CreateFunction("B", "magnetic field", equilibrium);
    const Real r = 2.0, z = 0.3, dz = 1e-4;
    const auto a_values =
      EvaluateAt(*vectors, a_eq, {r, 0.0, z + dz, r, 0.0, z - dz});
    const auto b_values = EvaluateAt(*vectors, b_eq, {r, 0.0, z});
    const Real dadz = (a_values(0, 1) - a_values(1, 1)) / (2.0 * dz);
    // The central difference is accurate to O(dz^2), far below dz.
    REQUIRE_THAT(b_values(0, 0), WithinRel(-dadz, dz));
  }

  SECTION("options reach Fusion-IO: the equilibrium part is axisymmetric")
  {
    auto psi_eq = scalars->CreateFunction(
      "psi_eq", "poloidal flux", {.part = FusionIOOptions::Part::Equilibrium});
    const std::vector<Real> points = {1.6, 0.1, 0.2, 1.6, 2.5, 0.2};
    const auto total = EvaluateAt(*scalars, psi, points);
    const auto equilibrium = EvaluateAt(*scalars, psi_eq, points);
    REQUIRE_THAT(equilibrium(1, 0),
                 WithinRel(equilibrium(0, 0), pcms::test::ExactTol));
    // The perturbation in the total varies with phi.
    REQUIRE_FALSE(total(1, 0) == total(0, 0));
  }

  SECTION("values match direct fusion-io evaluation")
  {
    const FusionIOOptions options{.timeslice = 1,
                                  .part = FusionIOOptions::Part::Perturbed};
    fio_option_list fio_options;
    source->Get().get_field_options(&fio_options);
    fio_options.set_option(FIO_TIMESLICE, 1);
    fio_options.set_option(FIO_PART, FIO_PERTURBED_ONLY);
    const std::vector<Real> points = {1.3, 0.0, 0.4, 1.6, 0.7, -0.6,
                                      2.0, 2.5, 0.1, 2.2, 4.0, -0.3};
    for (const auto& [quantity, space] :
         {std::pair{"poloidal flux", scalars},
          std::pair{"magnetic field", vectors}}) {
      CAPTURE(quantity);
      const auto values = EvaluateAt(
        *space, space->CreateFunction("f", quantity, options), points);
      const auto field =
        source->GetField(source->Quantities().at(quantity), fio_options);
      for (size_t i = 0; i < values.extent(0); ++i) {
        Real expected[3];
        REQUIRE(field->eval(&points[3 * i], expected) == FIO_SUCCESS);
        for (int c = 0; c < field->dimension(); ++c) {
          CAPTURE(i, c);
          // The absolute floor covers components near zero.
          REQUIRE_THAT(values(i, c),
                       WithinRel(expected[c], pcms::test::ExactTol) ||
                         WithinAbs(expected[c], pcms::test::ExactTol));
        }
      }
    }
  }
}

namespace
{

// A source reporting Cartesian geometry; nothing else is used.
class CartesianSource final : public fio_source
{
public:
  int open(const char*) override { return FIO_SUCCESS; }
  int close() override { return FIO_SUCCESS; }
  int get_field_options(fio_option_list*) const override { return FIO_SUCCESS; }
  int get_field(const field_type, fio_field**, const fio_option_list*) override
  {
    return FIO_UNSUPPORTED;
  }
  int get_available_fields(fio_field_list*) const override
  {
    return FIO_SUCCESS;
  }
  int get_int_parameter(const parameter_type t, int* p) const override
  {
    if (t != FIO_GEOMETRY) {
      return FIO_UNSUPPORTED;
    }
    *p = FIO_CARTESIAN;
    return FIO_SUCCESS;
  }
};

} // namespace

TEST_CASE("Fusion-IO evaluation errors")
{
  auto source = pcms::OpenFusionIOSource(FusionIOFormat::GEQDSK, GeqdskFile);
  auto scalars = pcms::FusionIOFunctionSpace::Create(source, 1);
  auto vectors = pcms::FusionIOFunctionSpace::Create(source, 3);
  auto psi = scalars->CreateFunction("psi", "poloidal flux");
  const auto quantities = pcms::AvailableQuantities(*source);
  REQUIRE(std::find(quantities.begin(), quantities.end(), "magnetic field") !=
          quantities.end());

  REQUIRE_THROWS_WITH(EvaluateAt(*scalars, psi, {50.0, 0.0, 0.0}),
                      ContainsSubstring("outside the source's domain"));
  REQUIRE_THROWS_WITH(scalars->CreateFunction("B", "magnetic field"),
                      ContainsSubstring("the space holds 1"));
  REQUIRE_THROWS_WITH(vectors->CreateFunction("psi", "poloidal flux"),
                      ContainsSubstring("the space holds 3"));
  REQUIRE_THROWS_WITH(
    scalars->CreateFunction("psi", "poloidal fluxx"),
    ContainsSubstring("\"poloidal fluxx\" is not a quantity") &&
      ContainsSubstring("poloidal flux,"));
  REQUIRE_THROWS_WITH(psi.GetDOFHolderData(),
                      ContainsSubstring("no DOFAccess"));
  REQUIRE_THROWS_WITH(
    pcms::OpenFusionIOSource(FusionIOFormat::GEQDSK, "missing"),
    ContainsSubstring("could not open"));
  REQUIRE_THROWS_WITH(
    pcms::FusionIOSource(pcms::FusionIOSourcePtr(new CartesianSource)),
    ContainsSubstring("cylindrical geometry"));

  // Fusion-IO writes one row per query point.
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    {1.8, 0.0, 0.1, 2.0, 0.0, 0.1}, scalars->GetCoordinateSystem(), 3);
  const auto evaluator = scalars->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(coords.coordinate_view));
  Kokkos::View<Real**, DeviceMemorySpace> short_out("short_out", 1, 1);
  REQUIRE_THROWS_WITH(
    evaluator->Evaluate(psi, pcms::test::TagLike(psi, short_out)),
    ContainsSubstring("but expected [2][1]"));

  auto other = pcms::FusionIOFunctionSpace::Create(
    pcms::OpenFusionIOSource(FusionIOFormat::GEQDSK, GeqdskFile), 1);
  REQUIRE_THROWS_WITH(EvaluateAt(*other, psi, {1.8, 0.0, 0.1}),
                      ContainsSubstring("field is not a quantity"));
}
