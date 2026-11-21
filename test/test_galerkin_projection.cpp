#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include "field_test_utils.hpp"
#include <cmath>
#include <map>
#include <random>
#include <vector>

namespace
{

using Host2D = Kokkos::View<pcms::Real**, pcms::HostMemorySpace>;
using Device2D = Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace>;

void SetDofs(pcms::Field<pcms::Real>& field, const Host2D& x)
{
  field.SetDOFHolderDataHost(
    pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(
      field.GetData().GetValueBasis(), pcms::MakeConstRank2View(x)));
}

Host2D RandomDofs(const pcms::FieldLayout& layout, std::mt19937& rng)
{
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  Host2D x("x", layout.GetNumDofHolder(), layout.GetNumComponents());
  for (std::size_t i = 0; i < x.size(); ++i) {
    x.data()[i] = uniform(rng);
  }
  return x;
}

// The generic adjoint test: <E x, y> == <x, E^T y> for random DOFs x of the
// space and random point values y, for any evaluator that SupportsTranspose.
// The space must be unghosted (serial), so that the holders are the owned ones.
// The values are in the evaluator's output basis (a bind with a coordinate
// map re-expresses vectors), and the fields have the basis `variance` gives.
void CheckAdjoint(const std::shared_ptr<const pcms::FunctionSpace>& space,
                  const pcms::PointEvaluator<pcms::Real>& evaluator,
                  unsigned seed,
                  pcms::VarianceSignature variance = pcms::values::Scalar)
{
  REQUIRE(evaluator.SupportsTranspose());
  const int num_components = space->GetLayout()->GetNumComponents();
  const int num_points = evaluator.NumPoints();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  const auto x = RandomDofs(*space->GetLayout(), rng);

  auto field = space->CreateFunction<pcms::Real>("x", variance);
  SetDofs(field, x);
  const auto basis = evaluator.OutputBasis(field.GetData().GetValueBasis());
  const int num_values =
    basis.Rank() == 0 ? num_components : basis.system->Dimension();
  Host2D y("y", num_points, num_values);
  for (std::size_t i = 0; i < y.size(); ++i) {
    y.data()[i] = uniform(rng);
  }
  Device2D e_x("Ex", num_points, num_values);
  evaluator.Evaluate(field,
                     pcms::ValueView<pcms::Real, pcms::DeviceMemorySpace>(
                       basis, pcms::MakeRank2View(e_x)));
  // the device layout is the one ValueView expects (LayoutLeft on a GPU)
  Device2D y_device("y_device", num_points, num_values);
  {
    auto y_mirror = Kokkos::create_mirror_view(y_device);
    for (int i = 0; i < num_points; ++i) {
      for (int c = 0; c < num_values; ++c) {
        y_mirror(i, c) = y(i, c);
      }
    }
    Kokkos::deep_copy(y_device, y_mirror);
  }
  const pcms::ValueView<const pcms::Real, pcms::DeviceMemorySpace> y_view(
    basis, pcms::MakeConstRank2View(y_device));

  auto functional = space->CreateFunction<pcms::Real>("functional", variance);
  evaluator.EvaluateTranspose(y_view, functional);
  const auto e_x_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, e_x);
  const auto transposed = functional.GetDOFHolderDataHost().GetValues();
  const std::vector<double> once(transposed.data_handle(),
                                 transposed.data_handle() + x.size());
  double lhs = 0, rhs = 0, scale = 0;
  for (std::size_t i = 0; i < y.size(); ++i) {
    lhs += e_x_host.data()[i] * y.data()[i];
    scale += std::abs(e_x_host.data()[i] * y.data()[i]);
  }
  for (std::size_t i = 0; i < x.size(); ++i) {
    rhs += x.data()[i] * once[i];
  }
  CAPTURE(lhs, rhs);
  CHECK(std::abs(lhs - rhs) <= 1e-12 * std::max(1.0, scale));

  // accumulates: a second call doubles the result
  evaluator.EvaluateTranspose(y_view, functional);
  const auto twice = functional.GetDOFHolderDataHost().GetValues();
  for (std::size_t i = 0; i < x.size(); ++i) {
    CHECK(twice.data_handle()[i] == Catch::Approx(2 * once[i]).margin(1e-14));
  }

  // the matrix form: E x from the CSR rows equals Evaluate
  if (evaluator.SupportsMatrix()) {
    const auto matrix = evaluator.GetMatrix();
    const auto offsets = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, matrix.row_offsets);
    const auto columns = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, matrix.columns);
    const auto values = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, matrix.values);
    const auto gids = space->GetLayout()->GetGidsHost();
    std::map<pcms::GO, int> holder_of;
    for (std::size_t h = 0; h < gids.extent(0); ++h) {
      holder_of[gids(h)] = static_cast<int>(h);
    }
    REQUIRE(offsets.extent(0) == static_cast<std::size_t>(num_points) + 1);
    for (int p = 0; p < num_points; ++p) {
      for (int c = 0; c < num_components; ++c) {
        double sum = 0;
        for (int k = offsets(p); k < offsets(p + 1); ++k) {
          sum += values(k) * x(holder_of.at(columns(k)), c);
        }
        CHECK(sum == Catch::Approx(e_x_host(p, c)).margin(1e-13));
      }
    }
  }
}

std::vector<pcms::Real> RandomPoints(int dim, int n, unsigned seed)
{
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uniform(0.01, 0.99);
  std::vector<pcms::Real> points(static_cast<std::size_t>(n) * dim);
  for (auto& p : points) {
    p = uniform(rng);
  }
  return points;
}

} // namespace

TEST_CASE("EvaluateTranspose is the adjoint of Evaluate (Omega_h Lagrange)",
          "[quadrature][adjoint]")
{
  Omega_h::Library lib;
  const int dim = GENERATE(2, 3);
  const int order = GENERATE(0, 1);
  const int num_components = GENERATE(1, 2);
  CAPTURE(dim, order, num_components);
  auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 3)
                       : pcms::test::BuildUnitCube(lib, 2);
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  const auto space = pcms::CreateLagrangeSpace(
    geometry, order, num_components,
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
  const auto bound = pcms::test::Bind(*geometry, RandomPoints(dim, 200, 3));
  const auto evaluator = space->CreatePointEvaluator<pcms::Real>(bound);
  CheckAdjoint(space, *evaluator, 11);
}

TEST_CASE("EvaluateTranspose through a coordinate map is the transpose of "
          "Evaluate",
          "[quadrature][adjoint]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitCube(lib, 2);
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  // query points in (R, theta, Z) whose images are in the cube
  const auto cartesian = RandomPoints(3, 50, 5);
  std::vector<pcms::Real> cylindrical;
  for (std::size_t i = 0; i < cartesian.size(); i += 3) {
    const pcms::Real x = cartesian[i];
    const pcms::Real y = cartesian[i + 1];
    cylindrical.insert(cylindrical.end(),
                       {std::hypot(x, y), std::atan2(y, x), cartesian[i + 2]});
  }
  const auto query = pcms::test::CreateDeviceCoordinateView(
    cylindrical, pcms::csys::CreateCylindricalRThetaZ(), 3);
  const auto bound =
    geometry->Bind(query.coordinate_view, pcms::OutOfBoundsPolicy{},
                   pcms::CreateCylindricalToCartesian());
  SECTION("a scalar field")
  {
    const auto space = pcms::test::MakeP1Space(geometry);
    CheckAdjoint(space, *space->CreatePointEvaluator<pcms::Real>(bound), 13);
  }
  SECTION("a vector field: Cartesian components come out cylindrical")
  {
    const auto space = pcms::CreateLagrangeSpace(
      geometry, 1, 3,
      pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
    const auto evaluator = space->CreatePointEvaluator<pcms::Real>(bound);
    const auto f = space->CreateFunction<pcms::Real>("v", pcms::values::Vector);
    CHECK(evaluator->OutputBasis(f.GetData().GetValueBasis()).system->Kind() ==
          "cylindrical_rthetaz");
    CheckAdjoint(space, *evaluator, 17, pcms::values::Vector);
  }
}

TEST_CASE("EvaluateTranspose checks its arguments", "[adjoint]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 1);
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  const auto space = pcms::test::MakeP1Space(geometry);
  const auto other_space = pcms::test::MakeP1Space(geometry);
  const auto bound = pcms::test::Bind(*geometry, RandomPoints(2, 4, 1));
  const auto evaluator = space->CreatePointEvaluator<pcms::Real>(bound);
  auto functional = space->CreateFunction<pcms::Real>();
  auto foreign = other_space->CreateFunction<pcms::Real>();
  const auto basis = functional.GetData().GetValueBasis();
  Device2D values("v", 4, 1);
  Device2D wrong_shape("d", 3, 1);
  using View = pcms::ValueView<const pcms::Real, pcms::DeviceMemorySpace>;
  CHECK_NOTHROW(evaluator->EvaluateTranspose(
    View(basis, pcms::MakeConstRank2View(values)), functional));
  CHECK_THROWS_AS(
    evaluator->EvaluateTranspose(
      View(basis, pcms::MakeConstRank2View(wrong_shape)), functional),
    std::invalid_argument);
  // a functional of another layout than the evaluator's
  CHECK_THROWS_AS(evaluator->EvaluateTranspose(
                    View(basis, pcms::MakeConstRank2View(values)), foreign),
                  pcms::pcms_error);
}

TEST_CASE("FunctionSpace::CreateQuadrature is exact for products of two "
          "functions of the space",
          "[quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  for (const int order : {0, 1}) {
    auto geometry =
      pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
    const auto space = pcms::CreateLagrangeSpace(
      geometry, order, 1,
      pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
    const auto quadrature = space->CreateQuadrature();
    const auto weights = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, quadrature.weights);
    double total = 0;
    for (std::size_t i = 0; i < weights.extent(0); ++i) {
      total += weights(i);
    }
    CHECK(total == Catch::Approx(1.0).epsilon(1e-13));
  }
}
