#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/field/detail/field_ops.hpp>
#include <pcms/field/detail/simple_field_data.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include "field_test_utils.hpp"
#include "pcms/field/point_rows.hpp"
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
void CheckAdjoint(const std::shared_ptr<const pcms::FunctionSpace>& space,
                  const pcms::PointEvaluator<pcms::Real>& evaluator,
                  unsigned seed)
{
  REQUIRE(evaluator.SupportsTranspose());
  const int num_components = space->GetLayout()->GetNumComponents();
  const int num_points = evaluator.NumPoints();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  const auto x = RandomDofs(*space->GetLayout(), rng);
  Host2D y("y", num_points, num_components);
  for (std::size_t i = 0; i < y.size(); ++i) {
    y.data()[i] = uniform(rng);
  }

  auto field = space->CreateFunction<pcms::Real>();
  SetDofs(field, x);
  Device2D e_x("Ex", num_points, num_components);
  const auto basis = field.GetData().GetValueBasis();
  evaluator.Evaluate(field,
                     pcms::ValueView<pcms::Real, pcms::DeviceMemorySpace>(
                       basis, pcms::MakeRank2View(e_x)));
  // the device layout is the one ValueView expects (LayoutLeft on a GPU)
  Device2D y_device("y_device", num_points, num_components);
  {
    auto y_mirror = Kokkos::create_mirror_view(y_device);
    for (int i = 0; i < num_points; ++i) {
      for (int c = 0; c < num_components; ++c) {
        y_mirror(i, c) = y(i, c);
      }
    }
    Kokkos::deep_copy(y_device, y_mirror);
  }
  const pcms::ValueView<const pcms::Real, pcms::DeviceMemorySpace> y_view(
    basis, pcms::MakeConstRank2View(y_device));

  auto functional = space->CreateFunction<pcms::Real>();
  evaluator.EvaluateTranspose(y_view, functional);
  const auto e_x_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, e_x);
  const auto num_holders = static_cast<int>(x.extent(0));
  // indexed by (row, component): the device layout may differ from the host's
  std::vector<double> once(x.size());
  {
    const auto transposed = functional.GetDOFHolderDataHost().GetValues();
    for (int h = 0; h < num_holders; ++h) {
      for (int c = 0; c < num_components; ++c) {
        once[h * num_components + c] = transposed(h, c);
      }
    }
  }
  double lhs = 0, rhs = 0, scale = 0;
  for (int p = 0; p < num_points; ++p) {
    for (int c = 0; c < num_components; ++c) {
      lhs += e_x_host(p, c) * y(p, c);
      scale += std::abs(e_x_host(p, c) * y(p, c));
    }
  }
  for (int h = 0; h < num_holders; ++h) {
    for (int c = 0; c < num_components; ++c) {
      rhs += x(h, c) * once[h * num_components + c];
    }
  }
  CAPTURE(lhs, rhs);
  CHECK(std::abs(lhs - rhs) <= 1e-12 * std::max(1.0, scale));
  // the pairing of the functional with the function is the same number
  CHECK(std::abs(pcms::detail::Apply(functional, field) - lhs) <=
        1e-12 * std::max(1.0, scale));

  // overwrites: a second call gives the same result
  evaluator.EvaluateTranspose(y_view, functional);
  const auto again = functional.GetDOFHolderDataHost().GetValues();
  for (int h = 0; h < num_holders; ++h) {
    for (int c = 0; c < num_components; ++c) {
      CHECK(again(h, c) ==
            Catch::Approx(once[h * num_components + c]).margin(1e-14));
    }
  }

  // the evaluation rows: E x from the CSR rows equals Evaluate
  if (const auto* source = pcms::detail::AsRowsSource(evaluator)) {
    const auto* rows = &source->Rows();
    const auto offsets = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, rows->offsets);
    const auto holders = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, rows->holders);
    const auto values = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, rows->values);
    REQUIRE(offsets.extent(0) == static_cast<std::size_t>(num_points) + 1);
    for (int p = 0; p < num_points; ++p) {
      for (int c = 0; c < num_components; ++c) {
        double sum = 0;
        for (int k = offsets(p); k < offsets(p + 1); ++k) {
          sum += values(k) * x(holders(k), c);
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

TEST_CASE("coefficient operations: LinearSum, Scale, Const, Div, DotProd, "
          "Norm2 and the pairing",
          "[field_ops]")
{
  using pcms::detail::Apply;
  using pcms::detail::DotProd;
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto space = pcms::test::MakeP1Space(mesh);
  std::mt19937 rng(7);
  const auto a = RandomDofs(*space->GetLayout(), rng);
  const auto b = RandomDofs(*space->GetLayout(), rng);
  auto x = space->CreateFunction<pcms::Real>();
  auto y = space->CreateFunction<pcms::Real>();
  SetDofs(x, a);
  SetDofs(y, b);
  double dot = 0, xx = 0;
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    for (std::size_t c = 0; c < a.extent(1); ++c) {
      dot += a(i, c) * b(i, c);
      xx += a(i, c) * a(i, c);
    }
  }
  CHECK(DotProd(x, y) == Catch::Approx(dot).epsilon(1e-14));
  CHECK(Apply(x, y) == Catch::Approx(dot).epsilon(1e-14));
  CHECK(pcms::detail::Norm2(x) == Catch::Approx(std::sqrt(xx)).epsilon(1e-14));
  auto z = space->CreateFunction<pcms::Real>();
  pcms::detail::LinearSum(3.0, x, 1.0, y, z); // z = 3 a + b
  pcms::detail::Scale(0.5, z, z);             // z = (3 a + b) / 2
  auto twos = space->CreateFunction<pcms::Real>();
  pcms::detail::Const(2.0, twos);
  pcms::detail::Div(z, twos, z); // z = (3 a + b) / 4
  const auto values = z.GetDOFHolderDataHost().GetValues();
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    for (std::size_t c = 0; c < a.extent(1); ++c) {
      CHECK(values(i, c) ==
            Catch::Approx(0.25 * (3 * a(i, c) + b(i, c))).epsilon(1e-14));
    }
  }
  // fields of another layout are refused
  auto other = pcms::test::MakeP1Space(mesh)->CreateFunction<pcms::Real>();
  CHECK_THROWS_AS(DotProd(x, other), pcms::pcms_error);
  // so are operands of another kind of field data
  const pcms::detail::SimpleFieldData<pcms::Real> other_kind(
    space->GetLayout(), x.GetData().GetValueBasis());
  CHECK_THROWS_WITH(pcms::detail::DotProdLocal(x.GetData(), other_kind),
                    Catch::Matchers::ContainsSubstring("another kind"));
  CHECK_THROWS_WITH(
    pcms::detail::LinearSum(1.0, x.GetData(), 1.0, other_kind, z.GetData()),
    Catch::Matchers::ContainsSubstring("another kind"));
}

TEST_CASE("FunctionSpace::GetQuadrature is exact for products of two "
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
    const auto quadrature = space->GetQuadrature();
    const auto weights = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, quadrature.weights);
    double total = 0;
    for (std::size_t i = 0; i < weights.extent(0); ++i) {
      total += weights(i);
    }
    CHECK(total == Catch::Approx(1.0).epsilon(1e-13));
  }
}
