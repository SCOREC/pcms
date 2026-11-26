#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/configuration.hpp>
#include <pcms/field/field_ops.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/transfer/forms.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/transfer/l2_inner_product.hpp>
#include <pcms/transfer/linear_solver.hpp>
#include "field_test_utils.hpp"
#ifndef PCMS_ENABLE_PETSC
#include <catch2/matchers/catch_matchers_string.hpp>
#endif
#include <cmath>
#include <cstdio>
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

  auto field_basis = pcms::detail::MakeValueBasis(
    space->GetCoordinateSystem(), variance, nullptr, std::nullopt);
  if (variance.Rank() == 0) {
    field_basis = pcms::test::ComponentsBasis(*space);
  }
  auto field = space->CreateFunction<pcms::Real>("x", field_basis);
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

  auto functional =
    space->CreateFunction<pcms::Real>("functional", field_basis);
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
  CHECK(std::abs(pcms::Apply(functional, field) - lhs) <=
        1e-12 * std::max(1.0, scale));

  // accumulates: a second call doubles the result
  evaluator.EvaluateTranspose(y_view, functional);
  const auto twice = functional.GetDOFHolderDataHost().GetValues();
  for (int h = 0; h < num_holders; ++h) {
    for (int c = 0; c < num_components; ++c) {
      CHECK(twice(h, c) ==
            Catch::Approx(2 * once[h * num_components + c]).margin(1e-14));
    }
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
  auto functional =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  auto foreign = other_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*other_space));
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

TEST_CASE("field operations: Scale, Axpy, Dot, norms and the pairing",
          "[field_ops]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto space = pcms::test::MakeP1Space(mesh);
  std::mt19937 rng(7);
  const auto a = RandomDofs(*space->GetLayout(), rng);
  const auto b = RandomDofs(*space->GetLayout(), rng);
  auto x =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  auto y =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  SetDofs(x, a);
  SetDofs(y, b);
  double dot = 0, xx = 0, inf = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    dot += a.data()[i] * b.data()[i];
    xx += a.data()[i] * a.data()[i];
    inf = std::max(inf, std::abs(a.data()[i]));
  }
  CHECK(pcms::Dot(x, y) == Catch::Approx(dot).epsilon(1e-14));
  CHECK(pcms::Apply(x, y) == Catch::Approx(dot).epsilon(1e-14));
  CHECK(pcms::Norm2(x) == Catch::Approx(std::sqrt(xx)).epsilon(1e-14));
  CHECK(pcms::NormInf(x) == Catch::Approx(inf).epsilon(1e-14));
  auto w =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  Host2D twos("w", a.extent(0), 1);
  Kokkos::deep_copy(twos, 2.0);
  SetDofs(w, twos);
  CHECK(pcms::WeightedRmsNorm(x, w) ==
        Catch::Approx(2 * std::sqrt(xx / a.size())).epsilon(1e-14));
  pcms::Axpy(3.0, x, y); // y = b + 3 a
  pcms::Scale(0.5, y);   // y = (b + 3 a) / 2
  const auto values = y.GetDOFHolderDataHost().GetValues();
  for (std::size_t i = 0; i < a.size(); ++i) {
    CHECK(values.data_handle()[i] ==
          Catch::Approx(0.5 * (b.data()[i] + 3 * a.data()[i])).epsilon(1e-14));
  }
  // fields of another layout are refused
  auto other = pcms::test::MakeP1Space(mesh)->CreateFunction<pcms::Real>();
  CHECK_THROWS_AS(pcms::Dot(x, other), pcms::pcms_error);
}

namespace
{

std::vector<double> FieldValues(const pcms::Field<pcms::Real>& field)
{
  const auto values =
    pcms::detail::FlattenToRank1View(field.GetDOFHolderDataHost().GetValues());
  return std::vector<double>(values.data_handle(),
                             values.data_handle() + values.size());
}

std::shared_ptr<const pcms::omega_h::OmegaHGeometry> MakeGeometry(
  Omega_h::Mesh& mesh)
{
  return std::dynamic_pointer_cast<const pcms::omega_h::OmegaHGeometry>(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()));
}

} // namespace

TEST_CASE("CreateProjection with a reference-rule quadrature reproduces fields "
          "of the target space exactly",
          "[quadrature][projection]")
{
  Omega_h::Library lib;
  const int dim = GENERATE(2, 3);
  const int target_order = GENERATE(0, 1);
  CAPTURE(dim, target_order);
  auto source_mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 5)
                              : pcms::test::BuildUnitCube(lib, 3);
  auto target_mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 3)
                              : pcms::test::BuildUnitCube(lib, 2);
  const auto source_geometry = MakeGeometry(source_mesh);
  const auto target_geometry = MakeGeometry(target_mesh);
  const auto source_space = pcms::test::MakeP1Space(source_geometry);
  const auto target_space = pcms::CreateLagrangeSpace(
    target_geometry, target_order, 1,
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
  auto source = source_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*source_space));
  auto target = target_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*target_space));
  pcms::test::SetFieldComponents(
    source, [](pcms::Real x, pcms::Real y, pcms::Real* out) { out[0] = 2.0; });
  pcms::ProjectionOptions options;
  options.quadrature =
    pcms::CreateGaussQuadrature(target_space->GetGeometry(), 2);
  auto projection = pcms::CreateProjection(source_space, target_space, options);
  projection->Apply(source, target);
  // a constant is in both spaces: the projection is the constant
  for (const double value : FieldValues(target)) {
    CHECK(value == Catch::Approx(2.0).margin(1e-10));
  }
}

TEST_CASE("CreateProjection projects several components at once",
          "[quadrature][projection]")
{
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildSquareMesh(lib, 5);
  auto target_mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto source_geometry = MakeGeometry(source_mesh);
  const auto target_geometry = MakeGeometry(target_mesh);
  const auto options =
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH};
  const auto source_space =
    pcms::CreateLagrangeSpace(source_geometry, 1, 2, options);
  const auto target_space =
    pcms::CreateLagrangeSpace(target_geometry, 1, 2, options);
  auto source = source_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*source_space));
  auto target = target_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*target_space));
  pcms::test::SetFieldComponents(
    source, [](pcms::Real x, pcms::Real y, pcms::Real* out) {
      out[0] = x + y;
      out[1] = 2 * x - y + 0.5;
    });
  auto projection = pcms::CreateProjection(source_space, target_space, {});
  projection->Apply(source, target);
  const auto values = FieldValues(target);
  const auto coords = pcms::test::CopyCoordinatesToHost(
    target_space->GetLayout()->GetDOFHolderCoordinates().GetValues());
  REQUIRE(values.size() == 2 * coords.extent(0));
  for (std::size_t i = 0; i < coords.extent(0); ++i) {
    CHECK(values[2 * i] ==
          Catch::Approx(coords(i, 0) + coords(i, 1)).margin(1e-10));
    CHECK(values[2 * i + 1] ==
          Catch::Approx(2 * coords(i, 0) - coords(i, 1) + 0.5).margin(1e-10));
  }
}

// ---------------------------------------------------------------------------
// The solvers: MatrixFree against Assembled
// ---------------------------------------------------------------------------

#ifdef PCMS_ENABLE_PETSC
namespace
{

// The max difference of the two solves of M x = b for a random functional b,
// and the CG iterations
struct SolveComparison
{
  double max_difference; // relative to the largest entry
  int iterations;
  int unknowns;
};

SolveComparison CompareSolvers(
  const std::shared_ptr<const pcms::FunctionSpace>& space)
{
  using namespace pcms::transfer::detail;
  const BilinearForm form(space);
  LinearSolver matrix_free(form, pcms::LinearSolveOptions{});
  LinearSolver assembled(
    form, pcms::LinearSolveOptions{
            .method = pcms::LinearSolveOptions::Method::Assembled});
  std::mt19937 rng(5);
  auto b =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  SetDofs(b, RandomDofs(*space->GetLayout(), rng));
  auto x_free =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  auto x_direct =
    space->CreateFunction<pcms::Real>("", pcms::test::ComponentsBasis(*space));
  const int iterations = matrix_free.Solve(b, x_free);
  assembled.Solve(b, x_direct);
  const double scale = pcms::NormInf(x_direct);
  pcms::Axpy(-1.0, x_direct, x_free);
  return {pcms::NormInf(x_free) / scale, iterations,
          space->GetLayout()->GetNumDofHolder() *
            space->GetLayout()->GetNumComponents()};
}

} // namespace

TEST_CASE("matrix-free Jacobi-PCG matches the assembled direct solve",
          "[quadrature][solver]")
{
  Omega_h::Library lib;
  for (const int n : {4, 8, 16}) {
    auto mesh = pcms::test::BuildSquareMesh(lib, n);
    const auto result = CompareSolvers(pcms::test::MakeP1Space(mesh));
    CHECK(result.max_difference < 1e-11);
    std::printf("[solver] P1 %dx%d square, %d unknowns: CG iterations = %d, "
                "relative max |x_cg - x_lu| = %.2e\n",
                n, n, result.unknowns, result.iterations,
                result.max_difference);
  }
  {
    auto mesh = pcms::test::BuildSquareMesh(lib, 8);
    const auto result = CompareSolvers(pcms::test::MakeP0Space(mesh));
    CHECK(result.max_difference < 1e-11);
    std::printf("[solver] P0 8x8 square, %d unknowns: CG iterations = %d, "
                "relative max |x_cg - x_lu| = %.2e\n",
                result.unknowns, result.iterations, result.max_difference);
  }
  {
    auto mesh = pcms::test::BuildSquareMesh(lib, 8);
    auto geometry =
      pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
    const auto space = pcms::CreateLagrangeSpace(
      geometry, 1, 3,
      pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
    const auto result = CompareSolvers(space);
    CHECK(result.max_difference < 1e-11);
    std::printf("[solver] P1 8x8 square, 3 components, %d unknowns: CG "
                "iterations = %d, relative max |x_cg - x_lu| = %.2e\n",
                result.unknowns, result.iterations, result.max_difference);
  }
  {
    auto mesh = pcms::test::BuildUnitCube(lib, 3);
    const auto result = CompareSolvers(pcms::test::MakeP1Space(mesh));
    CHECK(result.max_difference < 1e-11);
    std::printf("[solver] P1 3D cube, %d unknowns: CG iterations = %d, "
                "relative max |x_cg - x_lu| = %.2e\n",
                result.unknowns, result.iterations, result.max_difference);
  }
}
#else
TEST_CASE("the assembled solve needs PETSc", "[quadrature][solver]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  const auto space = pcms::test::MakeP1Space(mesh);
  const pcms::transfer::detail::BilinearForm form(space);
  CHECK_THROWS_WITH(pcms::transfer::detail::LinearSolver(
                      form,
                      pcms::LinearSolveOptions{
                        .method = pcms::LinearSolveOptions::Method::Assembled}),
                    Catch::Matchers::ContainsSubstring("PCMS_ENABLE_PETSC"));
}
#endif

TEST_CASE("projection: both solve methods and the lumped mass",
          "[quadrature][projection]")
{
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildSquareMesh(lib, 6);
  auto target_mesh = pcms::test::BuildSquareMesh(lib, 4);
  const auto source_geometry = MakeGeometry(source_mesh);
  const auto source_space = pcms::test::MakeP1Space(source_geometry);
  const auto target_space = pcms::test::MakeP1Space(target_mesh);
  auto source = source_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*source_space));
  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return std::sin(4 * x) * std::cos(3 * y);
    });
  const auto project = [&](pcms::ProjectionOptions options) {
    auto target = target_space->CreateFunction<pcms::Real>(
      "", pcms::test::ComponentsBasis(*target_space));
    pcms::CreateProjection(source_space, target_space, options)
      ->Apply(source, target);
    return target;
  };
  auto free = project({});
  auto lumped = project(
    pcms::ProjectionOptions{.mass_matrix_type = pcms::MassMatrixType::Lumped});
#ifdef PCMS_ENABLE_PETSC
  auto direct = project(pcms::ProjectionOptions{
    .solve = {.method = pcms::LinearSolveOptions::Method::Assembled}});
  auto difference = project({});
  pcms::Axpy(-1.0, direct, difference);
  const double max_difference = pcms::NormInf(difference);
  CHECK(max_difference < 1e-11);
  std::printf("[projection] MatrixFree vs Assembled: max difference %.2e\n",
              max_difference);
#else
  const auto& direct = free;
#endif
  // the lumped projection is another operator, but it has the integral of the
  // consistent one (1^T M_L x = 1^T b = 1^T M x)
  const pcms::L2Integral integral(target_space);
  CHECK(integral.Apply(lumped) ==
        Catch::Approx(integral.Apply(direct)).epsilon(1e-10));
  // MaxIterations too small to converge is an error
  CHECK_THROWS(pcms::CreateProjection(
                 source_space, target_space,
                 pcms::ProjectionOptions{
                   .solve = {.relative_tolerance = 1e-14, .max_iterations = 1}})
                 ->Apply(source, free));
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

// ---------------------------------------------------------------------------
// ApplyTranspose
// ---------------------------------------------------------------------------

namespace
{

// <P a, b> and <a, P^T b> for random coefficients a (of the source space) and
// b (of the target space)
std::pair<double, double> AdjointPair(
  const pcms::TransferOperator<pcms::Real>& op,
  const std::shared_ptr<const pcms::FunctionSpace>& source_space,
  const std::shared_ptr<const pcms::FunctionSpace>& target_space)
{
  REQUIRE(op.SupportsTranspose());
  std::mt19937 rng(21);
  auto a = source_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*source_space));
  auto b = target_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*target_space));
  SetDofs(a, RandomDofs(*source_space->GetLayout(), rng));
  SetDofs(b, RandomDofs(*target_space->GetLayout(), rng));
  auto pa = target_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*target_space));
  op.Apply(a, pa);
  auto ptb = source_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*source_space));
  op.ApplyTranspose(b, ptb);
  return {pcms::Dot(pa, b), pcms::Dot(a, ptb)};
}

} // namespace

TEST_CASE("ApplyTranspose of an interpolator is its adjoint, and preserves "
          "totals",
          "[transfer][adjoint]")
{
  Omega_h::Library lib;
  const int dim = GENERATE(2, 3);
  const int components = GENERATE(1, 2);
  auto source_mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 4)
                              : pcms::test::BuildUnitCube(lib, 3);
  auto target_mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 3)
                              : pcms::test::BuildUnitCube(lib, 2);
  const auto options =
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH};
  const auto source_space = pcms::CreateLagrangeSpace(MakeGeometry(source_mesh),
                                                      1, components, options);
  const auto target_space = pcms::CreateLagrangeSpace(MakeGeometry(target_mesh),
                                                      1, components, options);
  const auto interpolator =
    pcms::CreateInterpolator<pcms::Real>(source_space, target_space);
  const auto [lhs, rhs] =
    AdjointPair(*interpolator, source_space, target_space);
  CHECK(lhs == Catch::Approx(rhs).epsilon(1e-12));

  // loads: P1 reproduces constants, so the loads sum to the same total
  auto loads = target_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*target_space));
  Host2D ones("ones", target_space->GetLayout()->GetNumDofHolder(), components);
  Kokkos::deep_copy(ones, 1.0);
  SetDofs(loads, ones);
  auto distributed = source_space->CreateFunction<pcms::Real>(
    "", pcms::test::ComponentsBasis(*source_space));
  interpolator->ApplyTranspose(loads, distributed);
  const auto values = distributed.GetDOFHolderDataHost().GetValues();
  double total = 0;
  for (std::size_t i = 0;
       i < static_cast<std::size_t>(values.extent(0)) * values.extent(1); ++i) {
    total += values.data_handle()[i];
  }
  CHECK(total == Catch::Approx(ones.size()).epsilon(1e-12));
}

TEST_CASE("ApplyTranspose of a projection is its adjoint, for both solve "
          "methods",
          "[transfer][adjoint][projection]")
{
  Omega_h::Library lib;
  const int target_order = GENERATE(0, 1);
  const bool assembled = GENERATE(false, true);
  const bool lumped = GENERATE(false, true);
  CAPTURE(target_order, assembled, lumped);
  auto source_mesh = pcms::test::BuildSquareMesh(lib, 5);
  auto target_mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto options =
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH};
  const auto source_space =
    pcms::CreateLagrangeSpace(MakeGeometry(source_mesh), 1, 2, options);
  const auto target_space = pcms::CreateLagrangeSpace(MakeGeometry(target_mesh),
                                                      target_order, 2, options);
  pcms::ProjectionOptions projection_options;
  if (assembled) {
    projection_options.solve.method =
      pcms::LinearSolveOptions::Method::Assembled;
  }
  if (lumped) {
    projection_options.mass_matrix_type = pcms::MassMatrixType::Lumped;
  }
#ifndef PCMS_ENABLE_PETSC
  if (assembled && !lumped) {
    CHECK_THROWS(
      pcms::CreateProjection(source_space, target_space, projection_options));
    return;
  }
#endif
  const auto projection =
    pcms::CreateProjection(source_space, target_space, projection_options);
  const auto [lhs, rhs] = AdjointPair(*projection, source_space, target_space);
  CHECK(lhs == Catch::Approx(rhs).epsilon(1e-10));
  // with a control variate: B_VS^T z + I^T (v - B_VV z), z = M^-1 v
  const auto control_variate = pcms::CreateControlVariateProjection(
    source_space, target_space,
    {.mass_matrix_type = projection_options.mass_matrix_type});
  REQUIRE(control_variate->SupportsTranspose());
  const auto [cv_lhs, cv_rhs] =
    AdjointPair(*control_variate, source_space, target_space);
  CHECK(cv_lhs == Catch::Approx(cv_rhs).epsilon(1e-10));
}
