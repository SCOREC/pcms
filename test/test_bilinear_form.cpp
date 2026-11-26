#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/field/detail/field_ops.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/omega_h/lagrange_layout.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/transfer/forms.hpp>
#include <pcms/transfer/linear_solver.hpp>
#include <pcms/configuration.hpp>
#include "field_test_utils.hpp"
#include "forms_test_utils.hpp"
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include <pcms/transfer/detail/projection.hpp>

namespace pcms
{
using omega_h::CreateGeometry;
using omega_h::OmegaHLagrangeLayout;
} // namespace pcms

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace
{

using pcms::test::AssembleMass;
using pcms::test::Matrix;

// lumped = B 1
void Lump(const pcms::transfer::detail::BilinearForm& form,
          pcms::Field<pcms::Real>& lumped)
{
  auto ones = form.TestSpace()->CreateFunction<pcms::Real>();
  pcms::detail::Const(pcms::Real(1), ones);
  form.Apply(ones, lumped);
}

pcms::transfer::detail::BilinearForm MassForm(
  const std::shared_ptr<const pcms::FunctionSpace>& space)
{
  return pcms::transfer::detail::BilinearForm(space);
}

// The reference P1 mass matrix of a triangle mesh, summed over the elements:
// the element block is area / 12 * (1 + delta_ij)
Matrix BuildReferenceMassMap(Omega_h::Mesh& mesh,
                             const pcms::FunctionSpace& space)
{
  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space.GetLayout());
  REQUIRE(layout != nullptr);
  const auto areas =
    Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
  const auto faces2nodes =
    Omega_h::HostRead(mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b);
  const auto& gids = layout->GetGidsHost();
  Matrix ref;
  for (int e = 0; e < mesh.nelems(); ++e) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        ref[{static_cast<pcms::GO>(gids[faces2nodes[e * 3 + i]]),
             static_cast<pcms::GO>(gids[faces2nodes[e * 3 + j]])}] +=
          areas[e] / 12.0 * (i == j ? 2.0 : 1.0);
      }
    }
  }
  return ref;
}

} // namespace

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("BilinearForm: assembled matrix matches buildElementMassMatrix",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP1Space(mesh);
  const auto ref = BuildReferenceMassMap(mesh, *space);
  const auto form = MassForm(space);
  const auto matrix = AssembleMass(form);
  for (const auto& [key, ref_val] : ref) {
    CAPTURE(key.first, key.second, ref_val);
    REQUIRE(matrix.count(key) == 1);
    CHECK(matrix.at(key) == Catch::Approx(ref_val).epsilon(1e-12));
  }
}

TEST_CASE("BilinearForm: row sums are positive and total the domain area",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP1Space(mesh);
  const auto matrix = AssembleMass(MassForm(space));
  std::map<pcms::GO, pcms::Real> row_sums;
  for (const auto& [key, value] : matrix) {
    row_sums[key.first] += value;
  }
  pcms::Real grand_total = 0.0;
  for (const auto& [row, sum] : row_sums) {
    CAPTURE(row, sum);
    CHECK(sum > 0.0);
    grand_total += sum;
  }
  CHECK(grand_total == Catch::Approx(1.0).epsilon(1e-10));
}

TEST_CASE("BilinearForm: lumped mass equals the consistent row sums",
          "[bilinear_form]")
{
  // Row-sum lumping: M_L(i,i) = sum_j M_C(i,j) = integral(N_i dx). The trace
  // equals the domain area.
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP1Space(mesh);
  const auto form = MassForm(space);
  const auto matrix = AssembleMass(form);
  auto lumped = space->CreateFunction<pcms::Real>();
  Lump(form, lumped);
  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space->GetLayout());
  const auto& gids = layout->GetGidsHost();
  const auto values = lumped.GetDOFHolderDataHost().GetValues();
  pcms::Real trace = 0.0;
  for (int v = 0; v < mesh.nverts(); ++v) {
    pcms::Real row_sum = 0.0;
    for (const auto& [key, value] : matrix) {
      if (key.first == static_cast<pcms::GO>(gids[v])) {
        row_sum += value;
      }
    }
    CAPTURE(v, row_sum);
    CHECK(values(v, 0) == Catch::Approx(row_sum).epsilon(1e-12));
    trace += values(v, 0);
  }
  CHECK(trace == Catch::Approx(1.0).epsilon(1e-10));
}

TEST_CASE("BilinearForm: P0 mass is diagonal with the element areas, lumped "
          "or not",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto geometry = pcms::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  auto space = pcms::CreateLagrangeSpace(
    geometry, 0, 1,
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
  const auto form = MassForm(space);
  const auto matrix = AssembleMass(form);
  auto lumped = space->CreateFunction<pcms::Real>();
  Lump(form, lumped);
  const auto lumped_values = lumped.GetDOFHolderDataHost().GetValues();

  const auto elem_areas = Omega_h::measure_elements_real(&mesh);
  const auto elem_areas_h = Omega_h::HostRead<Omega_h::Real>(elem_areas);
  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space->GetLayout());
  const auto& gids = layout->GetGidsHost();
  for (const auto& [key, value] : matrix) {
    if (key.first != key.second) {
      CHECK(value == Catch::Approx(0.0).margin(1e-14));
    }
  }
  for (int e = 0; e < mesh.nelems(); ++e) {
    const pcms::GO gid = gids[e];
    CAPTURE(e, elem_areas_h[e]);
    CHECK(matrix.at({gid, gid}) ==
          Catch::Approx(elem_areas_h[e]).epsilon(1e-12));
    CHECK(lumped_values(e, 0) == Catch::Approx(elem_areas_h[e]).epsilon(1e-12));
  }
}

TEST_CASE("BilinearForm: a multi-component space has the same mass matrix, "
          "and the form needs a space with a natural quadrature",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::CreateLagrangeSpace(
    pcms::CreateGeometry(mesh, pcms::detail::DeferredCartesian()), 1, 2,
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
  const auto matrix = AssembleMass(MassForm(space));
  const auto scalar = AssembleMass(MassForm(pcms::test::MakeP1Space(mesh)));
  REQUIRE(matrix.size() == scalar.size());
  for (const auto& [key, value] : scalar) {
    CHECK(matrix.at(key) == Catch::Approx(value).epsilon(1e-12));
  }
}

TEST_CASE("BilinearForm: a trial on another mesh, applied, transposed and "
          "paired",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto coarse = pcms::test::BuildSquareMesh(lib, 2);
  auto fine = pcms::test::BuildSquareMesh(lib, 4);
  const auto test_space = pcms::test::MakeP1Space(coarse);
  const auto trial_space = pcms::test::MakeP1Space(fine);
  const auto linear = [](pcms::Real x, pcms::Real y, pcms::Real* out) {
    out[0] = 1 + x - 2 * y;
  };
  auto u = trial_space->CreateFunction<pcms::Real>();
  auto v = test_space->CreateFunction<pcms::Real>();
  auto w = test_space->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(u, linear);
  pcms::test::SetFieldComponents(w, linear);
  pcms::test::SetFieldComponents(
    v, [](pcms::Real x, pcms::Real y, pcms::Real* out) { out[0] = x + y; });
  // int (x + y) (1 + x - 2 y) over the unit square
  const pcms::Real exact = 5.0 / 12.0;
  const pcms::transfer::detail::BilinearForm form(test_space, trial_space);
  CHECK(form.Pair(v, u) == Catch::Approx(exact).epsilon(1e-12));
  auto functional = test_space->CreateFunction<pcms::Real>();
  form.Apply(u, functional);
  CHECK(pcms::detail::DotProd(v, functional) ==
        Catch::Approx(exact).epsilon(1e-12));
  auto transposed = trial_space->CreateFunction<pcms::Real>();
  form.ApplyTranspose(v, transposed);
  CHECK(pcms::detail::DotProd(transposed, u) ==
        Catch::Approx(exact).epsilon(1e-12));
  // the same points with the test space as the trial
  const auto square = form.WithTrial(test_space);
  CHECK(square.Pair(v, w) == Catch::Approx(exact).epsilon(1e-12));
  CHECK(square.IsSquare());
  CHECK_FALSE(form.IsSquare());
}

namespace
{

// Field data without the coefficient operations.
class NoOperationsData final : public pcms::FieldData<pcms::Real>
{
public:
  explicit NoOperationsData(pcms::ValueBasis basis) : basis_(std::move(basis))
  {
  }
  pcms::FieldValueType GetValueType() const override
  {
    return pcms::detail::ValueTypeOfRank(basis_.Rank());
  }
  const pcms::ValueBasis& GetValueBasis() const override { return basis_; }

private:
  pcms::ValueBasis basis_;
};

// A space on the layout and evaluator of another, whose fields' data has no
// coefficient operations.
class NoOperationsSpace final : public pcms::FunctionSpace
{
public:
  explicit NoOperationsSpace(std::shared_ptr<const pcms::FunctionSpace> space)
    : FunctionSpace(space->GetLayout(), space, space->GetGeometry()),
      space_(std::move(space))
  {
  }

protected:
  pcms::Quadrature CreateQuadrature() const override
  {
    return space_->GetQuadrature();
  }
  pcms::FieldDataVariant CreateFieldData(pcms::Type /*storage_type*/,
                                         pcms::ValueBasis basis) const override
  {
    return std::unique_ptr<pcms::FieldData<pcms::Real>>(
      std::make_unique<NoOperationsData>(std::move(basis)));
  }

private:
  std::shared_ptr<const pcms::FunctionSpace> space_;
};

} // namespace

TEST_CASE("CreateLinearSolver: a target without coefficient operations throws",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 2);
  const auto space =
    std::make_shared<const NoOperationsSpace>(pcms::test::MakeP1Space(mesh));
  const pcms::transfer::detail::BilinearForm form(space);
  CHECK_THROWS_WITH(
    pcms::transfer::detail::CreateLinearSolver(
      form, pcms::MassMatrixType::Consistent,
      pcms::detail::LinearSolveOptions{}),
    Catch::Matchers::ContainsSubstring("no coefficient operations"));
}

TEST_CASE("CreateLinearSolver: the assembled and matrix-free solves recover "
          "u, the lumped solve divides by B 1",
          "[bilinear_form]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  const auto space = pcms::test::MakeP1Space(mesh);
  const auto form = MassForm(space);
  auto u = space->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(
    u, [](pcms::Real x, pcms::Real y, pcms::Real* out) {
      out[0] = std::sin(3 * x) + y * y;
    });
  auto bu = space->CreateFunction<pcms::Real>();
  form.Apply(u, bu);
  const auto u_host = u.GetDOFHolderDataHost().GetValues();
  // the default is the assembled solve for a serial form with rows when PCMS
  // has PETSc, the matrix-free one otherwise
  const bool matrix_free_default =
    dynamic_cast<const pcms::transfer::detail::CgSolver*>(
      pcms::transfer::detail::CreateLinearSolver(
        form, pcms::MassMatrixType::Consistent, {})
        .get()) != nullptr;
#ifdef PCMS_ENABLE_PETSC
  CHECK_FALSE(matrix_free_default);
#else
  CHECK(matrix_free_default);
#endif
  for (const auto method :
       {pcms::detail::LinearSolveOptions::Method::MatrixFree,
#ifdef PCMS_ENABLE_PETSC
        pcms::detail::LinearSolveOptions::Method::Assembled,
#endif
        pcms::detail::LinearSolveOptions::Method::Automatic}) {
    auto x = space->CreateFunction<pcms::Real>();
    pcms::transfer::detail::CreateLinearSolver(
      form, pcms::MassMatrixType::Consistent,
      pcms::detail::LinearSolveOptions{.method = method})
      ->Solve(bu, x);
    const auto x_host = x.GetDOFHolderDataHost().GetValues();
    for (std::size_t h = 0; h < x_host.extent(0); ++h) {
      CHECK(x_host(h, 0) == Catch::Approx(u_host(h, 0)).margin(1e-8));
    }
  }
  auto lumped = space->CreateFunction<pcms::Real>();
  Lump(form, lumped);
  auto x = space->CreateFunction<pcms::Real>();
  pcms::transfer::detail::CreateLinearSolver(form, pcms::MassMatrixType::Lumped,
                                             pcms::detail::LinearSolveOptions{})
    ->Solve(bu, x);
  const auto x_host = x.GetDOFHolderDataHost().GetValues();
  const auto b_host = bu.GetDOFHolderDataHost().GetValues();
  const auto l_host = lumped.GetDOFHolderDataHost().GetValues();
  for (std::size_t h = 0; h < x_host.extent(0); ++h) {
    CHECK(x_host(h, 0) ==
          Catch::Approx(b_host(h, 0) / l_host(h, 0)).epsilon(1e-12));
  }
}

TEST_CASE("BilinearForm: an evaluate-only space pairs but does not apply",
          "[bilinear_form]")
{
  pcms::UniformGrid<1> grid;
  grid.bot_left = {0.0};
  grid.edge_length = {1.0};
  grid.divisions = {8};
  const auto geometry =
    pcms::uniform_grid::CreateGeometry(grid, pcms::csys::CreateCartesian(1));
  const auto spline = pcms::CreateSplineSpace(geometry);
  auto f = spline->CreateFunction<pcms::Real>();
  pcms::test::SetFieldComponents(
    f, [](pcms::Real x, pcms::Real* out) { out[0] = x * x; });
  const pcms::transfer::detail::BilinearForm form(
    spline, pcms::transfer::detail::FormOptions{
              .quadrature = pcms::CreateGaussQuadrature(geometry, 6)});
  CHECK(form.Pair(f, f) == Catch::Approx(0.2).epsilon(1e-12));
  CHECK_FALSE(form.TestEvaluator().SupportsTranspose());
  CHECK_FALSE(form.TrialEvaluator().SupportsTranspose());
  auto out = spline->CreateFunction<pcms::Real>();
  CHECK_THROWS_WITH(form.Apply(f, out),
                    Catch::Matchers::ContainsSubstring("not supported"));
  CHECK_THROWS_WITH(form.ApplyTranspose(f, out),
                    Catch::Matchers::ContainsSubstring("not supported"));
}
