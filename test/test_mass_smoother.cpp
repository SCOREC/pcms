#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <petscmat.h>

#include <pcms/transfer/mass_smoother.hpp>
#include <pcms/transfer/omega_h_mass_integrator.hpp>
#include <pcms/transfer/petsc_utils.hpp>
#include <pcms/utility/arrays.h>
#include <pcms/utility/assert.h>
#include "field_test_utils.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

// The smoother only sees a BilinearFormIntegrator and a FieldLayout, so the
// first tests drive it with hand-built matrices behind a test integrator and a
// permuted global-id tag; the Omega_h mass integrator only enters afterwards.

namespace
{

using pcms::Real;
using Dense = std::vector<std::vector<double>>;

class DenseIntegrator : public pcms::BilinearFormIntegrator
{
public:
  explicit DenseIntegrator(const Dense& a, bool diagonal = false)
    : diagonal_(diagonal)
  {
    const auto n = static_cast<PetscInt>(a.size());
    pcms::createSeqAIJMat(PETSC_COMM_SELF, n, n, n, nullptr, &mat_);
    for (PetscInt r = 0; r < n; ++r) {
      for (PetscInt c = 0; c < n; ++c) {
        if (a[r][c] != 0.0) {
          MatSetValue(mat_, r, c, a[r][c], INSERT_VALUES);
        }
      }
    }
    MatAssemblyBegin(mat_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat_, MAT_FINAL_ASSEMBLY);
  }
  ~DenseIntegrator() override
  {
    if (mat_ != nullptr) {
      MatDestroy(&mat_);
    }
  }
  Mat GetMatrix() const noexcept override { return mat_; }
  bool IsDiagonal() const noexcept override { return diagonal_; }

private:
  Mat mat_ = nullptr;
  bool diagonal_;
};

// Rows of a 4-node 1D linear chain mass matrix with h = 1.
Dense ChainMass()
{
  Dense a = {{2, 1, 0, 0}, {1, 4, 1, 0}, {0, 1, 4, 1}, {0, 0, 1, 2}};
  for (auto& row : a) {
    for (auto& v : row) {
      v /= 6.0;
    }
  }
  return a;
}

Dense ExtractDense(Mat M, int n)
{
  Dense a(n, std::vector<double>(n, 0.0));
  for (PetscInt r = 0; r < n; ++r) {
    for (PetscInt c = 0; c < n; ++c) {
      PetscScalar v = 0.0;
      MatGetValues(M, 1, &r, 1, &c, &v);
      a[r][c] = PetscRealPart(v);
    }
  }
  return a;
}

std::vector<double> RowSums(const Dense& a)
{
  std::vector<double> m(a.size(), 0.0);
  for (std::size_t r = 0; r < a.size(); ++r) {
    for (double v : a[r]) {
      m[r] += v;
    }
  }
  return m;
}

// One sweep of inv(M_L) M in row order.
std::vector<double> SweepDense(const Dense& a, const std::vector<double>& u)
{
  const auto m = RowSums(a);
  std::vector<double> out(u.size(), 0.0);
  for (std::size_t r = 0; r < a.size(); ++r) {
    for (std::size_t c = 0; c < a.size(); ++c) {
      out[r] += a[r][c] * u[c];
    }
    out[r] /= m[r];
  }
  return out;
}

std::vector<double> HolderValues(const pcms::Field<Real>& field)
{
  const auto flat = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  return std::vector<double>(flat.data_handle(),
                             flat.data_handle() + flat.size());
}

void SetHolderValues(pcms::Field<Real>& field, const std::vector<double>& v,
                     int num_components = 1)
{
  const auto n = static_cast<pcms::LO>(v.size()) / num_components;
  field.SetDOFHolderDataHost(pcms::Rank2View<const Real, pcms::HostMemorySpace>(
    v.data(), n, num_components));
}

std::vector<double> RandomValues(std::size_t n, unsigned seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<double> v(n);
  for (auto& x : v) {
    x = dist(gen);
  }
  return v;
}

Omega_h::Mesh BuildReferenceTriangle(Omega_h::Library& lib)
{
  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(
    &mesh, OMEGA_H_SIMPLEX, 2, Omega_h::LOs({0, 1, 2}),
    Omega_h::Reals({0.0, 0.0, 1.0, 0.0, 0.0, 1.0}));
  pcms::test::AddDefaultClassification(mesh);
  return mesh;
}

// Lumped mass per DOF holder from the consistent Omega_h mass matrix.
std::vector<double> HolderLumpedMass(const pcms::FunctionSpace& space)
{
  auto integrator = pcms::BuildOmegaHMassIntegrator(space);
  const auto& layout = *space.GetLayout();
  const auto n = layout.GetNumOwnedDofHolder();
  Vec row_sums = nullptr;
  pcms::createSeqVec(PETSC_COMM_SELF, n, &row_sums);
  MatGetRowSum(integrator->GetMatrix(), row_sums);
  const PetscScalar* vals = nullptr;
  VecGetArrayRead(row_sums, &vals);
  const auto perm = layout.GetGlobalToLocalPermutationHost();
  std::vector<double> m(n);
  for (pcms::LO i = 0; i < n; ++i) {
    m[i] = PetscRealPart(vals[perm(i)]);
  }
  VecRestoreArrayRead(row_sums, &vals);
  VecDestroy(&row_sums);
  return m;
}

} // namespace

TEST_CASE("MassSmoother: works on any integrator and honours the permutation",
          "[mass_smoother]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  mesh.add_tag<Omega_h::GO>(Omega_h::VERT, "perm_global", 1,
                            Omega_h::GOs({3, 0, 2, 1}));
  const auto chain = ChainMass();
  const std::vector<double> v = {0.1, 0.7, 0.3, 0.9};

  auto expected_for = [&](const std::vector<double>& holder_values,
                          const pcms::FieldLayout& layout) {
    const auto perm = layout.GetGlobalToLocalPermutationHost();
    std::vector<double> u_row(4);
    for (int i = 0; i < 4; ++i) {
      u_row[perm(i)] = holder_values[i];
    }
    const auto s_row = SweepDense(chain, u_row);
    std::vector<double> expected(4);
    for (int i = 0; i < 4; ++i) {
      expected[i] = s_row[perm(i)];
    }
    return std::make_pair(expected, u_row);
  };

  SECTION("scalar field")
  {
    auto space = pcms::test::MakeP1Space(mesh, "perm_global");
    const auto perm = space->GetLayout()->GetGlobalToLocalPermutationHost();
    REQUIRE(perm(0) == 3);
    REQUIRE(perm(1) == 0);
    REQUIRE(perm(2) == 2);
    REQUIRE(perm(3) == 1);

    pcms::MassSmoother smoother(*space,
                                std::make_unique<DenseIntegrator>(chain));
    auto field = space->CreateFunction<Real>();
    SetHolderValues(field, v);
    const auto [expected, u_row] = expected_for(v, *space->GetLayout());

    smoother.Apply(field);
    const auto got = HolderValues(field);
    for (int i = 0; i < 4; ++i) {
      CAPTURE(i);
      CHECK(got[i] == Catch::Approx(expected[i]).margin(1e-14));
      CHECK(got[i] >= *std::min_element(v.begin(), v.end()) - 1e-14);
      CHECK(got[i] <= *std::max_element(v.begin(), v.end()) + 1e-14);
    }

    const auto m = RowSums(chain);
    double before = 0.0;
    double after = 0.0;
    for (int i = 0; i < 4; ++i) {
      before += m[perm(i)] * v[i];
      after += m[perm(i)] * got[i];
    }
    CHECK(after == Catch::Approx(before).epsilon(1e-14));
  }

  SECTION("each component is smoothed independently")
  {
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 2, pcms::CoordinateSystem::Cartesian, "perm_global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    pcms::MassSmoother smoother(*space,
                                std::make_unique<DenseIntegrator>(chain));
    const std::vector<double> w = {0.9, 0.3, 0.7, 0.1};
    std::vector<double> interleaved(8);
    for (int i = 0; i < 4; ++i) {
      interleaved[2 * i] = v[i];
      interleaved[2 * i + 1] = w[i];
    }
    auto field = space->CreateFunction<Real>();
    SetHolderValues(field, interleaved, 2);
    const auto expected_v = expected_for(v, *space->GetLayout()).first;
    const auto expected_w = expected_for(w, *space->GetLayout()).first;

    smoother.Apply(field);
    const auto got = HolderValues(field);
    for (int i = 0; i < 4; ++i) {
      CAPTURE(i);
      CHECK(got[2 * i] == Catch::Approx(expected_v[i]).margin(1e-14));
      CHECK(got[2 * i + 1] == Catch::Approx(expected_w[i]).margin(1e-14));
    }
  }
}

TEST_CASE("MassSmoother: rejects matrices that cannot smooth",
          "[mass_smoother]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP1Space(mesh);
  const Dense identity = {
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

  SECTION("diagonal matrix")
  {
    CHECK_THROWS_AS(
      pcms::MassSmoother(*space,
                         std::make_unique<DenseIntegrator>(identity, true)),
      pcms::pcms_error);
  }
  SECTION("non-positive row sum")
  {
    const Dense zero_row = {
      {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 0}};
    CHECK_THROWS_AS(
      pcms::MassSmoother(*space, std::make_unique<DenseIntegrator>(zero_row)),
      pcms::pcms_error);
  }
  SECTION("negative entry with positive row sums")
  {
    const Dense stiffness_like = {
      {2, -1, 0, 0}, {-1, 2, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    CHECK_THROWS_AS(
      pcms::MassSmoother(*space,
                         std::make_unique<DenseIntegrator>(stiffness_like)),
      pcms::pcms_error);
  }
  SECTION("size mismatch")
  {
    const Dense three = {{1, 1, 0}, {1, 1, 1}, {0, 1, 1}};
    CHECK_THROWS_AS(
      pcms::MassSmoother(*space, std::make_unique<DenseIntegrator>(three)),
      pcms::pcms_error);
  }
  SECTION("P0 space has a diagonal mass matrix")
  {
    auto p0 = pcms::test::MakeP0Space(mesh);
    CHECK_THROWS_AS(pcms::BuildMassSmoother(*p0), pcms::pcms_error);
  }
  SECTION("field from another space")
  {
    auto other = pcms::test::MakeP1Space(mesh);
    auto smoother = pcms::BuildMassSmoother(*space);
    auto field = other->CreateFunction<Real>();
    SetHolderValues(field, {1.0, 2.0, 3.0, 4.0});
    CHECK_THROWS_AS(smoother->Apply(field), pcms::pcms_error);
  }
}

TEST_CASE("MassSmoother: single element closed form", "[mass_smoother]")
{
  Omega_h::Library lib;

  SECTION("reference triangle: 1/2 on the vertex, 1/4 on the others")
  {
    auto mesh = BuildReferenceTriangle(lib);
    auto space = pcms::test::MakeP1Space(mesh);
    auto smoother = pcms::BuildMassSmoother(*space);
    auto field = space->CreateFunction<Real>();
    SetHolderValues(field, {1.0, 0.0, 0.0});
    smoother->Apply(field);
    const auto got = HolderValues(field);
    CHECK(got[0] == Catch::Approx(0.5).margin(1e-14));
    CHECK(got[1] == Catch::Approx(0.25).margin(1e-14));
    CHECK(got[2] == Catch::Approx(0.25).margin(1e-14));
  }

  SECTION("reference tet: 2/5 on the vertex, 1/5 on the others")
  {
    auto mesh = pcms::test::BuildReferenceTet(lib);
    auto space = pcms::test::MakeP1Space(mesh);
    auto smoother = pcms::BuildMassSmoother(*space);
    auto field = space->CreateFunction<Real>();
    SetHolderValues(field, {1.0, 0.0, 0.0, 0.0});
    smoother->Apply(field);
    const auto got = HolderValues(field);
    CHECK(got[0] == Catch::Approx(0.4).margin(1e-14));
    for (int i = 1; i < 4; ++i) {
      CAPTURE(i);
      CHECK(got[i] == Catch::Approx(0.2).margin(1e-14));
    }
  }
}

TEST_CASE("MassSmoother: conserves the integral and never widens the range",
          "[mass_smoother]")
{
  Omega_h::Library lib;
  auto run = [&](Omega_h::Mesh& mesh) {
    auto space = pcms::test::MakeP1Space(mesh);
    auto smoother = pcms::BuildMassSmoother(*space);
    auto field = space->CreateFunction<Real>();
    const auto n = static_cast<std::size_t>(mesh.nverts());
    SetHolderValues(field, RandomValues(n, 7));

    const double integral = pcms::test::IntegrateP1Field(mesh, field);
    auto range = pcms::test::P1FieldRange(field);
    for (int sweep = 0; sweep < 10; ++sweep) {
      smoother->Apply(field);
      CAPTURE(sweep);
      CHECK(pcms::test::IntegrateP1Field(mesh, field) ==
            Catch::Approx(integral).epsilon(1e-13));
      const auto new_range = pcms::test::P1FieldRange(field);
      CHECK(new_range.first >= range.first - 1e-14);
      CHECK(new_range.second <= range.second + 1e-14);
      range = new_range;
    }
    CHECK(range.first >= 0.0);

    SetHolderValues(field, std::vector<double>(n, 0.75));
    smoother->Apply(field);
    for (double x : HolderValues(field)) {
      CHECK(x == Catch::Approx(0.75).epsilon(1e-14));
    }
  };

  SECTION("2D box")
  {
    auto mesh =
      Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 5, 5, 0);
    run(mesh);
  }
  SECTION("3D box")
  {
    auto mesh = pcms::test::BuildUnitCube(lib, 3);
    run(mesh);
  }
}

TEST_CASE("MassSmoother: repeated sweeps match the dense operator and "
          "dissipate energy",
          "[mass_smoother]")
{
  Omega_h::Library lib;

  SECTION("three sweeps equal S^3 u")
  {
    auto mesh = pcms::test::BuildUnitSquare(lib, 0);
    auto space = pcms::test::MakeP1Space(mesh);
    auto integrator = pcms::BuildOmegaHMassIntegrator(*space);
    const auto dense = ExtractDense(integrator->GetMatrix(), mesh.nverts());
    const auto perm = space->GetLayout()->GetGlobalToLocalPermutationHost();

    const auto v = RandomValues(4, 11);
    std::vector<double> u_row(4);
    for (int i = 0; i < 4; ++i) {
      u_row[perm(i)] = v[i];
    }
    for (int k = 0; k < 3; ++k) {
      u_row = SweepDense(dense, u_row);
    }

    auto smoother = pcms::BuildMassSmoother(*space);
    auto field = space->CreateFunction<Real>();
    SetHolderValues(field, v);
    for (int k = 0; k < 3; ++k) {
      smoother->Apply(field);
    }
    const auto got = HolderValues(field);
    for (int i = 0; i < 4; ++i) {
      CAPTURE(i);
      CHECK(got[i] == Catch::Approx(u_row[perm(i)]).margin(1e-14));
    }
  }

  SECTION("lumped-mass energy is non-increasing")
  {
    auto mesh = pcms::test::BuildUnitCube(lib, 3);
    auto space = pcms::test::MakeP1Space(mesh);
    const auto m = HolderLumpedMass(*space);
    auto smoother = pcms::BuildMassSmoother(*space);
    auto field = space->CreateFunction<Real>();
    SetHolderValues(field, RandomValues(m.size(), 3));

    auto energy = [&]() {
      double e = 0.0;
      const auto u = HolderValues(field);
      for (std::size_t i = 0; i < u.size(); ++i) {
        e += m[i] * u[i] * u[i];
      }
      return e;
    };
    double previous = energy();
    for (int sweep = 0; sweep < 5; ++sweep) {
      smoother->Apply(field);
      const double current = energy();
      CAPTURE(sweep, previous, current);
      CHECK(current <= previous * (1.0 + 1e-14));
      CHECK(current < previous);
      previous = current;
    }
  }
}

TEST_CASE("MassSmoother: does not reproduce linear fields", "[mass_smoother]")
{
  // A vertex is mapped to the mass centroid of its basis function, which is
  // not the vertex on asymmetric stars; the integral is still conserved.
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP1Space(mesh);
  auto smoother = pcms::BuildMassSmoother(*space);
  auto field = space->CreateFunction<Real>();
  pcms::test::SetField(field, OMEGA_H_LAMBDA(Real x, Real) { return x; });
  const double integral = pcms::test::IntegrateP1Field(mesh, field);

  smoother->Apply(field);

  CHECK(pcms::test::IntegrateP1Field(mesh, field) ==
        Catch::Approx(integral).epsilon(1e-13));
  const auto coords = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  const auto got = HolderValues(field);
  double max_error = 0.0;
  for (int i = 0; i < mesh.nverts(); ++i) {
    max_error = std::max(max_error, std::fabs(got[i] - coords[2 * i]));
  }
  CHECK(max_error > 0.2);
}
